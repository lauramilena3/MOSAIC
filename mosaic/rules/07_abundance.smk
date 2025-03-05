#ruleorder: mapReadsToContigsPE > mapReadsToContigsSE

rule subsampleReadsIllumina_PE_mapping:
	input:
		paired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}_read_count.txt",),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	params:
		n_reads=2000000
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/subsampleReadsIllumina_PE_mapping/{sample}_{sampling}.tsv"
	threads: 1
	# resources:
	# 	mem_mb=4000
	shell:
		"""
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads={params.n_reads} sampleseed=1
		"""

rule subsampleReadsIllumina_PE_mapping_7M:
	input:
		paired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}_read_count.txt",),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		forward_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/7M_{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/7M_{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	params:
		n_reads=7000000
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/subsampleReadsIllumina_PE_mapping/{sample}_{sampling}.tsv"
	threads: 1
	# resources:
	# 	mem_mb=4000
	shell:
		"""
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads={params.n_reads} sampleseed=1
		"""

rule subsampleReadsIllumina_PE_vOTU_mapping:
	input:
		viral_subsampling=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sub_sampling_reads.txt"),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
	output:
		forward_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.sub.fastq.gz"),
		reverse_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.sub.fastq.gz"),
		viral_subsampling=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sub_sampling_reads_final.txt"),
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/subsampleReadsIllumina_PE_vOTU_mapping/{sample}.tsv"
	threads: 1
	shell:
		"""
		reads=$(cat {input.viral_subsampling})
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads=$reads sampleseed=1
		echo $(( $(zgrep -Ec "$" {output.forward_paired}) / 4 )) > {output.viral_subsampling}
		"""

rule buildBowtieDB_assembly:
	input:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
	output:
		contigs_bt2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.1.bt2",
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_{sampling}_bowtie_assembly.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		bowtie2-build {input.scaffolds} {params.prefix} --threads {threads}
		"""

def input_mapping(wildcards):
	input_list = []
	if METAGENOME:
		input_list.append(dirs_dict["ASSEMBLY_TEST"] + f"/2M_{wildcards.sample}_forward_paired_clean.{wildcards.sampling}.fastq.gz")
		input_list.append(dirs_dict["ASSEMBLY_TEST"] + f"/2M_{wildcards.sample}_reverse_paired_clean.{wildcards.sampling}.fastq.gz")
	else:
		input_list.append(dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_forward_paired_clean.tot.fastq.gz")
		input_list.append(dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_reverse_paired_clean.tot.fastq.gz")
	return input_list
	
rule stat_mapReadsToAssembly:
	input:
		contigs_bt2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.1.bt2",
		forward_paired=lambda wildcards: input_mapping(wildcards)[0],
		reverse_paired=lambda wildcards: input_mapping(wildcards)[1],
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_{sample}_assembled_contigs.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_assembled_contigs.{sampling}.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs.{sampling}_covstats.txt",
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_{sampling}_assembly.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --all
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		"""

rule buildBowtieDB_viral:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		contigs_bt2=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}.1.bt2",
	params:
		prefix=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_{sampling}_bowtie_viral.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		bowtie2-build {input.positive_contigs} {params.prefix} --threads {threads}
		"""

rule stat_mapReadsToViral:
	input:
		contigs_bt2=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}.1.bt2",
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_{sample}_viral_contigs.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_viral_contigs.{sampling}.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs.{sampling}_covstats.txt",
	params:
		prefix=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_{sampling}_viral.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --all
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		"""

rule buildBowtieDB_derreplicated:
	input:
		derreplicated_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.fasta",
	output:
		contigs_bt2=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.1.bt2",
	params:
		prefix=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sampling}_bowtie_derreplicated.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		bowtie2-build {input.derreplicated_positive_contigs} {params.prefix} --threads {threads}
		"""

rule stat_mapReadsToDerreplicated:
	input:
		contigs_bt2=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.1.bt2",
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_{sample}_derreplicated_contigs.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_derreplicated_contigs.{sampling}.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs.{sampling}_covstats.txt",
	params:
		prefix=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_{sampling}_viral.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --all
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		"""


rule buildBowtieDB_unfiltered:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
	params:
		prefix=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sampling}_bowtie_unfiltered.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		bowtie2-build {input.representatives} {params.prefix} --threads {threads}
		"""

rule stat_mapReadsToUnfiltered:
	input:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_{sample}_unfiltered_contigs.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_unfiltered_contigs.{sampling}.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs.{sampling}_covstats.txt",
	params:
		prefix=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_{sampling}_unfiltered.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --all
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		"""

rule buildBowtieDB_filtered:
	input:
		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
	output:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
	params:
		prefix=dirs_dict["MAPPING_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sampling}_bowtie_filtered.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		bowtie2-build {input.filtered_representatives} {params.prefix} --threads {threads}
		"""

rule mapReadsToContigsPE:
	input:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.1.bt2",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}.sam"),
		bam=(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_{sample}.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_filtered_{sample}.{sampling}.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_unique_{sample}.{sampling}.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique_sorted.bam"),
		unique_sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique_sorted.bam.bai"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique_covstats.txt",
		basecov=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_basecov.txt",
		unique_basecov=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique_basecov.txt",
	params:
		prefix=dirs_dict["MAPPING_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_{sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --all
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		samtools view -@ 144 -hf 0x2 {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
		samtools view  -@ 144 -bS {output.unique_sam}> {output.unique_bam}
		samtools sort -@ 144 {output.unique_bam} -o {output.unique_sorted_bam}
		samtools index {output.unique_sorted_bam}
		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
		#genomecov
		bedtools genomecov -dz -ibam {output.filtered_bam} > {output.basecov}
		bedtools genomecov -dz -ibam {output.unique_sorted_bam}> {output.unique_basecov}
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
		"""

# rule mapReadsToContigsPE_sub:
# 	input:
# 		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.1.bt2",
# 		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_forward_paired_clean.tot.fastq.gz"),
# 		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_reverse_paired_clean.tot.fastq.gz"),
# 	output:
# 		sam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub.sam"),
# 		bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub.bam"),
# 		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_sorted.bam"),
# 		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_sorted.bam.bai"),
# 		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_filtered.bam"),
# 		flagstats=dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_{sample}_.sub.txt",
# 		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_filtered_{sample}.sub.txt",
# 		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_unique_{sample}.sub.txt",
# 		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_unique.sam"),
# 		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_unique.bam"),
# 		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_unique_sorted.bam"),
# 		covstats=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_covstats.txt",
# 		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_unique_covstats.txt",
# 		basecov=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_basecov.txt",
# 		unique_basecov=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_unique_basecov.txt",
# 	params:
# 		prefix=dirs_dict["MAPPING_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".sub",
# 	message:
# 		"Mapping reads to contigs"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_sub.tsv"
# 	threads: 8
# 	shell:
# 		"""
# 		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --all
# 		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
# 		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
# 		samtools index {output.sorted_bam}
# 		samtools flagstat {output.sorted_bam} > {output.flagstats}
# 		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
# 		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
# 		samtools view -@ 144 -hf 0x2 {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
# 		samtools view  -@ 144 -bS {output.unique_sam}> {output.unique_bam}
# 		samtools sort -@ 144 {output.unique_bam} -o {output.unique_sorted_bam}
# 		samtools index {output.unique_sorted_bam}
# 		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
# 		#genomecov
# 		bedtools genomecov -dz -ibam {output.filtered_bam} > {output.basecov}
# 		bedtools genomecov -dz -ibam {output.unique_sorted_bam}> {output.unique_basecov}
# 		#covstats
# 		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
# 		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
# 		"""

rule call_SNPs_sub:
	input:
		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		sorted_bam=(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_sorted.bam"),
	output:
		snp_temp=temp(dirs_dict["MAPPING_DIR"]+ "/{sample}_sub_SNP_calls.bcf"),
		snp=(dirs_dict["MAPPING_DIR"]+ "/{sample}_sub_SNP_calls.tsv"),
	message:
		"Calling SNPs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/SNP_calling/{sample}_sub.tsv"
	threads: 8
	shell:
		"""
		bcftools mpileup -Ou -I -d 1000 -f {input.filtered_representatives} {input.sorted_bam} | bcftools call -mv -Ob --ploidy 1 -o {output.snp_temp}
		bcftools view -i '%QUAL>=20' {output.snp_temp} | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" -o {output.snp}
		"""

rule buildBowtieDB_contaminants:
	input:
		contaminants=dirs_dict["CONTAMINANTS_DIR_POST"]+ "/{contaminant}.fasta",
	output:
		contigs_bt2=dirs_dict["CONTAMINANTS_DIR_POST"]+ "/{contaminant}.1.bt2",
	params:
		prefix=dirs_dict["CONTAMINANTS_DIR_POST"]+ "/{contaminant}",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{contaminant}_bowtie_contaminants.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 32
	shell:
		"""
		bowtie2-build {input.contaminants} {params.prefix} --threads {threads}
		"""


rule mapReads_contaminants:
	input:
		contigs_bt2=dirs_dict["CONTAMINANTS_DIR_POST"]+ "/{contaminant}.1.bt2",
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_reverse_paired_clean.tot.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_flagstats_{sample}_{contaminant}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_flagstats_filtered_{sample}.{contaminant}.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_flagstats_unique_{sample}.{contaminant}.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique_sorted.bam"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique_covstats.txt",
		basecov=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_basecov.txt",
		unique_basecov=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique_basecov.txt",
	params:
		prefix=dirs_dict["CONTAMINANTS_DIR_POST"]+ "/{contaminant}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_{contaminant}_contaminants.tsv"
	threads: 16
	shell:
		"""
		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --fast
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		samtools view -@ 144 -hf 0x2 {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
		samtools view  -@ 144 -bS {output.unique_sam}> {output.unique_bam}
		samtools sort -@ 144 {output.unique_bam} -o {output.unique_sorted_bam}
		samtools index {output.unique_sorted_bam}
		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
		#genomecov
		bedtools genomecov -dz -ibam {output.filtered_bam} > {output.basecov}
		bedtools genomecov -dz -ibam {output.unique_sorted_bam}> {output.unique_basecov}
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
		"""

rule buildBowtieDB_reference:
	input:
		contaminants=REFERENCE_DIR+ "/" + REFERENCE + ".fasta",
	output:
		contigs_bt2=REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2",
	params:
		prefix=REFERENCE_DIR+ "/" + REFERENCE + "",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/" + REFERENCE + "_bowtie_contaminants.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 32
	shell:
		"""
		bowtie2-build {input.contaminants} {params.prefix} --threads {threads}
		"""

rule buildBowtieDB_reference_long:
	input:
		contaminants=REFERENCE_DIR+ "/" + REFERENCE + ".fasta",
	output:
		contigs_bt2=REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2l",
	params:
		prefix=REFERENCE_DIR+ "/" + REFERENCE + "",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/" + REFERENCE + "_bowtie_contaminants.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 32
	shell:
		"""
		bowtie2-build {input.contaminants} {params.prefix} --threads {threads} --large-index
		"""

def input_bowtie_reference(wildcards):
	if LONG_INDEX:
		return(REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2l")
	else:
		return(REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2")

rule mapReads_reference:
	input:
		contigs_bt2=input_bowtie_reference,
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot.sam"),
		bam=(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_" + REFERENCE + "_{sample}.tot.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_filtered_" + REFERENCE + "_{sample}.tot.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_unique_" + REFERENCE + "_{sample}.tot.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique_sorted.bam"),
		unique_sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique_sorted.bam.bai"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique_covstats.txt",
		basecov=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_basecov.txt",
		unique_basecov=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique_basecov.txt",
	params:
		prefix=REFERENCE_DIR+ "/" + REFERENCE + "",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_tot_" + REFERENCE + "_contaminants.tsv"
	threads: 16
	shell:
		"""
		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --fast
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		samtools view -@ 144 -hf 0x2 {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
		samtools view  -@ 144 -bS {output.unique_sam}> {output.unique_bam}
		samtools sort -@ 144 {output.unique_bam} -o {output.unique_sorted_bam}
		samtools index {output.unique_sorted_bam}
		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
		#genomecov
		bedtools genomecov -dz -ibam {output.filtered_bam} > {output.basecov}
		bedtools genomecov -dz -ibam {output.unique_sorted_bam}> {output.unique_basecov}
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
		"""

rule mapReads_reference_sub:
	input:
		contigs_bt2=REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2",
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_reverse_paired_clean.tot.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub.sam"),
		bam=(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_" + REFERENCE + "_{sample}.sub.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_filtered_" + REFERENCE + "_{sample}.sub.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_unique_" + REFERENCE + "_{sample}.sub.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique_sorted.bam"),
		unique_sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique_sorted.bam.bai"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique_covstats.txt",
		basecov=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_basecov.txt",
		unique_basecov=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique_basecov.txt",
	params:
		prefix=REFERENCE_DIR+ "/" + REFERENCE + "",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_sub_" + REFERENCE + "_contaminants.tsv"
	threads: 16
	shell:
		"""
		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --fast
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		samtools view -@ 144 -hf 0x2 {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
		samtools view  -@ 144 -bS {output.unique_sam}> {output.unique_bam}
		samtools sort -@ 144 {output.unique_bam} -o {output.unique_sorted_bam}
		samtools index {output.unique_sorted_bam}
		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
		#genomecov
		bedtools genomecov -dz -ibam {output.filtered_bam} > {output.basecov}
		bedtools genomecov -dz -ibam {output.unique_sorted_bam}> {output.unique_basecov}
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
		"""

rule extract_mapped_reads:
	input:
		bam=(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}.bam"),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
	output:
		mapped_reads=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_mapped_reads.txt"),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_mapped.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_mapped.fastq.gz"),
	message:
		"Extracting mapped reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_tot_extract_mapped.tsv"
	threads: 16
	shell:
		"""
 		samtools view -F 4 {input.bam} | cut -f1 | sort | uniq > {output.mapped_reads}
		seqtk subseq {input.forward_paired} {output.mapped_reads} | gzip > {output.forward_paired}
		seqtk subseq {input.reverse_paired} {output.mapped_reads} | gzip > {output.reverse_paired}
		"""

# rule get_norm_RPKM:
# 	input:
# 		covstats_all=dirs_dict["MAPPING_DIR"]+ "/bbmap_covstats_{sample}.{sampling}_all.txt",
# 		covstats_toss=dirs_dict["MAPPING_DIR"]+ "/bbmap_covstats_{sample}.{sampling}_toss.txt",
# 	output:
# 		rpkm=dirs_dict["MAPPING_DIR"]+ "/norm_RPKM_{sample}_{sampling}.txt",
# 	params:
# 		ambiguous=config['ambiguous_mapping'],
# 	message:
# 		"Calculating normalised RPKM"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/normalise_RPKM/{sample}_{sampling}.tsv"
# 	threads: 1
# 	shell:
# 		"""
# 		perl ./scripts/Make-vOTU-RPKM-Norm.pl {input.covstats_toss} {input.covstats_all} {output.rpkm}
# 		"""

# rule plotWeesam:
# 	input:
# 		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/bbmap_{sample}_sorted.{sampling}.bam",
# 		weesam_dir=(config['weesam_dir']),
# 	output:
# 		html=directory(dirs_dict["MAPPING_DIR"]+ "/{sample}_weeSAM_{sampling}_html_results"),
# 		txt=dirs_dict["MAPPING_DIR"]+ "/{sample}_weeSAM_{sampling}.txt",
# 		dir=temp(directory("{sample}_weeSAM.{sampling}")),
# 	params:
# 		html=dirs_dict["MAPPING_DIR"]+ "/{sample}_weeSAM_{sampling}",
# 	message:
# 		"Plotting read mapping with weeSAM"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env5.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/weeSAM/{sample}_{sampling}.tsv"
# 	threads: 1
# 	shell:
# 		"""
# 		mkdir -p {output.dir}
# 		cd {output.dir}
# 		../{input.weesam_dir}/weeSAM --bam {input.bam_sorted} --html {params.html} --out {output.txt}
# 		"""

# rule getBreadthCoverage:
# 	input:
# 		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/bbmap_{sample}_sorted.{sampling}.bam",
# 	output:
# 		bam_cov=dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_genomecov.{sampling}.txt",
# 		cov_final=dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_coverage.{sampling}.txt",
# 		tmp_sort=temp(directory(dirs_dict["MAPPING_DIR"]+ "/temp_{sample}_{sampling}")),
# 	message:
# 		"Calculating breadth coverage contigs"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/getBreadthCoverage/{sample}_{sampling}.tsv"
# 	threads: 1
# 	priority: 1
# 	shell:
# 		"""
# 		bedtools genomecov -dz -ibam {input.bam_sorted} > {output.bam_cov}
# 		mkdir {output.tmp_sort}
# 		cut -f 1 {output.bam_cov} | sort -T {output.tmp_sort} | uniq -c | sort -nr -T {output.tmp_sort} | sed -e 's/^[[:space:]]*//' > {output.cov_final}
# 		"""

# rule tabletoBIOM:
# 	input:
# 		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_DB.{sampling}.txt",
# 	output:
# 		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_json.{sampling}.biom",
# 	message:
# 		"Getting vOTU tables"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/tabletoBIOM/{sampling}.tsv"
# 	threads: 1
# 	shell:
# 		"""
# 		biom convert -i {input.abundances} -o {output.abundances} --table-type="OTU table" --to-json
# 		"""

# rule getSummaryTable:
# 	input:
# #		hmm_results=dirs_dict["VIRAL_DIR"]+ "/hmm_parsed.{sampling}.out",
# #		table=dirs_dict["VIRAL_DIR"]+ "/viral_table.{sampling}.csv",
# #		genome_file=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/genome_by_genome_overview.csv",
# 	output:
# 		summary=dirs_dict["MAPPING_DIR"]+ "/vOTU_summary.{sampling}.txt",
# 	message:
# 		"Getting vOTU tables"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/getSummaryTable/{sampling}.tsv"
# 	threads: 1
# 	shell:
# 		"""
# 		touch {output.summary}
# 		"""
# rule createContigBBDb:
# 	input:
# 		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
# 	output:
# 		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
# 		contigs_info=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta.fai",
# 		contigs_lenght=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_lenght.{sampling}.txt",
# 	params:
# 		prefix=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
# 	message:
# 		"Creating contig DB with Bowtie2"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		bowtie2-build -f {input.representatives} {params.prefix}
# 		#Get genome file
# 		samtools faidx {input.representatives}
# 		bbmap.sh ref={input.representatives}
# 		awk -F' ' '{{print $1"	"$2}}' {output.contigs_info} > {output.contigs_lenght}
#
# 		"""

# rule getAbundancesPE_user:
# 	input:
# 		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
# 		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
# 		unpaired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{{sampling}}.txt", sample=SAMPLES),
# 		paired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.{{sampling}}.txt", sample=SAMPLES),
# 	output:
# 		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table.{sampling}.txt",
# 	message:
# 		"Getting vOTU tables"
# 	threads: 1
# 	run:
# 		import pandas as pd
# 		import numpy as np
#
# 		lenght=7000
# 		percentage=0.7
# 		min_bases=5000
# 		df_tpmean=pd.DataFrame()
# 		sampling=wildcards.sampling
# 		for sample in SAMPLES:
# 			#READ NUMBER
# 			paired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_paired_clean."+sampling+".txt")
# 			unpaired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_unpaired_clean."+sampling+".txt")
# 			paired=int(paired_size.readline())
# 			unpaired=int(unpaired_size.readline())
# 			reads=((paired*2)+unpaired)/1000000
# 			#NORMALIZE TP MEAN
# 			tpmean_file=dirs_dict["MAPPING_DIR"]+ "/" +sample+"_tpmean." + sampling + ".tsv"
# 			tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
# 			tpmean[sample] = tpmean[sample].apply(lambda x: x/reads)
# 			#REMOVE LOW COVERED CONTIGS
# 			breadth_file = dirs_dict["MAPPING_DIR"]+ "/" +sample+"_filtered_coverage." + sampling + ".txt"
# 			breadth = pd.read_csv(breadth_file, sep=" ", header=0, names=("breadth", "contig"))
# 			df=pd.merge(tpmean, breadth, on='contig', how='outer')
# 			#Divide dataframe in lenghts
# 			df['percentage']=df['breadth']/df['length']
# 			df=df.fillna(0)
# 			positive = df[(df['breadth']>7000) | (df['percentage']>percentage) ]
# 			if df_tpmean.empty:
# 				positive.drop("breadth", axis=1, inplace=True)
# 				positive.drop("length", axis=1, inplace=True)
# 				#positive.drop("percentage", axis=1, inplace=True)
# 				df_tpmean=positive
# 			else:
# 				positive.drop("length", axis=1, inplace=True)
# 				positive.drop("breadth", axis=1, inplace=True)
# 				#positive.drop("percentage", axis=1, inplace=True)
# 				df_tpmean=pd.merge(positive, df_tpmean, on='contig', how='outer')
# 		filename="vOTU_abundance_table." + sampling + ".txt"
# 		df_tpmean=df_tpmean.fillna(0)
# 		df_tpmean.rename(columns={'contig':'#OTU ID'}, inplace=True)
# 		df_tpmean.to_csv(dirs_dict["MAPPING_DIR"]+ "/" + filename, sep='\t', index=False, header=True)

# rule filterBAM:
# 	input:
# 		bam=dirs_dict["MAPPING_DIR"]+ "/bowtie_{sample}.{sampling}.bam",
# 	output:
# 		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam",
# 		bam_filtered=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted_filtered.{sampling}.bam",
# 	params:
# 		out_dir=dirs_dict["MAPPING_DIR"],
# 		temp_bam_filtered=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}_filtered.bam",
# 		p_ident=config['p_ident'],
# 	message:
# 		"Filtering reads in Bam file with BamM"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env2.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		samtools sort {input.bam} -o {output.bam_sorted}
# 		bamm filter --bamfile {output.bam_sorted} --percentage_id {params.p_ident} -o {params.out_dir}
# 		mv {params.temp_bam_filtered} {output.bam_filtered}
# 		"""

# rule getAbundancesPE:
# 	input:
# 		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_coverage.{{sampling}}.txt", sample=SAMPLES),
# 		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
# 		unpaired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{{sampling}}.txt", sample=SAMPLES),
# 		paired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.{{sampling}}.txt", sample=SAMPLES),
# 	output:
# 		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table.{sampling}.txt",
# 	message:
# 		"Getting vOTU tables"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/getAbundancesPE/{sampling}.tsv"
# 	threads: 1
# 	run:
# 		import pandas as pd
# 		import numpy as np
#
# 		lenght=7000
# 		percentage=0.7
# 		min_bases=5000
# 		df_tpmean=pd.DataFrame()
# 		sampling=wildcards.sampling
# 		for sample in SAMPLES:
# 			#READ NUMBER
# 			paired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_paired_clean."+sampling+".txt")
# 			unpaired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_unpaired_clean."+sampling+".txt")
# 			paired=int(paired_size.readline())
# 			unpaired=int(unpaired_size.readline())
# 			reads=((paired))/1000000
# 			#NORMALIZE TP MEAN
# 			tpmean_file=dirs_dict["MAPPING_DIR"]+ "/BamM_" +sample+"_tpmean." + sampling + ".tsv"
# 			tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
# 			tpmean[sample] = tpmean[sample].apply(lambda x: x/reads)
# 			#REMOVE LOW COVERED CONTIGS
# 			breadth_file = dirs_dict["MAPPING_DIR"]+ "/bedtools_" +sample+"_filtered_coverage." + sampling + ".txt"
# 			breadth = pd.read_csv(breadth_file, sep=" ", header=0, names=("breadth", "contig"))
# 			df=pd.merge(tpmean, breadth, on='contig', how='outer')
# 			#Divide dataframe in lenghts
# 			df['percentage' ]=df['breadth']/df['length']
# 			df=df.fillna(0)
# 			positive = df[(df['breadth']>7000) | (df['percentage']>percentage) ]
# 			if df_tpmean.empty:
# 				positive.drop("breadth", axis=1, inplace=True)
# 				positive.drop("length", axis=1, inplace=True)
# 				#positive.drop("percentage", axis=1, inplace=True)
# 				df_tpmean=positive
# 			else:
# 				positive.drop("length", axis=1, inplace=True)
# 				positive.drop("breadth", axis=1, inplace=True)
# 				#positive.drop("percentage", axis=1, inplace=True)
# 				df_tpmean=pd.merge(positive, df_tpmean, on='contig', how='outer')
# 		filename="vOTU_abundance_table." + sampling + ".txt"
# 		df_tpmean=df_tpmean.fillna(0)
# 		df_tpmean.rename(columns={'contig':'#OTU ID'}, inplace=True)
# 		df_tpmean.to_csv(dirs_dict["MAPPING_DIR"]+ "/" + filename, sep='\t', index=False, header=True)

# rule getAbundancesSE:
# 	input:
# 		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
# 		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
# 		unpaired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{{sampling}}.txt", sample=SAMPLES),
# 	output:
# 		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table.{sampling}.txt",
# 	message:
# 		"Getting vOTU tables"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/getAbundancesSE/{sampling}.tsv"
# 	threads: 1
# 	run:
# 		import pandas as pd
# 		import numpy as np
#
# 		lenght=7000
# 		percentage=0.7
# 		min_bases=5000
# 		SAMPLING_TYPE={sampl}
# 		for sampling in SAMPLING_TYPE:
# 			df_tpmean=pd.DataFrame()
# 			for sample in SAMPLES:
# 				#READ NUMBER
# 				unpaired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_unpaired_clean."+sampling+".txt")
# 				unpaired=int(unpaired_size.readline())
# 				reads=(unpaired)/1000000
# 				#NORMALIZE TP MEAN
# 				tpmean_file=dirs_dict["MAPPING_DIR"]+ "/" +sample+"_tpmean." + sampling + ".tsv"
# 				tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
# 				tpmean[sample] = tpmean[sample].apply(lambda x: x/reads)
# 				#REMOVE LOW COVERED CONTIGS
# 				breadth_file = dirs_dict["MAPPING_DIR"]+ "/" +sample+"_filtered_coverage." + sampling + ".txt"
# 				breadth = pd.read_csv(breadth_file, sep=" ", header=0, names=("breadth", "contig"))
# 				df=pd.merge(tpmean, breadth, on='contig', how='outer')
# 				#Divide dataframe in lenghts
# 				df['percentage']=df['breadth']/df['length']
# 				df=df.fillna(0)
# 				positive = df[(df['breadth']>7000) | (df['percentage']>percentage) ]
# 				if df_tpmean.empty:
# 					positive.drop("breadth", axis=1, inplace=True)
# 					positive.drop("length", axis=1, inplace=True)
# 					positive.drop("percentage", axis=1, inplace=True)
# 					df_tpmean=positive
# 				else:
# 					positive.drop("length", axis=1, inplace=True)
# 					positive.drop("breadth", axis=1, inplace=True)
# 					positive.drop("percentage", axis=1, inplace=True)
# 					df_tpmean=pd.merge(positive, df_tpmean, on='contig', how='outer')
# 			filename="vOTU_abundance_table." + sampling + ".txt"
# 			df_tpmean=df_tpmean.fillna(0)
# 			df_tpmean.rename(columns={'contig':'#OTU ID'}, inplace=True)
# 			df_tpmean.to_csv(dirs_dict["MAPPING_DIR"]+ "/" + filename, sep='\t', index=False, header=True)

# rule getAbundancesDB:
# 	input:
# 		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_coverage.{{sampling}}.txt", sample=SAMPLES),
# 		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
# 		counts=expand(dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_counts.{{sampling}}.tsv", sample=SAMPLES),
# 		paired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.{{sampling}}.txt", sample=SAMPLES),
# 	output:
# 		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_DB.{sampling}.txt",
# 		parsed_abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_DB_70.{sampling}.txt",
# 	message:
# 		"Getting vOTU tables"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/getAbundancesDB/{sampling}.tsv"
# 	threads: 1
# 	priority: 1
# 	run:
# 		import pandas as pd
# 		import numpy as np
# 		import os
# 		os.chdir(RESULTS_DIR)
# 		sampling=wildcards.sampling
# 		df_tpmean=pd.DataFrame()
# 		df_counts=pd.DataFrame()
#
# 		for sample in SAMPLES:
# 			#READ NUMBER
# 			paired_size=open("02_CLEAN_DATA"+ "/" +sample+"_paired_clean."+sampling+".txt")
# 			unpaired_size=open("02_CLEAN_DATA"+ "/" +sample+"_unpaired_clean."+sampling+".txt")
# 			paired=int(paired_size.readline())
# 			unpaired=int(unpaired_size.readline())
# 			#reads=((paired*2)+unpaired)/1000000
# 			#NORMALIZE TP MEAN
# 			tpmean_file="06_MAPPING"+ "/BamM_" +sample+"_tpmean." + sampling + ".tsv"
# 			tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
# 			tpmean[sample] = tpmean[sample].apply(lambda x: x*1000000/paired)
# 			tpmean["contig"] = tpmean["contig"].str.strip()
# 			#READ COUNTS
# 			print(paired)
# 			counts_file="06_MAPPING"+ "/BamM_" +sample+"_counts." + sampling + ".tsv"
# 			counts = pd.read_csv(counts_file, sep="\t", header=0, names=("contig", "length", sample))
# 			counts["contig"] = counts["contig"].str.strip()
#
# 			breadth_file = "06_MAPPING"+ "/bedtools_" +sample+"_coverage." + sampling + ".txt"
# 			#breadth = pd.read_csv(breadth_file, sep=" ", header=0, names=("breadth", "contig"))
# 			contig=[]
# 			brdth=[]
# 			with open(breadth_file) as fp:
# 				for line in fp:
# 					first_line = line.strip()
# 					if first_line=="":
# 						contig.append("")
# 						brdth.append("")
# 					else:
# 						brdth.append(first_line.split(" ", 1)[0])
# 						contig.append(first_line.split(" ", 1)[1])
# 				breadth = pd.DataFrame({'contig': contig,'breadth': brdth})
# 			# TPMEAN
# 			df=pd.merge(tpmean,breadth,left_on='contig',right_on='contig')
#
# 			df["breadth"] = pd.to_numeric(df["breadth"])
# 			df['percentage' ]=df['breadth']/df['length']
# 			df=df.fillna(0)
# 			df.drop("breadth", axis=1, inplace=True)
# 			df.drop("length", axis=1, inplace=True)
# 			df.columns = ['contig', sample + "_depth", sample + "_breadth" ]
# 			#REMOVE LOW COVERED CONTIGS
# 			df=df[df[sample + "_breadth"]>0]
# 			if df_tpmean.empty:
# 				df_tpmean=df
# 			else:
# 				df_tpmean=pd.merge(df, df_tpmean, on='contig', how='outer')
#
# 			#COUNTS
# 			df2=pd.merge(counts,breadth,left_on='contig',right_on='contig')
#
# 			df2["breadth"] = pd.to_numeric(df2["breadth"])
# 			df2['percentage']=df2['breadth']/df2['length']
# 			df2=df2.fillna(0)
# 			df2.drop("breadth", axis=1, inplace=True)
# 			df2.drop("length", axis=1, inplace=True)
# 			df2.columns = ['contig', sample + "_depth", sample + "_breadth" ]
# 			#REMOVE NON COVERED CONTIGS
# 			df2=df2[df2[sample + "_breadth"]>0]
# 			if df_counts.empty:
# 				df_counts=df2
# 			else:
# 				df_counts=pd.merge(df2, df_counts, on='contig', how='outer')
# 		#TP MEAN
# 		df_tpmean=df_tpmean.fillna(0)
# 		df_tpmean.rename(columns={'contig':'OTU'}, inplace=True)
#
# 		a_series = (df_tpmean != 0).any(axis=1)
# 		df_tpmean = df_tpmean.loc[a_series]
# 		df_tpmean.set_index('OTU', inplace=True)
#
# 		df_tpmean_70=df_tpmean.loc[(df_tpmean >= 0.7).any(axis=1)]
# 		cols_d = [c for c in df_tpmean_70.columns if c.lower()[-5:] == 'depth']
# 		cols_b = [c for c in df_tpmean_70.columns if c.lower()[-5:] != 'depth']
#
# 		df_tpmean_70_d=df_tpmean_70[cols_d]
# 		df_tpmean_70_b=df_tpmean_70[cols_b]
#
# 		filter_df=(df_tpmean_70_b >= 0.7)
# 		filter_df.columns=df_tpmean_70_d.columns
# 		filtered_df_tpmean_70_d=df_tpmean_70_d[filter_df].fillna(0)
#
# 		#COUNTS
# 		df_counts=df_counts.fillna(0)
# 		df_counts.rename(columns={'contig':'OTU'}, inplace=True)
#
# 		a_series = (df_counts != 0).any(axis=1)
# 		df_counts = df_counts.loc[a_series]
# 		df_counts.set_index('OTU', inplace=True)
#
# 		df_counts_70=df_counts.loc[(df_counts >= 0.7).any(axis=1)]
# 		cols_d = [c for c in df_counts_70.columns if c.lower()[-5:] == 'depth']
# 		cols_b = [c for c in df_counts_70.columns if c.lower()[-5:] != 'depth']
#
# 		df_counts_70_d=df_counts_70[cols_d]
# 		df_counts_70_b=df_counts_70[cols_b]
#
# 		filter_df2=(df_counts_70_b >= 0.7)
# 		filter_df2.columns=df_counts_70_d.columns
# 		filtered_df_counts_70_d=df_counts_70_d[filter_df2].fillna(0)
#
# 		#SAVE FILES
#
# 		filename="06_MAPPING/vOTU_abundance_table_DB." + sampling + ".txt"
# 		filename_70="06_MAPPING/vOTU_abundance_table_DB_70." + sampling + ".txt"
# 		filename_70_counts="06_MAPPING/vOTU_abundance_counts_DB_70." + sampling + ".txt"
#
# 		df_tpmean.to_csv(filename, sep='\t', header=True)
# 		filtered_df_tpmean_70_d.to_csv(filename_70, sep='\t', header=True)
# 		filtered_df_counts_70_d.to_csv(filename_70_counts, sep='\t', header=True)



# rule mapReadsToContigsSE:
# 	input:
# 		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
# 	output:
# 		sam=dirs_dict["MAPPING_DIR"]+ "/bbmap_{sample}.{sampling}.sam",
# 		bam=dirs_dict["MAPPING_DIR"]+ "/bbmap_{sample}.{sampling}.bam",
# 		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam",
# 		bam_indexed=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam.bai",
# 	params:
# 		contigs=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
# 	message:
# 		"Mapping reads to contigs"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsSE/{sample}_{sampling}.tsv"
# 	threads: 4
# 	shell:
# 		"""
# 		#Mapping reads to contigs
# 		#bowtie2 --non-deterministic -x {params.contigs} -U {input.unpaired} -S {output.sam} -p {threads}
# 		#Sam to Bam
# 		samtools view -b -S {output.sam} > {output.bam}
# 		samtools sort {output.bam} -o {output.bam_sorted}
# 		samtools index {output.bam_sorted}
# 		"""

# rule tpmeanPerConfidence:
# 	input:
# 		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam",
# 		bam_indexed=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam.bai",
# 	output:
# 		tpmean=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_tpmean.{sampling}.tsv",
# 		counts=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_counts.{sampling}.tsv",
# 	message:
# 		"Calculating tpmean depth coverage"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env2.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/createContigBowtieDb/{sample}_{sampling}.tsv"
# 	threads: 1
# 	shell:
# 		"""
# 		bamm parse -c {output.tpmean} -m tpmean -b {input.bam_sorted}
# 		bamm parse -c {output.counts} -m counts -b {input.bam_sorted}
# 		"""

# rule createContigBowtieDb:
# 	input:
# 		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
# 	output:
# 		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
# 		contigs_info=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta.fai",
# 		contigs_lenght=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_lenght.{sampling}.txt",
# 	params:
# 		prefix=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
# 	message:
# 		"Creating contig DB with Bowtie2"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/createContigBowtieDb/{sampling}.tsv"
# 	threads: 1
# 	shell:
# 		"""
# 		bowtie2-build -f {input.representatives} {params.prefix}
# 		#Get genome file
# 		samtools faidx {input.representatives}
# 		awk -F' ' '{{print $1"	"$2}}' {output.contigs_info} > {output.contigs_lenght}
# 		"""
