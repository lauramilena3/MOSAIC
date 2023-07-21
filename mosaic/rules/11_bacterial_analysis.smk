rule estimateBacterialGenomeCompletness:
	input:
		corrected2_racon=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
		checkm_db=(config['checkm_db']),
	output:
		checkMoutdir_temp=temp(directory(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}_temp")),
		checkMoutdir=directory(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}"),
	params:
		checkv_db=dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{sampling}",
	log:
		checkMoutdir=(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}_log"),
	message:
		"Estimating genome completeness with CheckM "
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/estimateGenomeCompletness/{sample}_{sampling}_checkm.tsv"
	threads: 4
	shell:
		"""
		mkdir {output.checkMoutdir_temp}
		cp {input.corrected2_racon} {output.checkMoutdir_temp}
		cd {output.checkMoutdir_temp}
		checkm lineage_wf -t {threads} -x fasta {output.checkMoutdir_temp} {output.checkMoutdir} 1> {log}
		"""

def input_vOTU_clustering(wildcards):
	input_list=[]
	input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.tot.fasta", sample=SAMPLES)),
	if len(config['additional_reference_contigs'])>0:
		input_list.append(config['additional_reference_contigs'])
	return input_list

rule merge_assembly:
	input:
		assembled_contigs=input_vOTU_clustering,
	output:
		merged_assemblies=dirs_dict["ASSEMBLY_DIR"] + "/merged_microbial_assembly.tot.fasta",
	message:
		"Merging microbial assemblies"
	shell:
		"""
		cat {input.assembled_contigs} > {output.merged_assemblies}
		"""

rule buildBowtieDB_microbial:
	input:
		merged_assemblies=dirs_dict["ASSEMBLY_DIR"] + "/merged_microbial_assembly.tot.fasta",
	output:
		contigs_bt2=dirs_dict["ASSEMBLY_DIR"] + "/merged_microbial_assembly.tot.1.bt2",
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/merged_microbial_assembly.tot",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sampling}_bowtie_microbial.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		bowtie2-build {input.merged_assemblies} {params.prefix} --threads {threads}
		"""

rule mapReadsToContigsPE:
	input:
		contigs_bt2=dirs_dict["ASSEMBLY_DIR"] + "/merged_microbial_assembly.{sampling}.1.bt2",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}.bam"),
		sorted_bam=(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_flagstats_{sample}_.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_flagstats_filtered_{sample}.{sampling}.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_flagstats_unique_{sample}.{sampling}.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique_sorted.bam"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique_covstats.txt",
		basecov=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_basecov.txt",
		unique_basecov=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique_basecov.txt",
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/merged_microbial_assembly.tot",
	message:
		"Mapping microbial reads to assembly"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_{sampling}_microbial.tsv"
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

rule bacterial_binning:
	input:
		merged_assemblies=dirs_dict["ASSEMBLY_DIR"] + "/merged_microbial_assembly.tot.fasta",
		sorted_bam=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam", sample=SAMPLES),
	output:
		metabat_outdir=directory(dirs_dict["MAPPING_DIR"] + "/MetaBAT_results/"),
	message:
		"Binning microbial contigs with MetaBAT"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/MetaBAT/binning.tsv"
	threads: 8
	shell:
		"""
		mkdir {output.metabat_outdir}
		cd {output.metabat_outdir}
		runMetaBat.sh {input.merged_assemblies} {input.sorted_bam}
		"""
