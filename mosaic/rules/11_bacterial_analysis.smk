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

def input_microbial_merge(wildcards):
	input_list=[]
	input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.tot.fasta", sample=SAMPLES)),
	if len(config['additional_reference_contigs'])>0:
		input_list.append(config['additional_reference_contigs'])
	return input_list

# rule merge_assembly:
# 	input:
# 		assembled_contigs=input_microbial_merge,
# 	output:
#		 combined_microbial=
# 		merged_assemblies=dirs_dict["ASSEMBLY_DIR"] + "/merged_microbial_assembly.tot.fasta",
# 	message:
# 		"Merging microbial assemblies"
# 	shell:
# 		"""
# 		cat {input.assembled_contigs} > {output.merged_assemblies}
# 		"""

rule derreplicate_microbial:
	input:
		assembled_contigs=input_microbial_merge,
	output:
		combined_positive_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial.tot.fasta",
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		derreplicated_tmp=directory(dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot_tmp"),
	params:
		rep_name="combined_microbial_derreplicated_tot_tmp",
		rep_name_full=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot_tmp_rep_seq.fasta",
		rep_temp="combined_microbial_derreplicated_tot_tmp",
		dir_assembly=dirs_dict["ASSEMBLY_DIR"],
	message:
		"Derreplicating assembled contigs with mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	threads: 16
	shell:
		"""
		cat {input.assembled_contigs} > {output.combined_positive_contigs}
		cd {params.dir_assembly}
		mmseqs easy-cluster --createdb-mode 1 --min-seq-id 1 -c 1 --cov-mode 1 {output.combined_positive_contigs} {params.rep_name} {params.rep_temp} --threads {threads}
		mv {params.rep_name_full} {output.derreplicated_microbial_contigs}
		"""


rule buildBowtieDB_microbial:
	input:
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
	output:
		contigs_bt2=dirs_dict["ASSEMBLY_DIR"] + "/combined_microbial_derreplicated_tot.1.bt2",
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/combined_microbial_derreplicated_tot",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/bowtie_microbial.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		bowtie2-build {input.derreplicated_microbial_contigs} {params.prefix} --threads {threads}
		"""

rule mapReadsToContigs_microbial:
	input:
		contigs_bt2=dirs_dict["ASSEMBLY_DIR"] + "/combined_microbial_derreplicated_tot.1.bt2",
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
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
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/combined_microbial_derreplicated_tot",
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

rule bacterial_binning_metabat:
	input:
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		sorted_bam=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam", sample=SAMPLES),
	output:
		metabat_outdir=directory(dirs_dict["MAPPING_DIR"] + "/MetaBAT_results/"),
	message:
		"Binning microbial contigs with MetaBAT"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/MetaBAT/binning.tsv"
	threads: 64
	shell:
		"""
		mkdir {output.metabat_outdir}
		cd {output.metabat_outdir}
		runMetaBat.sh -t {threads} {input.derreplicated_microbial_contigs} {input.sorted_bam}
		"""
rule bacterial_binning_MaxBin2:
	input:
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		sorted_bam=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam", sample=SAMPLES),
	output:
		maxbin_outdir=directory(dirs_dict["MAPPING_DIR"] + "/MaxBin2_results/"),
	message:
		"Binning microbial contigs with MaxBin2"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/MaxBin2/binning.tsv"
	threads: 32
	shell:
		"""
		mkdir {output.maxbin_outdir}
		cd {output.maxbin_outdir}
		run_MaxBin.pl -contig {input.derreplicated_microbial_contigs} -reads {input.sorted_bam} -out {output.maxbin_outdir} -thread {threads}
		"""

rule bacterial_binning_CONCOCT:
	input:
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		sorted_bam=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam", sample=SAMPLES),
		sorted_bam_index=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam.bai", sample=SAMPLES),
	output:
		CONCOCT_10k_fasta=temp(dirs_dict["MAPPING_DIR"] + "/CONCOCT_10K_contigs.fasta"),
		CONCOCT_10k_bed=temp(dirs_dict["MAPPING_DIR"] + "/CONCOCT_10K_contigs.bed"),
		CONCOCT_coverage=temp(dirs_dict["MAPPING_DIR"] + "/CONCOCT_coverage.txt"),
		# CONCOCT_outdir=(dirs_dict["MAPPING_DIR"] + "CONCOCT_results/"),
		CONCOCT_fasta=directory(dirs_dict["MAPPING_DIR"] + "/CONCOCT_results/CONCOCT_fasta_bins"),
		CONCOCT_clustering=(dirs_dict["MAPPING_DIR"] + "/CONCOCT_results/clustering_merged.csv"),
	params:
		CONCOCT_outdir=(dirs_dict["MAPPING_DIR"] + "/CONCOCT_results/"),
	message:
		"Binning microbial contigs with CONCOCT"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/CONCOCT_outdir/binning.tsv"
	threads: 32
	shell:
		"""
		mkdir {params.CONCOCT_outdir}
		cd {params.CONCOCT_outdir}
		cut_up_fasta.py {input.derreplicated_microbial_contigs} -c 10000 -o 0 --merge_last -b {output.CONCOCT_10k_bed} > {output.CONCOCT_10k_fasta}
		concoct_coverage_table.py {output.CONCOCT_10k_bed} {input.sorted_bam} > {output.CONCOCT_coverage}
		concoct --composition_file {output.CONCOCT_10k_fasta} --coverage_file {output.CONCOCT_coverage} -b {params.CONCOCT_outdir} -threads {threads}
		merge_cutup_clustering.py {params.CONCOCT_outdir}/clustering_gt1000.csv > {output.CONCOCT_clustering}
		extract_fasta_bins.py {input.derreplicated_microbial_contigs}  {output.CONCOCT_clustering} --output_path  {output.CONCOCT_fasta}
		"""

rule polish_bins:
	input:
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		metabat_outdir=(dirs_dict["MAPPING_DIR"] + "/MetaBAT_results/"),
		CONCOCT_clustering=(dirs_dict["MAPPING_DIR"] + "/CONCOCT_results/clustering_merged.csv"),
		maxbin_outdir=(dirs_dict["MAPPING_DIR"] + "/MaxBin2_results/"),
	output:
		scaffolds2bin_concoct=(dirs_dict["MAPPING_DIR"] + "/concoct_scaffolds2bin.tsv"),
		scaffolds2bin_metabat=(dirs_dict["MAPPING_DIR"] + "/metabat_scaffolds2bin.tsv"),
		scaffolds2bin_maxbin=(dirs_dict["MAPPING_DIR"] + "/maxbin_scaffolds2bin.tsv"),
		DAS_Tool_results=directory(dirs_dict["MAPPING_DIR"] + "/DAS_Tool_results/"),
	message:
		"Binning microbial contigs with MaxBin2"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/MaxBin2/binning.tsv"
	threads: 64
	shell:
		"""
		perl -pe "s/,/\tconcoct./g;" {input.CONCOCT_clustering} > {output.scaffolds2bin_concoct}
		Fasta_to_Scaffolds2Bin.sh -i {input.metabat_outdir}/*metabat-bins*/ -e fa > {output.scaffolds2bin_metabat}
		Fasta_to_Scaffolds2Bin.sh -i {input.metabat_outdir} -e fasta > {output.scaffolds2bin_maxbin}
		DAS_Tool-.sif -i {output.scaffolds2bin_concoct},{output.scaffolds2bin_metabat},{output.scaffolds2bin_maxbin} \
		 -l concoct,metabat,maxbin -c {input.derreplicated_microbial_contigs} -o {output.DAS_Tool_results} --search_engine diamond --threads {threads}
		"""


rule estimateBinningQuality:
	input:
		DAS_Tool_results=(dirs_dict["MAPPING_DIR"] + "/DAS_Tool_results/"),
		checkm_db=(config['checkm_db']),
	output:
		checkMoutdir_temp=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM_temp")),
		checkMoutdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM"),
		checkMoutplots=directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM_plots"),
	params:
		checkm_table=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM/tab_results_checkM.csv"),
		checkm_outfile=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM/output_results_checkM.txt"),
	log:
		checkMoutdir=(dirs_dict["vOUT_DIR"] + "/microbial_checkM_log"),
	message:
		"Estimating genome completeness with CheckM "
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/estimateGenomeCompletness/microbial_checkm.tsv"
	threads: 16
	shell:
		"""
		mkdir {output.checkMoutdir_temp}
		cp -r {input.DAS_Tool_results}/* {output.checkMoutdir_temp}
		cd {output.checkMoutdir_temp}
		checkm lineage_wf -t {threads} --tab_table {params.checkm_table} -f {params.checkm_outfile} -x fa {output.checkMoutdir_temp} {output.checkMoutdir}  1> {log}
		"""

rule taxonomy_binning:
	input:
		metabat_outdir=(dirs_dict["MAPPING_DIR"] + "/MetaBAT_results/"),
		gtdbtk_db=(config['gtdbtk_db']),
	output:
		GTDB_outdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_GTDB-Tk"),
	params:
		mash_outdir=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_GTDB-Tk_mash"),
	message:
		"Assigning microbial taxonomy with GTDB-Tk "
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/taxonomy_assignment/microbial_GTDB-Tk.tsv"
	threads: 64
	shell:
		"""
		gtdbtk classify_wf --genome_dir {input.metabat_outdir}/*metabat-bins* --out_dir {output.GTDB_outdir} --cpus {threads} --mash_db {params.mash_outdir} --extension fa
		"""


