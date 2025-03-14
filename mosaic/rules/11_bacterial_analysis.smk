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
		mkdir -p {output.checkMoutdir_temp}
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

rule merge_microbial:
	input:
		assembled_contigs=input_microbial_merge,
	output:
		combined_positive_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial.tot.fasta",
		combined_positive_contigs_2k=dirs_dict["ASSEMBLY_DIR"]+ "/2K_combined_microbial.tot.fasta",
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		derreplicated_tmp=directory(dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot_tmp"),
	params:
		rep_name="combined_microbial_derreplicated_tot_tmp",
		rep_name_full=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot_tmp_rep_seq.fasta",
		rep_temp="combined_microbial_derreplicated_tot_tmp",
		dir_assembly=dirs_dict["ASSEMBLY_DIR"],
		min_len=1000
	message:
		"Derreplicating assembled contigs with mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	threads: 16
	shell:
		"""
		cat {input.assembled_contigs} > {output.combined_positive_contigs}
		seqtk seq -L 2000 {output.combined_positive_contigs} > {output.combined_positive_contigs_2k}
		cd {params.dir_assembly}
		mmseqs easy-cluster --createdb-mode 1 --min-seq-id 1 -c 1 --cov-mode 1 {output.combined_positive_contigs} {params.rep_name} {params.rep_temp} --threads {threads}
		seqtk seq -L 2000 {params.rep_name_full} > {output.derreplicated_microbial_contigs}
		"""
		

rule buildBowtieDB_microbial:
	input:
		combined_positive_contigs_2k=dirs_dict["ASSEMBLY_DIR"]+ "/2K_combined_microbial.tot.fasta",
	output:
		contigs_bt2=dirs_dict["ASSEMBLY_DIR"] + "/2K_combined_microbial.tot.1.bt2",
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/2K_combined_microbial.tot",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/bowtie_microbial.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 64
	shell:
		"""
		bowtie2-build {input.combined_positive_contigs_2k} {params.prefix} --threads {threads}
		"""

rule mapReadsToContigs_microbial:
	input:
		contigs_bt2=dirs_dict["ASSEMBLY_DIR"] + "/2K_combined_microbial.tot.1.bt2",
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
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/2K_combined_microbial.tot",
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

# rule bacterial_binning_metabat_preprocess:
# 	input:
# 		sorted_bam=(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam"),
# 	output:
# 		cov=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted_bam_pileup_coverage.txt"),
# 		abundance=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted_bam_pileup_abundance.txt"),
# 	message:
# 		"Binning microbial contigs with MetaBAT"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/MetaBAT/binning_{sample}.tsv"
# 	threads: 1
# 	shell:
# 		"""
# 		pileup.sh in={input.sorted_bam} out={output.cov} 32bit=t
# 		awk '{{print $1"\t"$5}}' {output.cov} | grep -v '^#' > {output.abundance}
# 		"""

# rule bacterial_binning_metabat:
# 	input:
# 		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
# 		sorted_bam=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam", sample=SAMPLES),
# 	output:
# 		metabat_outdir=directory(dirs_dict["MAPPING_DIR"] + "/MetaBAT_results/"),
# 	message:
# 		"Binning microbial contigs with MetaBAT"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/MetaBAT/binning.tsv"
# 	threads: 32
# 	shell:
# 		"""
# 		mkdir -p {output.metabat_outdir}
# 		cd {output.metabat_outdir}
# 		runMetaBat.sh -t {threads} {input.derreplicated_microbial_contigs} {input.sorted_bam}
# 		"""

# rule bacterial_binning_MaxBin2:
# 	input:
# 		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
# 		abundances=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted_bam_pileup_abundance.txt", sample=SAMPLES),
# 	output:
# 		maxbin_outdir=directory(dirs_dict["MAPPING_DIR"] + "/MaxBin2_results/"),
# 		abund_list=temp(dirs_dict["MAPPING_DIR"] + "/MaxBin2_abundance_list.txt"),
# 	message:
# 		"Binning microbial contigs with MaxBin2"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/MaxBin2/binning.tsv"
# 	threads: 32
# 	shell:
# 		"""
# 		mkdir -p {output.maxbin_outdir}
# 		cd {output.maxbin_outdir}
# 		ls {input.abundances} > {output.abund_list}
# 		run_MaxBin.pl -contig {input.derreplicated_microbial_contigs} -abund_list {output.abund_list} -out {output.maxbin_outdir} -thread {threads}
# 		"""

# rule bacterial_binning_CONCOCT:
# 	input:
# 		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
# 		sorted_bam=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam", sample=SAMPLES),
# 		sorted_bam_index=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam.bai", sample=SAMPLES),
# 	output:
# 		CONCOCT_10k_fasta=temp(dirs_dict["MAPPING_DIR"] + "/CONCOCT_10K_contigs.fasta"),
# 		CONCOCT_10k_bed=temp(dirs_dict["MAPPING_DIR"] + "/CONCOCT_10K_contigs.bed"),
# 		CONCOCT_coverage=temp(dirs_dict["MAPPING_DIR"] + "/CONCOCT_coverage.txt"),
# 		# CONCOCT_outdir=(dirs_dict["MAPPING_DIR"] + "CONCOCT_results/"),
# 		CONCOCT_fasta=directory(dirs_dict["MAPPING_DIR"] + "/CONCOCT_results/CONCOCT_fasta_bins"),
# 		CONCOCT_clustering=(dirs_dict["MAPPING_DIR"] + "/CONCOCT_results/clustering_merged.csv"),
# 	params:
# 		CONCOCT_outdir=(dirs_dict["MAPPING_DIR"] + "/CONCOCT_results/"),
# 	message:
# 		"Binning microbial contigs with CONCOCT"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/CONCOCT_outdir/binning.tsv"
# 	threads: 32
# 	shell:
# 		"""
# 		mkdir -p {params.CONCOCT_outdir}
# 		cd {params.CONCOCT_outdir}
# 		cut_up_fasta.py {input.derreplicated_microbial_contigs} -c 10000 -o 0 --merge_last -b {output.CONCOCT_10k_bed} > {output.CONCOCT_10k_fasta}
# 		concoct_coverage_table.py {output.CONCOCT_10k_bed} {input.sorted_bam} > {output.CONCOCT_coverage}
# 		concoct --composition_file {output.CONCOCT_10k_fasta} --coverage_file {output.CONCOCT_coverage} -b {params.CONCOCT_outdir} -t {threads}
# 		merge_cutup_clustering.py {params.CONCOCT_outdir}/clustering_gt1000.csv > {output.CONCOCT_clustering}
# 		mkdir {output.CONCOCT_fasta}
# 		extract_fasta_bins.py {input.derreplicated_microbial_contigs}  {output.CONCOCT_clustering} --output_path  {output.CONCOCT_fasta}
# 		"""

rule bacterial_binning_VAMB:
	input:
		combined_positive_contigs_2k=dirs_dict["ASSEMBLY_DIR"]+ "/2K_combined_microbial.tot.fasta",
		sorted_bam=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam", sample=SAMPLES_NO_TECHNICAL),
		sorted_bam_index=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam.bai", sample=SAMPLES_NO_TECHNICAL),
	output:
		vamb_bins=directory(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results/all_bins"),
	params:
		vamb_outdir_temp=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results_temp"),
		vamb_outdir=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results"),
		min_votu_len=config['min_votu_length'],
	message:
		"Binning microbial contigs with vamb"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/VAMB_outdir/binning.tsv"
	threads: 64
	shell:
		"""
		rm -rf {params.vamb_outdir_temp}
		vamb -o "_" --outdir {params.vamb_outdir_temp} --fasta {input.combined_positive_contigs_2k}  \
				--bamfiles {input.sorted_bam} --minfasta {params.min_votu_len} -p {threads}
		mkdir {params.vamb_outdir}
		mv {params.vamb_outdir_temp}/* {params.vamb_outdir}
		rm -rf {params.vamb_outdir_temp}
		for d in {params.vamb_outdir}/bins/*/ ; do cp ${{d}}*fna {output.vamb_bins} & done
		"""

# rule polish_bins:
# 	input:
# 		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
# 		metabat_outdir=(dirs_dict["MAPPING_DIR"] + "/MetaBAT_results/"),
# 		CONCOCT_clustering=(dirs_dict["MAPPING_DIR"] + "/CONCOCT_results/clustering_merged.csv"),
# 		maxbin_outdir=(dirs_dict["MAPPING_DIR"] + "/MaxBin2_results/"),
# 	output:
# 		scaffolds2bin_concoct=(dirs_dict["MAPPING_DIR"] + "/concoct_scaffolds2bin.tsv"),
# 		scaffolds2bin_metabat=(dirs_dict["MAPPING_DIR"] + "/metabat_scaffolds2bin.tsv"),
# 		scaffolds2bin_maxbin=(dirs_dict["MAPPING_DIR"] + "/maxbin_scaffolds2bin.tsv"),
# 		DAS_Tool_results=directory(dirs_dict["MAPPING_DIR"] + "/DAS_Tool_results/"),
# 	params:
# 		DAS_Tool_results=(dirs_dict["MAPPING_DIR"] + "/DAS_Tool_results/DAS_Tool_results"),
# 	message:
# 		"Binning microbial contigs with MaxBin2"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/DAS_Tool/binning.tsv"
# 	threads: 64
# 	shell:
# 		"""
# 		perl -pe "s/,/\tconcoct./g;" {input.CONCOCT_clustering} | tail -n +2 > {output.scaffolds2bin_concoct}
# 		Fasta_to_Contig2Bin.sh -i {input.metabat_outdir}/*metabat-bins*/ -e fa > {output.scaffolds2bin_metabat}
# 		Fasta_to_Contig2Bin.sh -i {input.maxbin_outdir} -e fasta > {output.scaffolds2bin_maxbin}
# 		mkdir {output.DAS_Tool_results} 
# 		cd {output.DAS_Tool_results} 
# 		DAS_Tool -i {output.scaffolds2bin_concoct},{output.scaffolds2bin_metabat},{output.scaffolds2bin_maxbin} \
# 			-l concoct,metabat,maxbin -c {input.derreplicated_microbial_contigs} -o {params.DAS_Tool_results} \
# 			--search_engine diamond --threads {threads} --write_bins --write_bin_evals --score_threshold 0
# 		"""

rule predict_spacers:
	input:
		combined_positive_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		minced_dir=(os.path.join(workflow.basedir, config['minced_dir'])),
	output:
		spacers=(dirs_dict["ANNOTATION"] + "/minced_predicted_spacers.tsv"),
	message:
		"Getting CRISPR spacers with MinCED"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/CRISPR/minced.tsv"
	threads: 64
	shell:
		"""
		{input.minced_dir}/minced -spacers {input.combined_positive_contigs} {output.spacers}
		"""

rule estimateBinningQuality:
	input:
		vamb_bins=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results/all_bins/"),
		checkm_db=(config['checkm_db']),
	output:
		checkMoutdir_temp=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM_temp")),
		checkMoutdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM"),
		# checkMoutplots=directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM_plots"),
	params:
		checkm_table=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM/tab_results_checkM.csv"),
		checkm_outfile=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM/output_results_checkM.txt"),
		all_bins=dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM_temp/all_bins",
	log:
		checkMoutdir=(dirs_dict["vOUT_DIR"] + "/microbial_checkM_log"),
	message:
		"Estimating genome completeness with CheckM "
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/estimateGenomeCompletness/microbial_checkm.tsv"
	threads: 32
	shell:
		"""
		mkdir -p {output.checkMoutdir_temp}
		cp -r {input.vamb_bins} {output.checkMoutdir_temp}
		cd {params.all_bins}
		checkm lineage_wf --tab_table -t {threads} -f {params.checkm_outfile} -x fna {params.all_bins} {output.checkMoutdir} 1> {log}
		"""

rule taxonomy_binning:
	input:
		vamb_bins=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results/all_bins/"),
		gtdbtk_db=(config['gtdbtk_db']),
	output:
		GTDB_outdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_vamb_GTDB-Tk"),
		tempdir=temp(directory(dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_dir")),
	params:
		mash_outdir=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_vamb_GTDB-Tk_mash"),
	message:
		"Assigning microbial taxonomy with GTDB-Tk "
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/taxonomy_assignment/microbial_vamb_GTDB-Tk.tsv"
	threads: 64
	shell:
		"""
		mkdir {output.tempdir}
		cd {output.tempdir}
		conda env config vars set GTDBTK_DATA_PATH={input.gtdbtk_db}/release214/
		gtdbtk classify_wf --genome_dir {input.vamb_bins}/ --out_dir {output.GTDB_outdir} --cpus {threads} --mash_db {params.mash_outdir} --extension fna
		"""

rule taxonomy_binning_assembly:
	input:
		racoon_assembly=expand(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample_nanopore}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta", sample_nanopore=NANOPORE_SAMPLES, sampling=SAMPLING_TYPE),
		gtdbtk_db=(config['gtdbtk_db']),
	output:
		GTDB_outdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/assembly_racon_microbial_GTDB-Tk"),
	params:
		mash_outdir=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_racon_GTDB-Tk_mash"),
	message:
		"Assigning microbial taxonomy with GTDB-Tk "
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/taxonomy_assignment/assembly_microbial_GTDB-Tk.tsv"
	threads: 64
	shell:
		"""
		mkdir {output.GTDB_outdir}/input_assemblies
		cp {input.racoon_assembly} {output.GTDB_outdir}/input_assemblies
		conda env config vars set GTDBTK_DATA_PATH={input.gtdbtk_db}/release214/
		gtdbtk classify_wf --genome_dir {output.GTDB_outdir}/input_assemblies/ --out_dir {output.GTDB_outdir} --cpus {threads} --mash_db {params.mash_outdir} --extension fasta
		"""
# rule taxonomy_assembly:
# 	input:
# 		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
# 		gtdbtk_db=(config['gtdbtk_db']),
# 	output:
# 		GTDB_outdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/assembly_microbial_GTDB-Tk"),
# 		GTDB_temp=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/combined_microbial_derreplicated_singlefasta_GTDB-Tk")),
# 	params:
# 		mash_outdir=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_GTDB-Tk_mash"),
# 	message:
# 		"Assigning microbial taxonomy with GTDB-Tk "
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/taxonomy_assignment/assembly_microbial_GTDB-Tk.tsv"
# 	threads: 64
# 	shell:
# 		"""
# 		seqkit split --quiet -i {input.derreplicated_microbial_contigs} --out-dir {output.GTDB_temp}
# 		conda env config vars set GTDBTK_DATA_PATH={input.gtdbtk_db}/release214/
# 		gtdbtk classify_wf --genome_dir {output.GTDB_temp} --out_dir {output.GTDB_outdir} --cpus {threads} --mash_db {params.mash_outdir} --extension fasta
# 		"""

rule DRAM_microbial_annotation:
	input:
		vamb_bins=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results/all_bins/"),
		DRAM_db=config['DRAM_db'],
	output:
		DRAM_output=directory(dirs_dict["ANNOTATION"]+ "/DRAM_annotate_results_{sampling}"),
		DRAM_summary=directory(dirs_dict["ANNOTATION"]+ "/DRAM_distill_results_{sampling}"),
	params:
		DRAM_annotations=dirs_dict["ANNOTATION"]+ "/DRAM_annotate_results_{sampling}/annotations.tsv",
		# trna=directory(dirs_dict["vOUT_DIR"]+ "/DRAM_combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq_{sampling}/trnas.tsv"),
		# rrna=directory(dirs_dict["vOUT_DIR"]+ "/DRAM_combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq_{sampling}/rrnas.tsv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/DRAM/{sampling}.tsv"
	message:
		"Annotate contigs with DRAM"
	threads: 32
	shell:
		"""
		DRAM.py annotate -i '{input.DAS_Tool_bins}/*fna' -o {output.DRAM_output} --threads {threads}
		DRAM.py distill -i {params.DRAM_annotations} -o {output.DRAM_summary} 
		"""


# rule taxonomy_assembly:
# 	input:
# 		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
# 		gtdbtk_db=(config['gtdbtk_db']),
# 	output:
# 		GTDB_outdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/assembly_microbial_GTDB-Tk"),
# 		GTDB_temp=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/combined_microbial_derreplicated_singlefasta_GTDB-Tk")),
# 	params:
# 		mash_outdir=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_GTDB-Tk_mash"),
# 	message:
# 		"Assigning microbial taxonomy with GTDB-Tk "
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/taxonomy_assignment/assembly_microbial_GTDB-Tk.tsv"
# 	threads: 64
# 	shell:
# 		"""
# 		seqkit split --quiet -i {input.derreplicated_microbial_contigs} --out-dir {output.GTDB_temp}
# 		conda env config vars set GTDBTK_DATA_PATH={input.gtdbtk_db}/release214/
# 		gtdbtk classify_wf --genome_dir {output.GTDB_temp} --out_dir {output.GTDB_outdir} --cpus {threads} --mash_db {params.mash_outdir} --extension fasta
# 		"""

rule single_fasta_microbial:
	input:
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
	output:
		derreplicated_microbial_contigs_dir=temp(directory(dirs_dict["ASSEMBLY_DIR"]+ "/single_combined_microbial_derreplicated_tot")),
	message:
		"formating microbial contigs into single fasta"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	threads: 1
	shell:
		"""
		seqkit split --quiet -i {input.derreplicated_microbial_contigs} --out-dir {output.derreplicated_microbial_contigs_dir}
	 	"""

rule sourmash_sketch_microbial:
	input:
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		derreplicated_microbial_contigs_dir=((dirs_dict["ASSEMBLY_DIR"]+ "/single_combined_microbial_derreplicated_tot")),
	output:
		manysketch_csv=temp(dirs_dict["ANNOTATION"] + "/combined_microbial_derreplicated_tot_manysketch.csv"),
		sketch=temp(dirs_dict["ANNOTATION"] + "/combined_microbial_derreplicated_tot_sourmash.sig.zip"),
	params: 
		name="combined_microbial_derreplicated_tot"
	message:
		"Building sketches with sourmash"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/sourmash/combined_microbial_derreplicated_tot_sketch.tsv"
	threads: 64
	shell:
		"""
		echo name,genome_filename,protein_filename > {output.manysketch_csv}
		grep "^>" {input.derreplicated_microbial_contigs} | sed 's/^>//' | awk -v dir="{input.derreplicated_microbial_contigs_dir}/{params.name}.part_" '{{print $1 "," dir $1 ".fasta,"}}' >> {output.manysketch_csv}
		sourmash scripts manysketch {output.manysketch_csv} -p k=31,abund -o {output.sketch} -c {threads}
		"""

rule sourmash_gather_microbial:
	input:
		sketch=(dirs_dict["ANNOTATION"] + "/combined_microbial_derreplicated_tot_sourmash.sig.zip"),
		sourmash_rocksdb=config['sourmash_rocksdb'],
	output:
		gather=temp(dirs_dict["ANNOTATION"] + "/combined_microbial_derreplicated_tot_gather_sourmash.csv"),
	message:
		"Genome containtment with sourmash gather"
	params:
		threshold_bp="0"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/sourmash/combined_microbial_derreplicated_tot_gather.tsv"
	threads: 64
	shell:
		"""
		sourmash scripts fastmultigather {input.sketch} {input.sourmash_rocksdb} -c {threads} -o {output.gather} -t {params.threshold_bp} -s 1000
		"""

rule sourmash_tax_microbial:
	input:
		gather=(dirs_dict["ANNOTATION"] + "/combined_microbial_derreplicated_tot_gather_sourmash.csv"),
		sourmash_tax=config['sourmash_tax'],
	output:
		csv_report=(dirs_dict["ANNOTATION"] + "/sourmash_combined_microbial_derreplicated_tot.classifications.csv"),
	params:
		outdir=(dirs_dict["ANNOTATION"]),
		name="sourmash_combined_microbial_derreplicated_tot",
	message:
		"Assigning taxonomy with sourmash tax"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/sourmash/combined_microbial_derreplicated_tot_tax.tsv"
	threads: 1
	shell:
		"""
		sourmash tax genome --gather-csv {input.gather} -t {input.sourmash_tax}  -o {params.name}\
			--output-dir {params.outdir} -F csv_summary
		"""
		
rule defense_finder:
	input:
		aa=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_ORFs_tot.faa",
	output:
		defenseFinder_dir=directory(dirs_dict["ANNOTATION"]+ "/DefenseFinder_results_{sampling}/"),
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/DefenseFinder/{sampling}.tsv"
	message:
		"Detecting defense systems with DefenseFinder"
	threads: 32
	shell:
		"""
		python scripts/process_fasta_satellite_finder.py {params.faa} {output.faa_temp}
		defense-finder run -–db-type gembase -w {threads} --out-dir {output.defenseFinder_dir} {input.aa}
		"""


# rule getORFs_microbial:
# 	input:
# 		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
# 	output:
# 		coords=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_ORFs_tot.coords",
# 		aa=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_ORFs_tot.faa",
# 	message:
# 		"Calling ORFs with prodigal"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		prodigal -i {input.derreplicated_microbial_contigs} -o {output.coords} -a {output.aa} -p meta
# 		"""