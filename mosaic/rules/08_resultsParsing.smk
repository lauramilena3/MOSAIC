def benchmark_snapshot(wildcards):
	from pathlib import Path
	root=Path(dirs_dict["BENCHMARKS"])
	return [(str(path.relative_to(root)), path.stat().st_size, path.stat().st_mtime_ns)
		for path in sorted(root.rglob("*.tsv"))
		if path.relative_to(root).parts[0] != "benchmark_summary"]

def benchmark_registry(wildcards):
	registry=[]
	for rule in workflow.rules:
		if rule.benchmark is None or rule.name == "benchmark_summary":
			continue
		threads=rule.resources.get("_cores", 1)
		registry.append(dict(name=rule.name, pattern=str(rule.benchmark),
			threads_declared=threads if isinstance(threads, (int, float)) else None,
			inputs=[str(value) for value in rule.input if not callable(value)]))
	return registry

rule abundance_lifestyle_summary:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		summary=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_summary.tot.csv",
		rpkm=dirs_dict["MAPPING_DIR"] + "/filtered_RPKM_normalised_tot.txt",
		counts=dirs_dict["MAPPING_DIR"] + "/filtered_counts_normalised_tot.txt",
		coverage=dirs_dict["MAPPING_DIR"] + "/breadth_coverage_percent_tot.txt",
		composition=dirs_dict["ANNOTATION"] + "/nucleotide_content_viral_contigs_clustered_with_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.tsv",
		bacphlip=dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot_bacphlip.csv",
		phagcn=dirs_dict["ANNOTATION"] + "/PhaGCN_taxonomy_report_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.csv",
		taxmyphage=dirs_dict["ANNOTATION"] + "/taxmyphage_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot",
		iphop=dirs_dict["ANNOTATION"] + "/iphop_hostID_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot_resultsDir",
		refseq=dirs_dict["ANNOTATION"] + "/blast_output_ViralRefSeq_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.csv",
		metavr=[dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.tot/vOTU_summary.tsv"] if METAVR_blast else [],
		crispr=[dirs_dict["ANNOTATION"] + "/spacepharer_minced_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.tsv"] if config["microbial_spacers"] else [],
		crispr_hosts=[config["abundance_crispr_host_metadata"]] if config.get("abundance_crispr_host_metadata", "") else [],
		sample_metadata=[config["abundance_sample_metadata"]] if config.get("abundance_sample_metadata", "") else [],
	output:
		html=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot.html",
		metadata=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/vOTU_metadata_polished.csv",
		samples=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/sample_summary.csv",
		relative=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/relative_abundance.csv",
		presence=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/presence_absence.csv",
		prevalence=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/vOTU_prevalence.csv",
		categories=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/category_summary.csv",
		genomes=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/genome_summary.csv",
		accumulation=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/sample_accumulation.csv",
		groups=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/group_summary.csv",
		group_presence=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/group_ubiquitous_vOTUs.csv",
		figures=directory(dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/figures"),
	params:
		samples=SAMPLES,
		group_column=config.get("abundance_group_column", ""),
		common_prevalence=float(config.get("abundance_common_prevalence", 0.2)),
		core_prevalence=float(config.get("abundance_core_prevalence", 0.5)),
	message:
		"Summarizing final vOTU abundance, quality, lifestyle, hosts and prevalence"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/abundance_lifestyle_summary/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/10_abundance_lifestyle_summary.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/10_abundance_lifestyle_summary.py.ipynb"

rule microdiversity_summary:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		summaries=expand(MICRO_DIR + "/{sample}/summary.tsv", sample=SAMPLES),
	output:
		html=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot.html",
		combined=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/vOTU_sample_microdiversity.csv",
		samples=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/sample_summary.csv",
		pi=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/pi_callable.csv",
		callable_fraction=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/callable_fraction.csv",
		figures=directory(dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/figures"),
	params:
		samples=SAMPLES,
	message:
		"Summarizing callable-site nucleotide diversity across samples"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/microdiversity_summary/tot.tsv"
	threads: 1
	resources:
		mem_mb=4000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/11_microdiversity_summary.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/11_microdiversity_summary.py.ipynb"

rule benchmark_summary:
	output:
		jobs=dirs_dict["PLOTS_DIR"] + "/09_benchmark_jobs.csv",
		rules=dirs_dict["PLOTS_DIR"] + "/09_benchmark_rules.csv",
		figure=dirs_dict["PLOTS_DIR"] + "/09_benchmark_summary.png",
		html=dirs_dict["PLOTS_DIR"] + "/09_benchmark_summary.html",
	params:
		benchmark_dir=dirs_dict["BENCHMARKS"],
		registry=benchmark_registry,
		snapshot=benchmark_snapshot,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/benchmark_summary/tot.tsv"
	threads: 1
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/09_benchmark_summary.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/09_benchmark_summary.py.ipynb"

rule plot_assemblies:
	input:
		aa="{fasta}_ORFs.{sampling}.fasta",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
		corrected1=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		corrected2=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		corrected3=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_3_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		corrected4=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_4_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		corrected_medaka=dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		scaffolds_pilon1=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_1_{sampling}/pilon.faa"),
		scaffolds_pilon2=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_2_{sampling}/pilon.faa"),
		scaffolds_pilon3=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_3_{sampling}/pilon.faa"),
		scaffolds_pilon4=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_4_{sampling}/pilon.faa"),
	output:
		plot=(dirs_dict["CLEAN_DATA_DIR"] + "/protein_lengths_plot.{sampling}.png"),
		svg=(dirs_dict["CLEAN_DATA_DIR"] + "/protein_lengths_plot.{sampling}.svg"),
	message:
		"Plot unique reads with BBtools"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/plot_assemblies/sampling={sampling}.tsv"
	threads: 1
	run:
		import pandas as pd
		import seaborn as sns; sns.set()
		import matplotlib.pyplot as plt

		plt.figure(figsize=(12,12))
		sns.set(font_scale=2)
		sns.set_style("whitegrid")

		read_max=0

		for h in input.histograms:
			df=pd.read_csv(h, sep="\t")
			df.columns=["count", "percent", "c", "d", "e", "f", "g", "h", "i", "j"]
			df=df[["count", "percent"]]
			ax = sns.lineplot(x="count", y="percent", data=df,err_style='band', label=h.split("/")[-1].split("_kmer")[0])
			read_max=max(read_max,df["count"].max())

		ax.set(ylim=(0, 100))
		ax.set(xlim=(0, read_max*1.2))

		ax.set_xlabel("Read count",fontsize=20)
		ax.set_ylabel("New k-mers (%)",fontsize=20)
		ax.figure.savefig(output.plot)
		ax.figure.savefig(output.svg, format="svg")

rule QC_parsing:
	input:
		inputReadsCount,
		histograms=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.{{sampling}}.csv", sample=SAMPLES),
		preqc_txt=dirs_dict["QC_DIR"]+ "/preQC_illumina_report_data/multiqc_fastqc.txt",
		postqc_txt=dirs_dict["QC_DIR"]+ "/postQC_illumina_report_data/multiqc_fastqc.txt",
		read_count_raw_forward=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + "_read_count.txt", sample=SAMPLES),
		read_count_raw_reverse=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + "_read_count.txt", sample=SAMPLES),
		read_count_trimmed_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_read_count.txt", sample=SAMPLES),
		read_count_trimmed_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_read_count.txt", sample=SAMPLES),
		read_count_duk_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot_read_count.txt", sample=SAMPLES),
		read_count_duk_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot_read_count.txt", sample=SAMPLES),
    	read_count_norm_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.tot_read_count.txt", sample=SAMPLES),
		read_count_norm_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.tot_read_count.txt", sample=SAMPLES),
		supper_dedup=expand(dirs_dict["QC_DIR"] + "/{sample}_stats_pcr_duplicates.log", sample=SAMPLES),
		histogram_kmer_pre=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_count_histogram_pre.tot.txt", sample=SAMPLES),
		histogram_kmer_post=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_count_histogram_post.tot.txt", sample=SAMPLES),
		peak_kmer=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_count_peaks.tot.txt", sample=SAMPLES),
	output:
		kmer_png=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_plot.{sampling}.png"),
		kmer_svg=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_plot.{sampling}.svg"),
		kmer_fit_png=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_plot_fitted.{sampling}.png"),
		kmer_fit_svg=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_plot_fitted.{sampling}.svg"),
		kmer_fit_html=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_fitted.{sampling}.html"),
		kmer_dist_pre_png=(dirs_dict["PLOTS_DIR"] + "/01_kmer_distribution_plot_pre.{sampling}.png"),
		kmer_dist_pre_svg=(dirs_dict["PLOTS_DIR"] + "/01_kmer_distribution_plot_pre.{sampling}.svg"),
		kmer_dist_post_png=(dirs_dict["PLOTS_DIR"] + "/01_kmer_distribution_plot_post.{sampling}.png"),
		kmer_dist_post_svg=(dirs_dict["PLOTS_DIR"] + "/01_kmer_distribution_plot_post.{sampling}.svg"),
		qc_summary_html=(dirs_dict["PLOTS_DIR"] + "/01_post_qc_read_summary.{sampling}.html"),
		percentage_kept_reads_png=(dirs_dict["PLOTS_DIR"] + "/01_percentage_kept_reads.{sampling}.png"),
		percentage_kept_reads_svg=(dirs_dict["PLOTS_DIR"] + "/01_percentage_kept_reads.{sampling}.svg"),
		percentage_kept_Mbp_png=(dirs_dict["PLOTS_DIR"] + "/01_percentage_kept_Mbp.{sampling}.png"),
		percentage_kept_Mbp_svg=(dirs_dict["PLOTS_DIR"] + "/01_percentage_kept_Mbp.{sampling}.svg"),
		step_qc_reads_html=(dirs_dict["PLOTS_DIR"] + "/01_multistep_qc_report.{sampling}.html"),
		steps_qc_reads_png=(dirs_dict["PLOTS_DIR"] + "/01_qc_bystep_counts.{sampling}.png"),
		steps_qc_reads_svg=(dirs_dict["PLOTS_DIR"] + "/01_qc_bystep_counts.{sampling}.svg"),
		steps_qc_percentage_png=(dirs_dict["PLOTS_DIR"] + "/01_qc_bystep_percentage.{sampling}.png"),
		steps_qc_percentage_svg=(dirs_dict["PLOTS_DIR"] + "/01_qc_bystep_percentage.{sampling}.svg"),
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.{sampling}.csv",
		supperdedupper_html=(dirs_dict["PLOTS_DIR"] + "/01_superdedupper_PCR.{sampling}.html"),
		supperdedupper_png=(dirs_dict["PLOTS_DIR"] + "/01_superdedupper_PCR.{sampling}.png"),
		supperdedupper_svg=(dirs_dict["PLOTS_DIR"] + "/01_superdedupper_PCR.{sampling}.svg"),
	params:
		results_dir=RESULTS_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		samples=SAMPLES,
		forward_tag=config['forward_tag'],
		reverse_tag=config['reverse_tag'],
		raw_dir=dirs_dict["RAW_DATA_DIR"],
		qc_dir=dirs_dict["QC_DIR"],
		remove_euk=REMOVE_EUK,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/QC_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/01_QC.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/01_QC.py.ipynb"


rule assembly_parsing_short:
	input:
		quast_report_dir=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}",
	output:
		log_number_contigs_png=(dirs_dict["PLOTS_DIR"] + "/03_log_number_contigs_plot.{sampling}.png"),
		log_number_contigs_svg=(dirs_dict["PLOTS_DIR"] + "/03_log_number_contigs_plot.{sampling}.svg"),
		contig_length_bp_png=(dirs_dict["PLOTS_DIR"] + "/03_contig_length_bp_plot.{sampling}.png"),
		contig_length_bp_svg=(dirs_dict["PLOTS_DIR"] + "/03_contig_length_bp_plot.{sampling}.svg"),
		contig_number_total_png=(dirs_dict["PLOTS_DIR"] + "/03_contig_number_total_plot.{sampling}.png"),
		contig_number_total_svg=(dirs_dict["PLOTS_DIR"] + "/03_contig_number_total_plot.{sampling}.svg"),
		contig_length_total_png=(dirs_dict["PLOTS_DIR"] + "/03_contig_length_total_plot.{sampling}.png"),
		contig_length_total_svg=(dirs_dict["PLOTS_DIR"] + "/03_contig_length_total_plot.{sampling}.svg"),
	params:
		input_quast_report=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}/transposed_report.tsv",
		results_dir=RESULTS_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		samples=SAMPLES,
		forward_tag=config['forward_tag'],
		reverse_tag=config['reverse_tag'],
		raw_dir=dirs_dict["RAW_DATA_DIR"],
		qc_dir=dirs_dict["QC_DIR"],
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assembly_parsing_short/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/03_assembly_short.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/03_assembly_short.py.ipynb"

rule assembly_parsing_short_RNA:
	input:
		assemblies=expand(RNA_DIR + "/{sample}/{assembler}.fasta", sample=SAMPLES, assembler=RNA_ASSEMBLERS),
		combined=expand(RNA_DIR + "/{sample}/combined.fasta", sample=SAMPLES),
		provenance=expand(RNA_DIR + "/{sample}/assembly_provenance.tsv", sample=SAMPLES),
	output:
		summary=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_summary.tot.csv",
		fate=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_provenance.tot.csv",
		contig_counts_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_contig_counts.tot.png",
		contig_counts_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_contig_counts.tot.svg",
		total_length_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_total_length.tot.png",
		total_length_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_total_length.tot.svg",
		length_distribution_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_length_distribution.tot.png",
		length_distribution_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_length_distribution.tot.svg",
		provenance_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_provenance.tot.png",
		provenance_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_provenance.tot.svg",
	params:
		samples=SAMPLES,
		assemblers=RNA_ASSEMBLERS,
		min_length=int(config.get("rna_min_contig_length", 500)),
	message:
		"Summarizing RNA assemblies across samples and assemblers"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assembly_parsing_short_RNA/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/03_assembly_short_RNA.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/03_assembly_short_RNA.py.ipynb"

rule assembly_parsing_long:
	input:
		caudovirales=("db/caudovirales_orf_lengths_09_05_2023.txt"),
		hybrid=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds_ORFs_length.tot.txt"),
		canu=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER +"_ORFs_length.tot.txt"),
		medaka=(dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		racon1=(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_1_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		racon2=(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		scaffolds_pilon1_final=(dirs_dict["ASSEMBLY_DIR"] + "/pilon_1_polished_{sample}_contigs_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		scaffolds_pilon2_final=(dirs_dict["ASSEMBLY_DIR"] + "/pilon_2_polished_{sample}_contigs_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		scaffolds_pilon3_final=(dirs_dict["ASSEMBLY_DIR"] + "/pilon_3_polished_{sample}_contigs_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		scaffolds_pilon4_final=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon_ORFs_length.tot.txt"),
	output:
		orf_length_png=(dirs_dict["PLOTS_DIR"] + "/03_ORF_length_{sample}.png"),
		orf_length_svg=(dirs_dict["PLOTS_DIR"] + "/03_ORF_length_{sample}.svg"),
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assembly_parsing_long/sample={sample}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/03_assembly_long_{sample}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/03_assembly_long.py.ipynb"

def inputAssemblyContigs(wildcards):
	inputs=[]
	inputs.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES))
	if NANOPORE:
		inputs.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample_nanopore}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{{sampling}}.fasta", sample_nanopore=NANOPORE_SAMPLES))
	if CROSS_ASSEMBLY:
		inputs.append(dirs_dict["ASSEMBLY_DIR"] + "/ALL_spades_filtered_scaffolds.{sampling}.fasta")
	if SUBASSEMBLY:
		inputs.extend(expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES, subsample=subsample_test)),
	return inputs

rule viralID_parsing:
	input:
		viral_sequences=input_vOTU_clustering,
		assembled_sequences=inputAssemblyContigs,
	output:
		viral_sequences_count_plot_png=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_count_{sampling}.png"),
		viral_sequences_count_plot_svg=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_count_{sampling}.svg"),
		viral_sequences_count_table_html=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_count_{sampling}.html"),
		viral_sequences_length_plot_png=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_length_{sampling}.png"),
		viral_sequences_length_plot_svg=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_length_{sampling}.svg"),
		viral_sequences_length_table_html=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_length_{sampling}.html"),
	params:
		samples=SAMPLES,
		contig_dir=dirs_dict["ASSEMBLY_DIR"],
		viral_dir=dirs_dict['VIRAL_DIR'],
		subassembly=SUBASSEMBLY,
		cross_assembly=CROSS_ASSEMBLY,
		long_assembler=LONG_ASSEMBLER
	benchmark:
		dirs_dict["BENCHMARKS"] + "/viralID_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/04_viral_ID_{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/04_viral_ID.py.ipynb"

def input_assembly_flagstats(wildcards):
	inputs=[]
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_filtered_{sample}.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT)),
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_assembled_contigs.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT)),
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_viral_contigs.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT)),
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_unfiltered_contigs.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT)),
	return inputs

rule mapping_statistics_parsing:
	input:
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.{sampling}.csv",
		assembled_sequences=inputAssemblyContigs,
		assembly_flagstats=input_assembly_flagstats,
		all_assembled_mapped_pairs=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_mapped_pairs_filtered_AllAssembled_{sample}.tot.txt", sample=SAMPLES) if MAP_TO_ALL_ASSEMBLED else [],
	output:
		mapping_stats_html=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_{sampling}.html"),
		filtered_viral_png=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_filtered_viral_{sampling}.png"),
		filtered_viral_svg=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_filtered_viral_{sampling}.svg"),
		filtered_unfiltered_png=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_filtered_unfiltered_{sampling}.png"),
		filtered_unfiltered_svg=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_filtered_unfiltered_{sampling}.svg"),
	params:
		samples=SAMPLES,
		mapping_dir=dirs_dict["MAPPING_DIR"],
		sampling="{sampling}",
		map_to_all_assembled=MAP_TO_ALL_ASSEMBLED,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapping_statistics_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_mapping_statistics_{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_mapping_statistics.py.ipynb"


def input_phage_isolates_assembled_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs.{sampling}_covstats.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_viral_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs.{sampling}_covstats.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_unfiltered_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs.{sampling}_covstats.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_filtered_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_covstats.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_filtered_flagstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_filtered_{sample}.{sampling}.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_assembled_flagstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_assembled_contigs.{sampling}.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_viral_flagstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_viral_contigs.{sampling}.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_unfiltered_flagstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_unfiltered_contigs.{sampling}.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_host_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_filtered_{sample}_vs_{host}_covstats.txt", host=HOSTS, sample=SAMPLES)


def input_phage_isolates_host_masked_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_filtered_{sample}_vs_{host}_masked_prophages_covstats.txt", host=HOSTS, sample=SAMPLES)


def input_phage_isolates_host_blast(wildcards):
	return expand(dirs_dict["vOUT_DIR"] + "/blastn_out_assembly_{host}.{sampling}.csv", host=HOSTS, sampling=wildcards.sampling)


def input_phage_isolates_host_fastas(wildcards):
	return expand(dirs_dict["HOST_DIR"] + "/{host}.fasta", host=HOSTS)


rule phage_isolates_summary:
	input:
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.{sampling}.csv",
		quast=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}/transposed_report.tsv",
		checkv=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		vibrant_positive=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_list.{sampling}.csv",
		virsorter_positive=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/positive_VS_list_{sampling}.txt",
		viral_refseq_blast=dirs_dict["ANNOTATION"] + "/blast_output_ViralRefSeq_combined_positive_viral_contigs.{sampling}.csv",
		nucleotide_content=dirs_dict["ANNOTATION"]+ "/nucleotide_content_combined_positive_viral_contigs.{sampling}.tsv",
		clusters=dirs_dict["vOUT_DIR"]+ "/new_references_clusters.{sampling}.csv",
		combined_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		aai_distance_matrix=dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_distance_matrix_AAI.txt",
		rpkm=dirs_dict["MAPPING_DIR"] + "/RPKM_normalised_{sampling}.txt",
		filtered_covstats=input_phage_isolates_filtered_covstats,
		assembled_covstats=input_phage_isolates_assembled_covstats,
		viral_covstats=input_phage_isolates_viral_covstats,
		unfiltered_covstats=input_phage_isolates_unfiltered_covstats,
		filtered_flagstats=input_phage_isolates_filtered_flagstats,
		assembled_flagstats=input_phage_isolates_assembled_flagstats,
		viral_flagstats=input_phage_isolates_viral_flagstats,
		unfiltered_flagstats=input_phage_isolates_unfiltered_flagstats,
		host_covstats=input_phage_isolates_host_covstats,
		host_masked_covstats=input_phage_isolates_host_masked_covstats,
		host_blast=input_phage_isolates_host_blast,
		host_fastas=input_phage_isolates_host_fastas
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_summary.{sampling}.html",
		summary_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_summary.{sampling}.csv",
		contig_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_contigs.{sampling}.csv",
		closest_relatives_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_closest_relatives.{sampling}.csv",
		host_mapping_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_mapping.{sampling}.csv",
		cluster_members_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_members.{sampling}.csv",
		cluster_summary_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_summary.{sampling}.csv",
		single_contig_samples_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_single_contig_samples.{sampling}.csv",
		host_blast_summary_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_blast_summary.{sampling}.csv",
		cluster_coverage_matrix_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_coverage_matrix.{sampling}.csv",
		host_blast_coverage_matrix_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_blast_coverage_matrix.{sampling}.csv",
		percent_covered_matrix_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_percent_covered_matrix.{sampling}.csv",
		covstats_rpkm_matrix_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_covstats_rpkm_matrix.{sampling}.csv",
		aai_cluster_index_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_aai_cluster_index.{sampling}.csv",
		aai_cluster_dir=directory(dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_aai_clusters.{sampling}"),
		decisions_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_decisions.{sampling}.png",
		decisions_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_decisions.{sampling}.svg",
		remaining_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_remaining_contigs.{sampling}.png",
		remaining_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_remaining_contigs.{sampling}.svg",
		host_viral_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_viral_contigs.{sampling}.png",
		host_viral_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_viral_contigs.{sampling}.svg",
		completeness_coverage_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_completeness_vs_coverage.{sampling}.png",
		completeness_coverage_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_completeness_vs_coverage.{sampling}.svg",
		dominant_signal_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_dominant_signal.{sampling}.png",
		dominant_signal_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_dominant_signal.{sampling}.svg",
		rpkm_heatmap_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_votu_rpkm_heatmap.{sampling}.png",
		rpkm_heatmap_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_votu_rpkm_heatmap.{sampling}.svg",
		assembly_fragmentation_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_assembly_fragmentation.{sampling}.png",
		assembly_fragmentation_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_assembly_fragmentation.{sampling}.svg",
		cluster_coverage_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_coverage.{sampling}.png",
		cluster_coverage_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_coverage.{sampling}.svg",
		host_blast_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_blast_clustermap.{sampling}.png",
		host_blast_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_blast_clustermap.{sampling}.svg",
		percent_covered_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_percent_covered.{sampling}.png",
		percent_covered_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_percent_covered.{sampling}.svg",
		covstats_rpkm_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_covstats_rpkm.{sampling}.png",
		covstats_rpkm_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_covstats_rpkm.{sampling}.svg"
	params:
		samples=SAMPLES,
		sampling="{sampling}",
		results_dir=RESULTS_DIR,
		isolates=ISOLATES,
		metagenome=METAGENOME,
		microbial=MICROBIAL,
		remove_euk=REMOVE_EUK,
		sourmash=SOURMASH
	benchmark:
		dirs_dict["BENCHMARKS"] + "/phage_isolates_summary/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/08_phage_isolates_summary.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_phage_isolates_summary.py.ipynb"


rule subsample_reads:
	input:
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.tot.csv",
		flagstats=expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_filtered_{sample}.tot.txt", sample=SAMPLES),
	output:
		viral_subsampling=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sub_sampling_reads.txt", sample=SAMPLES),
	params:
		samples=SAMPLES,
		mapping_dir=dirs_dict["MAPPING_DIR"],
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="tot",
		key_samples=SAMPLES_key
	benchmark:
		dirs_dict["BENCHMARKS"] + "/subsample_reads/tot.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_subsampling.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_subsampling.py.ipynb"

rule normalise_reads:
	input:
		postqc_txt=dirs_dict["QC_DIR"]+ "/postQC_illumina_report_data/multiqc_fastqc.txt",
		covstats=expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{{sampling}}_covstats.txt", sample=SAMPLES),
		covstats_unique=expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{{sampling}}_unique_covstats.txt", sample=SAMPLES),
	output:
		raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/RPKM_raw_{sampling}.txt",
		norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/RPKM_normalised_{sampling}.txt",
		raw_count_file=dirs_dict["MAPPING_DIR"] + "/counts_raw_{sampling}.txt",
		norm_count_file=dirs_dict["MAPPING_DIR"] + "/counts_normalised_{sampling}.txt",
		coverage_RPKM_file=dirs_dict["MAPPING_DIR"] + "/breadth_coverage_percent_{sampling}.txt",
		coverage_bases_RPKM_file=dirs_dict["MAPPING_DIR"] + "/breadth_coverage_bases_{sampling}.txt",
		filtered_raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/filtered_RPKM_raw_{sampling}.txt",
		filtered_norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/filtered_RPKM_normalised_{sampling}.txt",
		filtered_raw_count_file=dirs_dict["MAPPING_DIR"] + "/filtered_counts_raw_{sampling}.txt",
		filtered_norm_count_file=dirs_dict["MAPPING_DIR"] + "/filtered_counts_normalised_{sampling}.txt",
		filtered_75_raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/filtered_75_RPKM_raw_{sampling}.txt",
		filtered_75_norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/filtered_75_RPKM_normalised_{sampling}.txt",
	params:
		samples=SAMPLES,
		mapping_dir=dirs_dict["MAPPING_DIR"],
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="{sampling}",
		threshold_bases=200,
		threshold_RPKM=0.1,
		reference="",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalise_reads/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_Normalise.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_Normalise.py.ipynb"

rule normalise_reads_reference:
	input:
		postqc_txt=dirs_dict["QC_DIR"]+ "/postQC_illumina_report_data/multiqc_fastqc.txt",
		covstats=expand(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_{{sampling}}_covstats.txt", sample=SAMPLES),
		covstats_unique=expand(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_{{sampling}}_unique_covstats.txt", sample=SAMPLES),
	output:
		raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_RPKM_raw_{sampling}.txt",
		norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_RPKM_normalised_{sampling}.txt",
		raw_count_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_counts_raw_{sampling}.txt",
		norm_count_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_counts_normalised_{sampling}.txt",
		coverage_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_breadth_coverage_percent_{sampling}.txt",
		coverage_bases_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_breadth_coverage_bases_{sampling}.txt",
		filtered_raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_" + REFERENCE + "_RPKM_raw_{sampling}.txt",
		filtered_norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_" + REFERENCE + "_RPKM_normalised_{sampling}.txt",
		filtered_raw_count_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_" + REFERENCE + "_counts_raw_{sampling}.txt",
		filtered_norm_count_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_" + REFERENCE + "_counts_normalised_{sampling}.txt",
		filtered_75_raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_75_" + REFERENCE + "_RPKM_raw_{sampling}.txt",
		filtered_75_norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_75_" + REFERENCE + "_RPKM_normalised_{sampling}.txt",
	params:
		samples=SAMPLES,
		mapping_dir=dirs_dict["MAPPING_DIR"],
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="{sampling}",
		threshold_bases=200,
		threshold_RPKM=0.1,
		reference=REFERENCE,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalise_reads_reference/" + REFERENCE + "/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_Normalise_" + REFERENCE + ".{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_Normalise.py.ipynb"


# Reuse the reference-normalisation notebook and criteria, with separate inputs/outputs.
rule normalise_reads_RefSeq:
	input:
		postqc_txt=dirs_dict["QC_DIR"] + "/postQC_illumina_report_data/multiqc_fastqc.txt",
		covstats=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_covstats.txt", sample=SAMPLES),
		covstats_unique=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_covstats.txt", sample=SAMPLES),
	output:
		raw_RPKM_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_RPKM_raw_tot.txt",
		norm_RPKM_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_RPKM_normalised_tot.txt",
		raw_count_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_counts_raw_tot.txt",
		norm_count_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_counts_normalised_tot.txt",
		coverage_RPKM_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_breadth_coverage_percent_tot.txt",
		coverage_bases_RPKM_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_breadth_coverage_bases_tot.txt",
		mean_coverage_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_mean_depth_tot.txt",
		filtered_raw_RPKM_file=REFSEQ_MAPPING_DIR + "/filtered_RefSeqViral_RPKM_raw_tot.txt",
		filtered_norm_RPKM_file=REFSEQ_MAPPING_DIR + "/filtered_RefSeqViral_RPKM_normalised_tot.txt",
		filtered_raw_count_file=REFSEQ_MAPPING_DIR + "/filtered_RefSeqViral_counts_raw_tot.txt",
		filtered_norm_count_file=REFSEQ_MAPPING_DIR + "/filtered_RefSeqViral_counts_normalised_tot.txt",
		filtered_75_raw_RPKM_file=REFSEQ_MAPPING_DIR + "/filtered_75_RefSeqViral_RPKM_raw_tot.txt",
		filtered_75_norm_RPKM_file=REFSEQ_MAPPING_DIR + "/filtered_75_RefSeqViral_RPKM_normalised_tot.txt",
	params:
		samples=SAMPLES,
		mapping_dir=REFSEQ_MAPPING_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="tot",
		threshold_bases=200,
		threshold_RPKM=0.1,
		reference="RefSeqViral",
		index_label="RefSeq_accession",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalise_reads_RefSeq/tot.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_Normalise_RefSeqViral.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_Normalise.py.ipynb"

rule normalise_reads_all_assembled:
	input:
		postqc_txt=dirs_dict["QC_DIR"] + "/postQC_illumina_report_data/multiqc_fastqc.txt",
		covstats=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_AllAssembled_{sample}_tot_covstats.txt", sample=SAMPLES),
		covstats_unique=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_AllAssembled_{sample}_tot_unique_covstats.txt", sample=SAMPLES),
	output:
		raw_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_RPKM_raw_tot.txt",
		norm_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_RPKM_normalised_tot.txt",
		raw_count_file=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_counts_raw_tot.txt",
		norm_count_file=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_counts_normalised_tot.txt",
		coverage_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_breadth_coverage_percent_tot.txt",
		coverage_bases_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_breadth_coverage_bases_tot.txt",
		mean_coverage_file=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_mean_depth_tot.txt",
		filtered_raw_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_AllAssembled_RPKM_raw_tot.txt",
		filtered_norm_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_AllAssembled_RPKM_normalised_tot.txt",
		filtered_raw_count_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_AllAssembled_counts_raw_tot.txt",
		filtered_norm_count_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_AllAssembled_counts_normalised_tot.txt",
		filtered_75_raw_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_75_AllAssembled_RPKM_raw_tot.txt",
		filtered_75_norm_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_75_AllAssembled_RPKM_normalised_tot.txt",
	params:
		samples=SAMPLES,
		mapping_dir=ALL_ASSEMBLED_MAPPING_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="tot",
		threshold_bases=200,
		threshold_RPKM=0.1,
		reference="AllAssembled",
		index_label="assembled_contig",
		plot_max_points=5000,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalise_reads_all_assembled/tot.tsv"
	resources:
		mem_mb=64000
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_Normalise_AllAssembled.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_Normalise.py.ipynb"

rule all_assembled_mapping_summary:
	input:
		qc=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.tot.csv",
		mapped_pairs=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_mapped_pairs_filtered_AllAssembled_{sample}.tot.txt", sample=SAMPLES),
	output:
		tsv=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_mapping_summary_tot.tsv",
	params:
		samples=SAMPLES,
	message:
		"Summarizing reads mapping to all assembled contigs"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/all_assembled_mapping_summary/tot.tsv"
	threads: 1
	run:
		import csv

		with open(input.qc) as handle:
			cleaned_reads={row["sample"]: int(float(row["bbduk"])) for row in csv.DictReader(handle)}
		with open(output.tsv, "w") as handle:
			writer=csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(["sample", "cleaned_read_pairs", "subsampled_read_pairs", "properly_mapped_pairs", "properly_mapped_percent"])
			for sample, path in zip(params.samples, input.mapped_pairs):
				with open(path) as counts:
					mapped=int(counts.read().strip())
				cleaned=cleaned_reads[sample]
				subsampled=min(2000000, cleaned)
				writer.writerow([sample, cleaned, subsampled, mapped, round(100 * mapped / subsampled, 2) if subsampled else 0])


rule select_all_assembled_top_contigs:
	input:
		fasta=ALL_ASSEMBLED_DIR + "/all_assembled_contigs_derreplicated_rep_seq.tot.fasta",
		clusters=ALL_ASSEMBLED_DIR + "/all_assembled_contigs_derreplicated_cluster.tot.tsv",
		provenance=ALL_ASSEMBLED_DIR + "/all_assembled_contigs_provenance.tot.tsv",
		rpkm=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_RPKM_raw_tot.txt",
		counts=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_counts_raw_tot.txt",
		breadth=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_breadth_coverage_percent_tot.txt",
		depth=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_mean_depth_tot.txt",
		unique=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_AllAssembled_{sample}_tot_unique_covstats.txt", sample=SAMPLES),
	output:
		fasta=ALL_ASSEMBLED_TOP_PREFIX + ".fasta",
		ranking=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_top_contigs_tot.tsv",
		membership=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_top_contigs_membership_tot.tsv",
	params:
		samples=SAMPLES,
		top_n=int(config.get("all_assembled_top_n", 100)),
		negative_control=str(config.get("negative_control", "")).strip(),
	message:
		"Selecting the most abundant assembled contigs by mean raw RPKM"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/select_all_assembled_top_contigs/tot.tsv"
	threads: 1
	resources:
		mem_mb=16000,
	run:
		import pandas as pd

		rpkm=pd.read_csv(input.rpkm, index_col=0).reindex(columns=params.samples).fillna(0)
		ranking_samples=[sample for sample in params.samples if sample != params.negative_control]
		ranked=rpkm[ranking_samples]
		top=pd.DataFrame(index=rpkm.index)
		top["mean_RPKM_raw"]=ranked.mean(axis=1).fillna(0)
		top["max_RPKM_raw"]=ranked.max(axis=1).fillna(0)
		top["highest_abundance_sample"]=ranked.idxmax(axis=1) if ranking_samples else "not reported"
		top["contig_id"]=top.index
		top=top.loc[top["mean_RPKM_raw"] > 0].sort_values(
			["mean_RPKM_raw", "contig_id"], ascending=[False, True], kind="stable").head(max(0, params.top_n))
		top.insert(0, "rank", range(1, len(top) + 1))
		provenance=pd.read_csv(input.provenance, sep="\t")
		origins=provenance.set_index("contig_id")
		for field in ["sample", "assembler", "original_id", "length_bp"]:
			top["assembly_sample" if field == "sample" else field]=origins[field].reindex(top.index)
		clusters=pd.read_csv(input.clusters, sep="\t", header=None, names=["representative_id", "member_id"])
		clusters=clusters.drop_duplicates()
		top["cluster_size"]=clusters.groupby("representative_id").size().reindex(top.index).fillna(1).astype(int)
		for path, prefix in [(input.rpkm, "RPKM_raw"), (input.counts, "mapped_reads"),
				(input.breadth, "breadth_percent"), (input.depth, "mean_depth")]:
			matrix=pd.read_csv(path, index_col=0).reindex(index=top.index, columns=params.samples).fillna(0)
			for sample in params.samples:
				top[prefix + "_" + sample]=matrix[sample]
		for sample, path in zip(params.samples, input.unique):
			unique=pd.read_csv(path, sep="\t", usecols=[0, 4]).set_index("Contig")
			top["unique_mapped_reads_" + sample]=unique.iloc[:, 0].reindex(top.index).fillna(0)
			top["unique_mapping_ratio_" + sample]=(top["unique_mapped_reads_" + sample] /
				top["mapped_reads_" + sample].replace(0, float("nan"))).fillna(0)
		membership=clusters.loc[clusters["representative_id"].isin(top.index)].merge(
			provenance.rename(columns={"contig_id": "member_id", "sample": "assembly_sample"}), on="member_id", how="left")
		membership["rank"]=membership["representative_id"].map(top["rank"])
		membership=membership.sort_values(["rank", "member_id"], kind="stable")
		membership.to_csv(output.membership, sep="\t", index=False)

		def records(path):
			name, chunks=None, []
			with open(path) as handle:
				for line in handle:
					if line.startswith(">"):
						if name is not None:
							yield name, "".join(chunks)
						name, chunks=line[1:].split()[0], []
					else:
						chunks.append(line.strip())
			if name is not None:
				yield name, "".join(chunks)

		selected={name: sequence for name, sequence in records(input.fasta) if name in top.index}
		with open(output.fasta, "w") as handle:
			for name in top.index:
				sequence=selected[name]
				handle.write(f">{name}\n{sequence}\n")
				top.loc[name, "gc_percent"]=100 * sum(sequence.upper().count(base) for base in "GC") / len(sequence)
		top.to_csv(output.ranking, sep="\t", index=False)

rule collect_all_assembled_top_evidence:
	input:
		ranking=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_top_contigs_tot.tsv",
		genomad_dna=expand(dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_tot", sample=SAMPLES),
		genomad_rna=expand(RNA_DIR + "/{sample}/{assembler}_genomad", sample=SAMPLES, assembler=RNA_ASSEMBLERS) if RNA_MODE else [],
		circularity=expand(dirs_dict["VIRAL_DIR"] + "/{sample}_{assembler}_circularity.tot.tsv", sample=SAMPLES, assembler=["spades"] + (RNA_ASSEMBLERS if RNA_MODE else [])),
	output:
		tsv=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_top_contigs_existing_annotations_tot.tsv",
	params:
		viral_dir=dirs_dict["VIRAL_DIR"],
		rna_dir=RNA_DIR,
	message:
		"Collecting original-assembly geNomad and terminal-repeat evidence without reclassifying"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/collect_all_assembled_top_evidence/tot.tsv"
	threads: 1
	run:
		import csv
		from pathlib import Path

		def rows(path):
			if path.is_file() and path.stat().st_size:
				with path.open() as handle:
					return list(csv.DictReader(handle, delimiter="\t"))
			return []

		circles={}
		for path in input.circularity:
			for row in rows(Path(path)):
				circles[(row["sample"], row["assembler"], row["contig_id"])]=row
		circle_fields=["dtr_bp", "itr_bp", "dtr_sequence", "dtr_left_start", "dtr_left_end", "dtr_right_start", "dtr_right_end",
			"itr_left_sequence", "itr_right_sequence", "itr_left_start", "itr_left_end", "itr_right_start", "itr_right_end",
			"terminal_repeat_type", "repeat_warning"]
		fields=["contig_id", "genomad_assessed", "genomad_classification", "genomad_virus_score", "genomad_plasmid_score",
			"genomad_taxonomy", "genomad_topology", "genomad_evidence_scope", "genomad_reported_sequence", "genomad_source"] + circle_fields
		cache={}
		with open(output.tsv, "w") as handle:
			writer=csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
			writer.writeheader()
			for origin in rows(Path(input.ranking)):
				name, sample, assembler=origin["original_id"], origin["assembly_sample"], origin["assembler"]
				if assembler == "spades":
					folder=Path(params.viral_dir) / (sample + "_geNomad_tot")
					prefix=sample + "_spades_filtered_scaffolds.tot"
				else:
					folder=Path(params.rna_dir) / sample / (assembler + "_genomad")
					prefix=assembler
				if folder not in cache:
					summary=folder / (prefix + "_summary")
					virus=summary / (prefix + "_virus_summary.tsv")
					plasmid=summary / (prefix + "_plasmid_summary.tsv")
					cache[folder]=(virus.is_file(), rows(virus), rows(plasmid))
				assessed, viruses, plasmids=cache[folder]
				exact=[row for row in viruses if row["seq_name"] == name]
				regional=[row for row in viruses if row["seq_name"].split("|provirus_", 1)[0] == name and row["seq_name"] != name]
				plasmids=[row for row in plasmids if row["seq_name"] == name]
				hits=exact or regional or plasmids
				result=dict.fromkeys(fields, "not reported")
				result.update(contig_id=origin["contig_id"], genomad_assessed="yes" if assessed else "not reported", genomad_source=str(folder))
				if hits:
					result["genomad_classification"]="virus" if exact or regional else "plasmid"
					result["genomad_evidence_scope"]="provirus region" if regional and not exact else "whole contig"
					for column, source in [("genomad_virus_score", "virus_score"), ("genomad_plasmid_score", "plasmid_score"),
							("genomad_taxonomy", "taxonomy"), ("genomad_topology", "topology"), ("genomad_reported_sequence", "seq_name")]:
						values=[row[source] for row in hits if row.get(source)]
						result[column]=" | ".join(values) if values else "not reported"
				circle=circles.get((sample, assembler, name), {})
				for column in circle_fields:
					result[column]=circle.get(column) or "not reported"
				writer.writerow(result)

rule all_assembled_top_metadata:
	input:
		ranking=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_top_contigs_tot.tsv",
		existing=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_top_contigs_existing_annotations_tot.tsv",
		virsorter=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + ALL_ASSEMBLED_TOP_NAME + "_tot/final-viral-score.tsv",
		vibrant_quality=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + ALL_ASSEMBLED_TOP_NAME + "_positive_quality.tot.csv",
		vibrant_summary=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + ALL_ASSEMBLED_TOP_NAME + "_summary_results.tot.csv",
		checkv=ALL_ASSEMBLED_TOP_PREFIX + "_checkV/quality_summary.tsv",
		refseq=dirs_dict["ANNOTATION"] + "/blast_output_ViralRefSeq_" + ALL_ASSEMBLED_TOP_NAME + ".tot.csv",
		metavr=dirs_dict["ANNOTATION"] + "/blast_output_METAVR_" + ALL_ASSEMBLED_TOP_NAME + ".tot.csv",
		metavr_metadata=dirs_dict["PLOTS_DIR"] + "/" + ALL_ASSEMBLED_TOP_NAME + "_METAVR.tot/METAVR_main_table_for_hits.tsv",
		pharokka=ALL_ASSEMBLED_TOP_PREFIX + "_pharokka",
	output:
		metadata=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_top_contigs_metadata_tot.tsv",
		pharokka_cds=ALL_ASSEMBLED_MAPPING_DIR + "/AllAssembled_top_contigs_pharokka_cds_tot.tsv",
	message:
		"Combining abundance, classification and annotation of the top assembled contigs"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/all_assembled_top_metadata/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	run:
		import pandas as pd
		from pathlib import Path

		def table(path, **kwargs):
			path=Path(path)
			return pd.read_csv(path, sep="\t", **kwargs) if path.is_file() and path.stat().st_size else pd.DataFrame()

		metadata=table(input.ranking).set_index("contig_id")
		existing=table(input.existing).set_index("contig_id")
		metadata=metadata.join(existing)
		for path, key, prefix in [(input.virsorter, "seqname", "VirSorter2"), (input.checkv, "contig_id", "CheckV")]:
			frame=table(path)
			if key in frame:
				frame=frame.drop_duplicates(key).set_index(key).add_prefix(prefix + "_")
				metadata=metadata.join(frame.reindex(metadata.index))

		def original_name(name):
			if name in metadata.index:
				return name
			for separator in ["_fragment_", "|provirus_", "_provirus_"]:
				if name.split(separator, 1)[0] in metadata.index:
					return name.split(separator, 1)[0]
			return name

		for path, prefix in [(input.vibrant_quality, "VIBRANT"), (input.vibrant_summary, "VIBRANT_summary")]:
			frame=table(path)
			if "scaffold" in frame:
				frame["contig_id"]=frame["scaffold"].map(original_name)
				frame["evidence_scope"]=frame.apply(lambda row: "whole contig" if row["scaffold"] == row["contig_id"] else "region", axis=1)
				frame=frame.groupby("contig_id", sort=False).agg(lambda values: " | ".join(dict.fromkeys(values.dropna().astype(str))))
				frame.columns=[prefix + "_" + column.replace(" ", "_") for column in frame.columns]
				metadata=metadata.join(frame.reindex(metadata.index))

		blast_columns=["qseqid", "sseqid", "description", "qstart", "qend", "qlen", "slen", "qcovs", "evalue", "alignment_length", "pident"]
		for path, prefix in [(input.refseq, "RefSeq"), (input.metavr, "METAVR")]:
			frame=table(path, header=None, names=blast_columns)
			frame=frame.loc[frame["qseqid"] != "qseqid"] if "qseqid" in frame else pd.DataFrame(columns=blast_columns)
			for field in blast_columns[3:]:
				frame[field]=pd.to_numeric(frame[field], errors="coerce")
			frame=frame.sort_values(["evalue", "qcovs", "pident", "alignment_length", "sseqid"], ascending=[True, False, False, False, True], kind="stable")
			frame=frame.drop_duplicates("qseqid").set_index("qseqid").add_prefix(prefix + "_")
			metadata=metadata.join(frame.reindex(metadata.index))
		metavr=table(input.metavr_metadata)
		if "uvig" in metavr:
			metavr=metavr.drop_duplicates("uvig").set_index("uvig")
			for field in metavr.columns:
				metadata["METAVR_" + field]=metadata["METAVR_sseqid"].fillna("").str.split("|").str[0].map(metavr[field])

		cds=table(Path(input.pharokka) / "pharokka_cds_final_merged_output.tsv")
		if "contig" in cds:
			cds=cds.loc[cds["contig"].isin(metadata.index)].copy()
			metadata["Pharokka_CDS_count"]=cds.groupby("contig").size().reindex(metadata.index, fill_value=0)
			if "annot" in cds:
				known=cds.loc[cds["annot"].notna() & ~cds["annot"].str.contains("hypothetical|unknown", case=False, na=False)]
				metadata["Pharokka_annotated_CDS_count"]=known.groupby("contig").size().reindex(metadata.index, fill_value=0)
			if "category" in cds:
				metadata["Pharokka_function_categories"]=cds.groupby("contig")["category"].agg(lambda values: "; ".join(sorted(set(values.dropna()))))
		else:
			cds=pd.DataFrame(columns=["contig", "gene", "annot", "category"])
			metadata["Pharokka_CDS_count"]="not reported"
			metadata["Pharokka_annotated_CDS_count"]="not reported"
			metadata["Pharokka_function_categories"]="not reported"
		cds.to_csv(output.pharokka_cds, sep="\t", index=False)
		metadata=metadata.replace("", "not reported").fillna("not reported")
		metadata.reset_index().to_csv(output.metadata, sep="\t", index=False)


if MAP_TO_REFSEQ:
	rule refseq_detection:
		input:
			fasta=config["RefSeqViral_db"],
			covstats=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_covstats.txt", sample=SAMPLES),
			unique_covstats=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_covstats.txt", sample=SAMPLES),
			basecov=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_basecov.txt", sample=SAMPLES),
			unique_basecov=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_basecov.txt", sample=SAMPLES),
			rpkm=REFSEQ_MAPPING_DIR + "/RefSeqViral_RPKM_raw_tot.txt",
			counts=REFSEQ_MAPPING_DIR + "/RefSeqViral_counts_raw_tot.txt",
			breadth=REFSEQ_MAPPING_DIR + "/RefSeqViral_breadth_coverage_percent_tot.txt",
			depth=REFSEQ_MAPPING_DIR + "/RefSeqViral_mean_depth_tot.txt",
		output:
			long=REFSEQ_MAPPING_DIR + "/RefSeqViral_detection_long_tot.tsv",
			matrix=REFSEQ_MAPPING_DIR + "/RefSeqViral_detection_matrix_tot.tsv",
			metadata=REFSEQ_MAPPING_DIR + "/RefSeqViral_detection_metadata_tot.tsv",
			high=REFSEQ_MAPPING_DIR + "/RefSeqViral_detected_high_confidence_tot.tsv",
			candidates=REFSEQ_MAPPING_DIR + "/RefSeqViral_detected_candidates_tot.tsv",
			contaminants=REFSEQ_MAPPING_DIR + "/RefSeqViral_detected_contaminants_tot.tsv",
			figures=directory(dirs_dict["PLOTS_DIR"] + "/07_RefSeq_detection"),
		params:
			samples=SAMPLES,
			metadata_cache=REFSEQ_MAPPING_DIR + "/RefSeqViral_NCBI_metadata_cache.sqlite",
			metadata_refresh=config_bool("refseq_metadata_refresh", False),
			negative_control=str(config.get("negative_control", "")).strip(),
			high_min_length_bp=int(config.get("refseq_high_min_length_bp", 6667)),
			high_min_covered_bases=int(config.get("refseq_high_min_covered_bases", 5000)),
			high_min_breadth_percent=float(config.get("refseq_high_min_breadth_percent", 75)),
			high_min_unique_reads=int(config.get("refseq_high_min_unique_reads", 5)),
			candidate_min_unique_reads=int(config.get("refseq_candidate_min_unique_reads", 10)),
			candidate_min_breadth_percent=float(config.get("refseq_candidate_min_breadth_percent", 25)),
			ambiguous_min_reads=int(config.get("refseq_ambiguous_min_reads", 5)),
			ambiguous_max_unique_ratio=float(config.get("refseq_ambiguous_max_unique_ratio", 0.2)),
			negative_control_min_reads=int(config.get("refseq_negative_control_min_reads", 5)),
			negative_control_max_enrichment=float(config.get("refseq_negative_control_max_enrichment", 2)),
			enrichment_pseudocount=float(config.get("refseq_enrichment_pseudocount", 0.01)),
			plot_max_accessions=int(config.get("refseq_plot_max_accessions", 60)),
		message:
			"Classifying RefSeq viral read-mapping evidence"
		benchmark:
			dirs_dict["BENCHMARKS"] + "/refseq_detection/tot.tsv"
		threads: 1
		resources:
			mem_mb=8000,
		log:
			notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_RefSeq_detection.py.ipynb"
		notebook:
			dirs_dict["RAW_NOTEBOOKS"] + "/07_RefSeq_detection.py.ipynb"

def input_QC_long_only_nanopore_pre(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanostats_preQC.html", sample_nanopore=NANOPORE_SAMPLES))
	return(input_list)

def input_QC_long_only_nanopore_post(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanostats_postQC.html", sample_nanopore=NANOPORE_SAMPLES))
	return(input_list)

def input_QC_long_only_pacbio_pre(wildcards):
	input_list=[]
	if PACBIO:
		input_list.extend(expand(dirs_dict["QC_DIR"] + "/{sample_pacbio}_pacbio_nanostats_preQC.html", sample_pacbio=PACBIO_SAMPLES))
	return(input_list)

def input_QC_long_only_pacbio_post(wildcards):
	input_list=[]
	if PACBIO:
		input_list.extend(expand(dirs_dict["QC_DIR"] + "/{sample_pacbio}_pacbio_nanostats_postQC_{sampling}.html", sample_pacbio=PACBIO_SAMPLES, sampling=wildcards.sampling))
	return(input_list)


rule QC_long_only_parsing:
	input:
		nanopore_pre=input_QC_long_only_nanopore_pre,
		nanopore_post=input_QC_long_only_nanopore_post,
		pacbio_pre=input_QC_long_only_pacbio_pre,
		pacbio_post=input_QC_long_only_pacbio_post
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_summary.{sampling}.html",
		read_count_png=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_count.{sampling}.png",
		read_count_svg=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_count.{sampling}.svg",
		total_bases_png=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_total_bases.{sampling}.png",
		total_bases_svg=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_total_bases.{sampling}.svg",
		read_length_png=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_length.{sampling}.png",
		read_length_svg=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_length.{sampling}.svg",
		read_quality_png=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_quality.{sampling}.png",
		read_quality_svg=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_quality.{sampling}.svg"
	params:
		sampling="{sampling}",
		samples_nanopore=NANOPORE_SAMPLES,
		samples_pacbio=PACBIO_SAMPLES
	benchmark:
		dirs_dict["BENCHMARKS"] + "/QC_long_only_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/01_QC_long_only.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/01_QC_long_only.py.ipynb"

rule assembly_long_only_parsing:
	input:
		quast_report_dir=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}"
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_summary.{sampling}.html",
		contig_number_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_contig_number.{sampling}.png",
		contig_number_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_contig_number.{sampling}.svg",
		total_length_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_total_length.{sampling}.png",
		total_length_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_total_length.{sampling}.svg",
		n50_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_N50.{sampling}.png",
		n50_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_N50.{sampling}.svg"
	params:
		input_quast_report=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}/transposed_report.tsv",
		sampling="{sampling}"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assembly_long_only_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/03_assembly_long_only.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/03_assembly_long_only.py.ipynb"

def input_bacterial_results_checkm(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if ISOLATES:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}", sample=SAMPLES, sampling=wildcards.sampling))
	return(input_list)

def input_bacterial_results_sourmash(wildcards):
	input_list=[]
	if NANOPORE & NANOPORE_ONLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore.classifications.csv", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if NANOPORE & (not NANOPORE_ONLY) & PAIRED:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_hybrid.classifications.csv", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO & PACBIO_ONLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio.classifications.csv", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if PACBIO & PACBIO_HYBRID:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_hybrid.classifications.csv", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if ISOLATES:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}.classifications.csv", sample=SAMPLES, sampling=wildcards.sampling))
	return(input_list)

def input_bacterial_results_coverage(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["MAPPING_DIR"] + "/{sample}_{sampling}_long_read_contig_coverage.tsv", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO:
		input_list.extend(expand(dirs_dict["MAPPING_DIR"] + "/{sample}_{sampling}_long_read_contig_coverage.tsv", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	return(input_list)


def input_bacterial_results_genomad(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"] + "/{sample}_long_geNomad_{sampling}", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"] + "/{sample}_long_geNomad_{sampling}", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if ISOLATES:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}", sample=SAMPLES, sampling=wildcards.sampling))
	return(input_list)


rule selectMetaVRMetadata:
	input:
		blast=lambda wc: dirs_dict["ANNOTATION"] + "/blast_output_METAVR_" + ("filtered_" + REPRESENTATIVE_CONTIGS_BASE if wc.report == "08_METAVR_analysis" else wc.report.removesuffix("_METAVR")) + "." + wc.sampling + ".csv",
		metadata=os.path.join(config["METAVR_db"], "METAVR_main_table.parquet"),
	output:
		selected=dirs_dict["PLOTS_DIR"] + "/{report}.{sampling}/METAVR_main_table_for_hits.tsv",
	message:
		"Selecting MetaVR metadata for the matched viral genomes"
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/selectMetaVRMetadata/report={report}__sampling={sampling}.tsv"
	threads: 1
	wildcard_constraints:
		report="[^/]+",
		sampling="tot|sub"
	resources:
		mem_mb=8000,
	shell:
		"""
		set -euo pipefail
		mkdir -p $(dirname {output.selected:q})
		python - {input.blast:q} {input.metadata:q} {output.selected:q} <<'PYCODE'
import csv
import sys
from fastparquet import ParquetFile

columns = ["uvig", "taxon_oid", "quality", "completeness", "genome_type",
           "ictv_taxonomy", "ictv_taxonomy_method", "host_taxonomy", "host_taxonomy_method"]
with open(sys.argv[1], newline="") as handle:
    subjects = set()
    for row in csv.reader(handle, delimiter="\t"):
        if len(row) > 1 and row[1] != "sseqid":
            subjects.add(row[1].split("|")[0])
with open(sys.argv[3], "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(columns)
    if subjects:
        parquet = ParquetFile(sys.argv[2])
        for group in parquet.iter_row_groups(columns=columns):
            selected = group.loc[group["uvig"].isin(subjects), columns]
            writer.writerows(selected.itertuples(index=False, name=None))
PYCODE
		"""


rule METAVR_analysis:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		blast=dirs_dict["ANNOTATION"] + "/blast_output_METAVR_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.csv",
		checkv=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		uvig_metadata=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/METAVR_main_table_for_hits.tsv",
		source_metadata=os.path.join(config["METAVR_db"], "IMG_full_metadata.tsv.gz"),
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}.html",
		hit_pairs=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/hit_pairs.tsv",
		query_summary=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/vOTU_summary.tsv",
		sample_summary=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/sample_summary.tsv",
		environment_counts=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/environment_counts.tsv",
		viral_taxonomy_counts=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/viral_taxonomy_counts.tsv",
		host_taxonomy_counts=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/host_taxonomy_counts.tsv",
		host_method_counts=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/host_method_counts.tsv",
		metadata_subset=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/METAVR_metadata_for_hits.tsv",
		metadata_manifest=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/METAVR_metadata_for_hits.manifest.json",
		run_summary=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/run_summary.tsv",
		provenance=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/provenance.json",
		composition_png=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/coverage_composition.png",
		composition_svg=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/coverage_composition.svg",
		samples_png=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/sample_composition.png",
		samples_svg=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/sample_composition.svg",
	params:
		samples=SAMPLES,
		sampling="{sampling}",
		metadata_chunksize=250000,
	message:
		"Summarizing MetaVR hit environments, viral taxonomy and recorded host taxonomy"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/METAVR_analysis/sampling={sampling}.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/08_METAVR_analysis.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_METAVR_analysis.py.ipynb"


rule bacterial_results_parsing:
	input:
		checkm=input_bacterial_results_checkm,
		sourmash=input_bacterial_results_sourmash,
		gtdbtk=dirs_dict["ASSEMBLY_DIR"] + "/assembly_bacteria_GTDB-Tk_{sampling}",
		quast=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}/transposed_report.tsv",
		coverage=input_bacterial_results_coverage,
		genomad=input_bacterial_results_genomad,
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_summary.{sampling}.html",
		summary_csv=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_summary.{sampling}.csv",
		checkm_png=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_checkm.{sampling}.png",
		checkm_svg=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_checkm.{sampling}.svg",
		taxonomy_png=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_taxonomy.{sampling}.png",
		taxonomy_svg=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_taxonomy.{sampling}.svg",
		quality_taxonomy_png=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_quality_taxonomy.{sampling}.png",
		quality_taxonomy_svg=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_quality_taxonomy.{sampling}.svg"
	params:
		sampling="{sampling}",
		long_assembler=LONG_ASSEMBLER,
		long_assembler_pacbio=LONG_ASSEMBLER_PACBIO,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/bacterial_results_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/06_bacterial_results.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/06_bacterial_results.py.ipynb"
