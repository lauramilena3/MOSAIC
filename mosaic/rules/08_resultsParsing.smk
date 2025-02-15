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
		read_count_euk_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_noEuk.tot_read_count.txt", sample=SAMPLES),
		read_count_euk_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_noEuk.tot_read_count.txt", sample=SAMPLES),
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
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/03_assembly_short.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/03_assembly_short.py.ipynb"

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
		assembly_flagstats=input_assembly_flagstats
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
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_mapping_statistics_{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_mapping_statistics.py.ipynb"

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
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_Normalise_" + REFERENCE + ".{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_Normalise.py.ipynb"