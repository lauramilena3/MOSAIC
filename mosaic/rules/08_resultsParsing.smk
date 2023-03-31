# rule plot_kmer:
# 	input:
# 		histograms=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.{{sampling}}.csv", sample=SAMPLES),
# 	output:
# 		plot=(dirs_dict["CLEAN_DATA_DIR"] + "/kmer_rarefraction_plot.{sampling}.png"),
# 		svg=(dirs_dict["CLEAN_DATA_DIR"] + "/kmer_rarefraction_plot.{sampling}.svg"),

# 	message:
# 		"Plot unique reads with BBtools"
# 	threads: 1
# 	run:
# 		import pandas as pd
# 		import seaborn as sns; sns.set()
# 		import matplotlib.pyplot as plt

# 		plt.figure(figsize=(12,12))
# 		sns.set(font_scale=2)
# 		sns.set_style("whitegrid")

# 		read_max=0

# 		for h in input.histograms:
# 			df=pd.read_csv(h, sep="\t")
# 			df.columns=["count", "percent", "c", "d", "e", "f", "g", "h", "i", "j"]
# 			df=df[["count", "percent"]]
# 			ax = sns.lineplot(x="count", y="percent", data=df,err_style='band', label=h.split("/")[-1].split("_kmer")[0])
# 			read_max=max(read_max,df["count"].max())

# 		ax.set(ylim=(0, 100))
# 		ax.set(xlim=(0, read_max*1.2))

# 		ax.set_xlabel("Read count",fontsize=20)
# 		ax.set_ylabel("New k-mers (%)",fontsize=20)
# 		ax.figure.savefig(output.plot)
# 		ax.figure.savefig(output.svg, format="svg")

rule copy_notebooks:
	input:
		notebook="notebooks/{notebook}.py.ipynb"
	output:
		notebooks=dirs_dict["NOTEBOOKS_DIR"] +"{notebook}.py.ipynb"
	shell:
		"""
		cp {input.notebook} {output.notebook}
		"""

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
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/01_QC.py.ipynb",
		histograms=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.{{sampling}}.csv", sample=SAMPLES),
		preqc_txt=dirs_dict["QC_DIR"]+ "/preQC_illumina_report_data/multiqc_fastqc.txt",
		postqc_txt=dirs_dict["QC_DIR"]+ "/postQC_illumina_report_data/multiqc_fastqc.txt",
	output:
		kmer_png=(dirs_dict["PLOTS_DIR"] + "/kmer_rarefraction_plot.{sampling}.png"),
		kmer_svg=(dirs_dict["PLOTS_DIR"] + "/kmer_rarefraction_plot.{sampling}.svg"),
		kmer_fit_png=(dirs_dict["PLOTS_DIR"] + "/kmer_rarefraction_plot_fitted.{sampling}.png"),
		kmer_fit_svg=(dirs_dict["PLOTS_DIR"] + "/kmer_rarefraction_plot_fitted.{sampling}.svg"),
		kmer_fit_html=(dirs_dict["PLOTS_DIR"] + "/kmer_rarefraction_fitted.{sampling}.html"),
		qc_summary_html=(dirs_dict["PLOTS_DIR"] + "/post_qc_read_summary.{sampling}.html"),
		percentage_kept_reads_png=(dirs_dict["PLOTS_DIR"] + "/percentage_kept_reads.{sampling}.png"),
		percentage_kept_reads_svg=(dirs_dict["PLOTS_DIR"] + "/percentage_kept_reads.{sampling}.svg"),
		percentage_kept_Mbp_png=(dirs_dict["PLOTS_DIR"] + "/percentage_kept_Mbp.{sampling}.png"),
		percentage_kept_Mbp_svg=(dirs_dict["PLOTS_DIR"] + "/percentage_kept_Mbp.{sampling}.svg"),
		step_qc_reads_html=(dirs_dict["PLOTS_DIR"] + "/multistep_qc_report.{sampling}.html"),
		steps_qc_reads_png=(dirs_dict["PLOTS_DIR"] + "/qc_bystep_counts.{sampling}.png"),
		steps_qc_reads_svg=(dirs_dict["PLOTS_DIR"] + "/qc_bystep_counts.{sampling}.svg"),
		steps_qc_percentage_png=(dirs_dict["PLOTS_DIR"] + "/qc_bystep_percentage.{sampling}.png"),
		steps_qc_percentage_svg=(dirs_dict["PLOTS_DIR"] + "/qc_bystep_percentage.{sampling}.svg"),
		supperdedupper_html=(dirs_dict["PLOTS_DIR"] + "/superdedupper_PCR.{sampling}.html"),
		supperdedupper_png=(dirs_dict["PLOTS_DIR"] + "/superdedupper_PCR.{sampling}.png"),
		supperdedupper_svg=(dirs_dict["PLOTS_DIR"] + "/superdedupper_PCR.{sampling}.svg"),
	params:
		results_dir=RESULTS_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		samples=SAMPLES,
	notebook:
		dirs_dict["NOTEBOOKS_DIR"] + "/01_QC.py.ipynb"



