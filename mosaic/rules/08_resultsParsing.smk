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
		histograms=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.{{sampling}}.csv", sample=SAMPLES),
	output:
		kmer_png=(dirs_dict["PLOTS_DIR"] + "/kmer_rarefraction_plot.{sampling}.png"),
		kmer_svg=(dirs_dict["PLOTS_DIR"] + "/kmer_rarefraction_plot.{sampling}.svg"),
	params:
		results_dir=RESULTS_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		samples=SAMPLES,

	notebook:
		dirs_dict["NOTEBOOKS_DIR"] + "/01_QC.py.ipynb"



