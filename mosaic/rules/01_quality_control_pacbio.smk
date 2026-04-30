rule cleanPacbioReads:
	input:
		pacbio=RAW_DATA_DIR + "/{sample_pacbio}_" + str(config['pacbio_tag']) + ".fastq.gz"
	output:
		pacbio=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_pacbio}_pacbio_clean.{sampling}.fastq.gz"
	params:
		min_length=config['pacbio_min_length'],
		min_quality=config['pacbio_min_quality']
	message:
		"Filtering PacBio HiFi reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/pacbio.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/cleanPacbioReads/{sample_pacbio}_{sampling}.tsv"
	threads: 4
	shell:
		"""
		seqkit seq \
			-j {threads} \
			-m {params.min_length} \
			-Q {params.min_quality} \
			{input.pacbio} \
			-o {output.pacbio}
		"""

rule nanoPlotPacbio:
	input:
		pacbio=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_pacbio}_pacbio_clean.{sampling}.fastq.gz"
	output:
		html=dirs_dict["QC_DIR"] + "/pacbio_{sample_pacbio}_{sampling}/NanoPlot-report.html",
		stats=dirs_dict["QC_DIR"] + "/pacbio_{sample_pacbio}_{sampling}/NanoStats.txt"
	params:
		outdir=dirs_dict["QC_DIR"] + "/pacbio_{sample_pacbio}_{sampling}"
	message:
		"Running NanoPlot on PacBio HiFi reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/pacbio.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/nanoPlotPacbio/{sample_pacbio}_{sampling}.tsv"
	threads: 4
	shell:
		"""
		NanoPlot \
			--fastq {input.pacbio} \
			--outdir {params.outdir} \
			--threads {threads}
		"""