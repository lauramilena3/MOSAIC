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
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
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

rule preQualityCheckPacbio:
	input:
		pacbio=RAW_DATA_DIR + "/{sample_pacbio}_" + str(config['pacbio_tag']) + ".fastq.gz"
	output:
		html=dirs_dict["QC_DIR"] + "/{sample_pacbio}_pacbio_report_preQC.html",
		stats=dirs_dict["QC_DIR"] + "/{sample_pacbio}_pacbio_nanostats_preQC.html"
	message:
		"Performing PacBio HiFi pre-QC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/preQualityCheckPacbio/{sample_pacbio}.tsv"
	threads: 4
	shell:
		"""
		NanoPlot --fastq {input.pacbio} --outdir {wildcards.sample_pacbio}_pacbio_preQC_temp --threads {threads}
		mv {wildcards.sample_pacbio}_pacbio_preQC_temp/NanoPlot-report.html {output.html}
		mv {wildcards.sample_pacbio}_pacbio_preQC_temp/NanoStats.txt {output.stats}
		rm -rf {wildcards.sample_pacbio}_pacbio_preQC_temp
		"""

rule postQualityCheckPacbio:
	input:
		pacbio=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_pacbio}_pacbio_clean.{sampling}.fastq.gz"
	output:
		html=dirs_dict["QC_DIR"] + "/{sample_pacbio}_pacbio_report_postQC_{sampling}.html",
		stats=dirs_dict["QC_DIR"] + "/{sample_pacbio}_pacbio_nanostats_postQC_{sampling}.html"
	message:
		"Performing PacBio HiFi post-QC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/postQualityCheckPacbio/{sample_pacbio}_{sampling}.tsv"
	threads: 4
	shell:
		"""
		NanoPlot --fastq {input.pacbio} --outdir {wildcards.sample_pacbio}_pacbio_postQC_{wildcards.sampling}_temp --threads {threads}
		mv {wildcards.sample_pacbio}_pacbio_postQC_{wildcards.sampling}_temp/NanoPlot-report.html {output.html}
		mv {wildcards.sample_pacbio}_pacbio_postQC_{wildcards.sampling}_temp/NanoStats.txt {output.stats}
		rm -rf {wildcards.sample_pacbio}_pacbio_postQC_{wildcards.sampling}_temp
		"""