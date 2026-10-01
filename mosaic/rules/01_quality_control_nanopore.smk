rule preQualityCheckNanopore:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"]+"/{sample_nanopore}_nanopore.fastq.gz",
	output:
		nanoqc_dir=temp(directory(dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanoplot")),
		nanoqc=dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_preQC.html",
	message:
		"Performing nanoQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/preQualityCheckNanopore/sample_nanopore={sample_nanopore}.tsv"
#	threads: 1
	shell:
		"""
		nanoQC -o {output.nanoqc_dir} {input.raw_fastq}
		mv {output.nanoqc_dir}/nanoQC.html {output.nanoqc}
		"""

rule remove_adapters_quality_nanopore:
	input:
		raw_data=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore.fastq.gz",
	output:
		porechopped=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_porechopped.fastq.gz"),
		trimmed_data=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq.gz"),
	params:
		headcrop=config['nanopore_headcrop'],
		tailcrop=config['nanopore_tailcrop'],
		quality=config['nanopore_quality'],
		minlen=config['nanopore_minlen'],
		maxlen=config['nanopore_maxlen'],
	message:
		"Trimming Nanopore Adapters with Porechop"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/remove_adapters_quality_nanopore/sample_nanopore={sample_nanopore}.tsv"
	threads: 16
	shell:
		"""
		porechop -i {input.raw_data} -o {output.porechopped} --threads {threads} --discard_middle
		gunzip -c {output.porechopped} | NanoFilt -q {params.quality} -l {params.minlen} --maxlength {params.maxlen} \
				--headcrop {params.headcrop} --tailcrop {params.tailcrop} | gzip > {output.trimmed_data}
		"""

rule remove_contaminants_nanopore:
	input:
		trimmed_data=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq.gz",
		contaminants_fasta=expand(dirs_dict["CONTAMINANTS_DIR"] +"/{contaminants}.fasta",contaminants=CONTAMINANTS),
	output:
		fastq=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.fastq.gz"),
		size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean_read_count.tot.txt",
		phix_contaminants_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/{sample_nanopore}_nanopore_contaminants.fasta",
	params:
		has_contaminants=bool(CONTAMINANTS),
	message:
		"Removing configured contaminants or passing Nanopore reads through unchanged"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/remove_contaminants_nanopore/sample_nanopore={sample_nanopore}.tsv"
	threads: 2
	shell:
		"""
		if [ "{params.has_contaminants}" = "True" ]; then
			cat {input.contaminants_fasta} > {output.phix_contaminants_fasta}
			minimap2 -ax map-ont {output.phix_contaminants_fasta} {input.trimmed_data} | samtools fastq -c 6 -n -f 4 - > {output.fastq}
		else
			: > {output.phix_contaminants_fasta:q}
			cp -- {input.trimmed_data:q} {output.fastq:q}
		fi
		echo $(( $(zgrep -Ec "$" {output.fastq}) / 4 )) > {output.size}
		"""

if CONTAMINANTS==["GCF_000819615.1"]:
	ruleorder: get_size_nanopore > remove_contaminants_nanopore

	rule get_size_nanopore:
		benchmark:
			dirs_dict["BENCHMARKS"] + "/get_size_nanopore/sample_nanopore={sample_nanopore}.tsv"
		input:
			trimmed_data=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq.gz",
		output:
			fastq=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.fastq.gz"),
			size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean_read_count.tot.txt",
		message:
			"Calculating number of reads"
		conda:
			dirs_dict["ENVS_DIR"]+ "/env1.yaml"
		threads: 1
		shell:
			"""
			cp {input.trimmed_data} {output.fastq}
			echo $(( $(zgrep -Ec "$" {output.fastq}) / 4 )) > {output.size}
			"""


rule postQualityCheckNanopore:
	input:
		fastq=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.fastq.gz"),
	output:
		nanoqc_dir=temp(directory(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanoqc_post")),
		nanoqc=dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_postQC.html",
	message:
		"Performing nanoQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/postQualityCheckNanopore/sample_nanopore={sample_nanopore}.tsv"
#	threads: 1
	shell:
		"""
		nanoQC -o {output.nanoqc_dir} {input.fastq}
		mv {output.nanoqc_dir}/nanoQC.html {output.nanoqc}
		"""

rule qualityStatsNanopore_pre:
	input:
		raw_data=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore.fastq.gz",
	output:
		nanostats=dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanostats_preQC.html",
	message:
		"Performing nanoQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/qualityStatsNanopore_pre/sample_nanopore={sample_nanopore}.tsv"
#	threads: 1
	shell:
		"""
		NanoStat --fastq {input.raw_data} > {output.nanostats}
		"""

rule qualityStatsNanopore_post:
	input:
		fastq=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.fastq.gz"),
	output:
		nanostats=dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanostats_postQC.html",
	message:
		"Performing nanoQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/qualityStatsNanopore_post/sample_nanopore={sample_nanopore}.tsv"
#	threads: 1
	shell:
		"""
		NanoStat --fastq {input.fastq} > {output.nanostats}
		"""

rule subsampleReadsNanopore:
	input:
		nano_sizes=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean_read_count.tot.txt", sample_nanopore=NANOPORE_SAMPLES),
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.fastq",
	output:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.sub.fastq",
		size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean_read_count.sub.txt",
	message:
		"Subsampling Nanopore reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/subsampleReadsNanopore/sample_nanopore={sample_nanopore}.tsv"
	params:
		sizes=dirs_dict["CLEAN_DATA_DIR"] + "/*_nanopore_clean_read_count.tot.txt"
	threads: 1
	resources:
		mem_mb=4000
	shell:
		"""
		nanopore=$( cat {params.sizes} | sort -n | head -1 )
		reformat.sh in={input.nanopore} out={output.nanopore} reads=$nanopore ignorebadquality
		echo $(( $(grep -Ec "$" {output.nanopore}) / 4 )) > {output.size}
		"""
