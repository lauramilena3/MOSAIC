#ruleorder: trim_adapters_quality_illumina_PE > trim_adapters_quality_illumina_SE
#ruleorder: listContaminants_PE > listContaminants_SE
#ruleorder: removeContaminants_PE > removeContaminants_SE
#ruleorder: subsampleReadsIllumina_PE > subsampleReadsIllumina_SE
#ruleorder: normalizeReads_PE > normalizeReads_SE
#ruleorder: postQualityCheckIlluminaPE > postQualityCheckIlluminaSE

rule download_SRA:
	input:
		sratoolkit="tools/sratoolkit.2.10.0-ubuntu64"
	output:
		forward_file=(dirs_dict["RAW_DATA_DIR"] + "/{SRA}_pass_1.fastq"),
		reverse_file=(dirs_dict["RAW_DATA_DIR"] + "/{SRA}_pass_2.fastq"),
	params:
		SRA_dir=dirs_dict["RAW_DATA_DIR"],
	message:
		"Downloading SRA run"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
#	threads: 1
	shell:
		"""
		{input.sratoolkit}/bin/fastq-dump --outdir {params.SRA_dir} --skip-technical --readids --read-filter pass \\
		--dumpbase --split-files --clip -N 0 -M 0 {wildcards.SRA}
		"""

rule countReads_gz:
	input:
		fastq="{fastq_name}.fastq.gz",
	output:
		counts="{fastq_name}_read_count.txt",
	message:
		"Counting reads on fastq file"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	group:
		"read_counts_gz"
	resources:
		runtime_min= 5,
		mem_mb= 1000,
	shell:
		"""
		echo $(( $(zgrep -Ec "$" {input.fastq}) / 4 )) > {output.counts} 
		"""

rule countReads:
	input:
		fastq="{fastq_name}.fastq",
	output:
		counts="{fastq_name}_read_count.txt",
	message:
		"Counting reads on fastq file"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	group:
		"read_counts"
	resources:
		runtime_min= 5,
		mem_mb= 1000,
	shell:
		"""
		echo $(( $(grep -Ec "$" {input.fastq}) / 4 )) > {output.counts} 
		"""

rule fastQC_pre:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"] + "/{fastq_name}.fastq.gz"
	output:
		html=temp(dirs_dict["RAW_DATA_DIR"] + "/{fastq_name}_fastqc.html"),
		zipped=(dirs_dict["RAW_DATA_DIR"] + "/{fastq_name}_fastqc.zip")
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/qualityCheckIllumina/{fastq_name}_pre_qc.tsv"
	resources:
		runtime_min= 412,
		mem_mb= 500,
	shell:
		"""
		fastqc {input}
		"""

rule fastQC_post:
	input:
		raw_fastq=dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}.fastq.gz"
	output:
		html=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}_fastqc.html"),
		zipped=(dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}_fastqc.zip")
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/qualityCheckIllumina/{fastq_name}_post_qc.tsv"
	resources:
		runtime_min= 412,
		mem_mb= 500,
	shell:
		"""
		fastqc {input}
		"""


rule superDeduper_pcr:
	input:
		forward_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq.gz",
		reverse_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + ".fastq.gz",
	output:
		duplicate_stats=(dirs_dict["QC_DIR"] + "/{sample}_stats_pcr_duplicates.log"),
		deduplicate=temp(dirs_dict["QC_DIR"] + "/{sample}_stats_pcr_duplicates.out"),
	message:
		"Detect PCR duplicates"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/SuperDeduper/{sample}_pcr_duplicates.tsv"
	resources:
		runtime_min= 30,
		mem_mb= 1000,
	shell:
		"""
		hts_SuperDeduper -L {output.duplicate_stats} -1 {input.forward_file} -2 {input.reverse_file} > {output.deduplicate}
		"""

rule trim_adapters_quality_illumina_PE:
	input:
		forward_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq.gz",
		reverse_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + ".fastq.gz",
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		forward_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq.gz"),
		reverse_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq.gz"),
		merged_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.tot.fastq.gz"),
	params:
		adapters=dirs_dict["ADAPTERS_DIR"] + "/" + config['adapters_file']
	message:
		"Trimming Illumina Adapters with Trimmomatic"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/trim_adapters_quality_illumina_PE/{sample}.tsv"
	resources:
		runtime_min= 350,
		mem_mb= 1500,
	threads: 8
	shell:
		"""
		trimmomatic PE -threads {threads} -phred33 {input.forward_file} {input.reverse_file} \
			{output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
			ILLUMINACLIP:{params.adapters}:2:30:10:1:true LEADING:{config[trimmomatic_leading]} TRAILING:{config[trimmomatic_trailing]} \
			SLIDINGWINDOW:{config[trimmomatic_window_size]}:{config[trimmomatic_window_quality]} MINLEN:{config[trimmomatic_minlen]}
		cat {output.forward_unpaired} {output.reverse_unpaired} > {output.merged_unpaired}
		"""

rule sourmash_sketch_trim:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
	output:
		manysketch_csv=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_manysketch.csv"),
		sketch=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sourmash.sig.zip"),
	params:
		sample="{sample}"
	message:
		"Building sketches with sourmash"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/sourmash/{sample}_sketch.tsv"
	threads: 8
	shell:
		"""
		echo name,read1,read2 > {output.manysketch_csv}
		echo {params.sample},{input.forward_paired},{input.reverse_paired} >> {output.manysketch_csv}
		sourmash scripts manysketch {output.manysketch_csv} -p k=31,k=51,abund,scaled=1000,DNA -o {output.sketch} -c {threads}
		"""

rule sourmash_gather:
    input:
        sketch=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sourmash.sig.zip",
        sourmash_rocksdb=config["sourmash_rocksdb"],
    output:
        gather=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_gather_sourmash.csv",
    params:
        threshold_bp=50000,
    message:
        "Metagenome containment with sourmash fastmultigather"
    conda:
        dirs_dict["ENVS_DIR"] + "/sourmash.yaml"
    benchmark:
        dirs_dict["BENCHMARKS"] + "/sourmash/{sample}_gather.tsv"
    threads: 8
    shell:
        """
        sourmash scripts fastmultigather \
            {input.sketch} \
            {input.sourmash_rocksdb} \
            -c {threads} \
            -o {output.gather} \
            -t {params.threshold_bp} \
            -s 1000
        """

rule sourmash_tax:
	input:
		gather=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_gather_sourmash.csv"),
		sourmash_tax=config['sourmash_tax'],
	output:
		kreport=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sourmash.kreport.txt"),
	params:
		sample="{sample}_sourmash",
		outdir=(dirs_dict["CLEAN_DATA_DIR"]),
	message:
		"Assigning taxonomy with sourmash tax"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/sourmash/{sample}_tax.tsv"
	threads: 1
	shell:
		"""
		sourmash tax metagenome --gather-csv {input.gather} -t {input.sourmash_tax}  -o {params.sample} \
			--output-format kreport --rank species -f --output-dir {params.outdir}
		"""

rule contaminants_KRAKEN:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		merged_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.tot.fastq.gz"),
		kraken_db=(config['kraken_db']),
		kraken_tools=(config['kraken_tools']),
	output:
		kraken_output_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_tot.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_tot.csv"),
		kraken_domain=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_domains_tot.csv"),
		kraken_output_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_unpaired_tot.csv"),
		kraken_report_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_unpaired_tot.csv"),
	params:
		kraken_db=config['kraken_db'],
	message:
		"Assesing contamination with kraken2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kraken/{sample}_preliminary.tsv"
	threads: 16
	resources:
		runtime_min= 15,
		mem_mb= 18000,
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output_paired} --report {output.kraken_report_paired}
		grep -P 'D\t' {output.kraken_report_paired} | sort -r > {output.kraken_domain}
		#UNPAIRED
		kraken2 --db {params.kraken_db} --threads {threads} {input.merged_unpaired}  \
			--output {output.kraken_output_unpaired} --report {output.kraken_report_unpaired}
		"""
		# python {input.kraken_tools}/combine_kreports.py \
		# 	-r {output.kraken_report_paired} {output.kraken_report_unpaired} \
		# 	-o {output.kraken_report_combined}
		# ""

rule contaminants_KRAKEN_microbial:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		kraken_db=(config['kraken_db_nt']),
	output:
		kraken_output_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_microbial.tot.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_microbial.tot.csv"),
	params:
		kraken_db=config['kraken_db_nt'],
	message:
		"Assesing contamination with kraken2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kraken/{sample}_preliminary_microbial.tsv"
	threads: 32
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output_paired} --report {output.kraken_report_paired} \
			--report-minimizer-data
		"""

rule remove_euk:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		merged_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.tot.fastq.gz"),
		kraken_output_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_tot.csv"),
		kraken_output_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_unpaired_tot.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_tot.csv"),
		kraken_report_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_unpaired_tot.csv"),
		kraken_tools=(config['kraken_tools']),
	output:
		forward_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_noEuk.tot.fastq"),
		reverse_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_noEuk.tot.fastq"),
		unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_noEuk.tot.fastq"),
	message:
		"Removing eukaryotic reads with Kraken"
	params:
		# unclassified_name_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken_paired_R#.tot.fastq",
		host_taxid=config["contaminants_taxid"] 
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1_kraken.yaml"
	threads: 4
	benchmark:
		dirs_dict["BENCHMARKS"] +"/remove_euk_PE/{sample}.tsv"
	resources:
		mem_mb=40000,
		runtime_min= 1100,
	shell:
		"""
		python {input.kraken_tools}/extract_kraken_reads.py -k {input.kraken_output_paired} \
			-s1 {input.forward_paired} -s2 {input.reverse_paired} \
			-o {output.forward_paired} -o2 {output.reverse_paired} \
			--exclude -t {params.host_taxid} --include-children -r {input.kraken_report_paired} --fastq-output
		python {input.kraken_tools}/extract_kraken_reads.py -k {input.kraken_output_unpaired} \
			-s {input.merged_unpaired} -o {output.unpaired} --exclude -t {params.host_taxid} --include-children \
			-r {input.kraken_report_unpaired} --fastq-output
		"""

def remove_user_contaminants_forward(wildcards):
	if REMOVE_EUK:
		return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_forward_paired_noEuk.tot.fastq"
	return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_forward_paired.fastq.gz"

def remove_user_contaminants_reverse(wildcards):
	if REMOVE_EUK:
		return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_reverse_paired_noEuk.tot.fastq"
	return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_reverse_paired.fastq.gz"

def remove_user_contaminants_unpaired(wildcards):
	if REMOVE_EUK:
		return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_unpaired_noEuk.tot.fastq"
	return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_merged_unpaired.tot.fastq.gz"

rule remove_user_contaminants_PE:
	input:
		forward_paired=remove_user_contaminants_forward,
		reverse_paired=remove_user_contaminants_reverse,
		unpaired=remove_user_contaminants_unpaired,
		contaminants_fasta=expand(dirs_dict["CONTAMINANTS_DIR_DB"] +"/{contaminants}.fasta",contaminants=CONTAMINANTS),
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz"),
		phix_contaminants_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/{sample}_contaminants.fasta",
		stats=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_contaminant_stats_bbduk.tot.txt"),
	message:
		"Removing phiX174 and user given contaminants with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/remove_contaminants_PE/{sample}.tsv"
	threads: 4
	resources:
		mem_mb=40000,
		runtime_min= 15,
	shell:
		"""
		cat {input.contaminants_fasta} > {output.phix_contaminants_fasta}
		#PE
		#PAIRED
		bbduk.sh -Xmx{resources.mem_mb}m in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
			ref={output.phix_contaminants_fasta} k=31 hdist=1 threads={threads} stats={output.stats}
		#UNPAIRED
		bbduk.sh -Xmx{resources.mem_mb}m in={input.unpaired} out={output.unpaired} ref={output.phix_contaminants_fasta} k=31 hdist=1 threads={threads}
		"""

rule contaminants_KRAKEN_clean:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz",
		kraken_db=(config['kraken_db']),
		kraken_tools=(config['kraken_tools']),
	output:
		kraken_output_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_clean_tot.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_clean_tot.csv"),
	params:
		kraken_db=config['kraken_db'],
	message:
		"Assesing taxonomy with kraken2 on clean reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kraken/{sample}_clean.tsv"
	priority: 1
	threads: 8
	resources:
		runtime_min= 15,
		mem_mb= 18000,
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output_paired} --report {output.kraken_report_paired} \
			--report-minimizer-data
		"""

rule preMultiQC:
	input:
		#html=expand(dirs_dict["RAW_DATA_DIR"]+"/{sample}_{reads}_fastqc.html", sample=SAMPLES, reads=READ_TYPES),
		zipped=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{reads}_fastqc.zip", sample=SAMPLES, reads=READ_TYPES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html",
		multiqc_txt=dirs_dict["QC_DIR"]+ "/preQC_illumina_report_data/multiqc_fastqc.txt",
	params:
		fastqc_dir=dirs_dict["RAW_DATA_DIR"],
		html_name="preQC_illumina_report.html",
		multiqc_dir=dirs_dict["QC_DIR"],
	message:
		"Generating MultiQC report"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/multiQC/multiqc_pre.tsv"
	resources:
		runtime_min= 5,
		mem_mb= 4000,
	shell:
		"""
		multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule postMultiQC:
	input:
		# html_forward=expand(dirs_dict["CLEAN_DATA_DIR"]  + "/{sample}_forward_paired_clean.tot_fastqc.html", sample=SAMPLES),
		zipped_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot_fastqc.zip", sample=SAMPLES),
		# html_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot_fastqc.html", sample=SAMPLES),
		zipped_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot_fastqc.zip", sample=SAMPLES),
		# html_unpaired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot_fastqc.html", sample=SAMPLES),
		zipped_unpaired=expand(dirs_dict["CLEAN_DATA_DIR"]  + "/{sample}_unpaired_clean.tot_fastqc.zip", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/postQC_illumina_report.html",
		multiqc_txt=dirs_dict["QC_DIR"]+ "/postQC_illumina_report_data/multiqc_fastqc.txt",
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
		html_name="postQC_illumina_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	priority: 1
	benchmark:
		dirs_dict["BENCHMARKS"] +"/multiQC/multiqc_post.tsv"
	resources:
		runtime_min= 5,
		mem_mb= 4000,
	shell:
		"""
		multiqc -f {params.fastqc_dir}/*zip -o {params.multiqc_dir} -n {params.html_name}
		"""

rule prekrakenMultiQC:
	input:
		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_tot.csv", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/pre_decontamination_kraken_multiqc_report.html"
		# 		multiqc=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html",
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
		html_name="pre_decontamination_kraken_multiqc_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report kraken pre"
	priority: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/multiQC/multiqc_kraken_pre.tsv"
	resources:
		runtime_min= 5,
		mem_mb= 4000,
	shell:
		"""
		multiqc -f {input} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule postkrakenMultiQC:
	input:
		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_clean_tot.csv", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/post_decontamination_kraken_multiqc_report.html"
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
		html_name="post_decontamination_kraken_multiqc_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report kraken post"
	priority: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/multiQC/multiqc_kraken_post.tsv"
	resources:
		runtime_min= 5,
		mem_mb= 4000,
	shell:
		"""
		multiqc -f {input} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule krakenMicrobialMultiQC:
	input:
		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_microbial.tot.csv", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/microbial_kraken_multiqc_report.html"
		# 		multiqc=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html",
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
		html_name="microbial_kraken_multiqc_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report kraken pre"
	priority: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/multiQC/multiqc_kraken_microbial.tsv"
	shell:
		"""
		multiqc -f {input} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule normalizeReads_PE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq.gz",
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq.gz"),
		histogram_pre=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_count_histogram_pre.{sampling}.txt"),
		histogram_post=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_count_histogram_post.{sampling}.txt"),
		peaks=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_count_peaks.{sampling}.txt"),
	message:
		"Normalizing reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/normalizeReads_PE/{sample}_{sampling}.tsv"
	params:
		min_depth=config['min_norm'],
		max_depth=config['max_norm']
	threads: 16
	priority: 1
	wildcard_constraints:
		sampling="tot|sub"  
	resources:
		mem_mb=MEMORY_ECORR,
		runtime_min= 125,
	shell:
		"""
		#PE
		#paired
		bbnorm.sh -Xmx{resources.mem_mb}m in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
			target={params.max_depth} mindepth={params.min_depth} t={threads} khist={output.histogram_pre} peaks={output.peaks} khistout={output.histogram_post}
		#unpaired
		bbnorm.sh -Xmx{resources.mem_mb}m in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth} threads={threads}
		"""

rule concatenate_subassembly:
	input:
		forward_paired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",sample=SAMPLES),
		reverse_paired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",sample=SAMPLES),
		unpaired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz",sample=SAMPLES),
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/ALL_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/ALL_reverse_paired_clean.tot.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/ALL_unpaired_clean.tot.fastq.gz",
	message:
		"Concatenating clean reads for cross assembly"
	shell:
		"""
		cat {input.forward_paired} > {output.forward_paired}
		cat {input.reverse_paired} > {output.reverse_paired}
		cat {input.unpaired} > {output.unpaired}
		"""
rule kmer_rarefraction:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		histogram=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.{sampling}.csv"),
	message:
		"Counting unique reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kmer_rarefraction/{sample}_{sampling}.tsv"
	threads: 1
	resources:
		mem_mb=MEMORY_ECORR,
		runtime_min= 412,
	shell:
		"""
		bbcountunique.sh -Xmx{resources.mem_mb}m in1={input.forward_paired} in2={input.reverse_paired} out={output.histogram} interval={config[kmer_window]}
		"""
