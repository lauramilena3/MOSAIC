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
		trimmomatic_values=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_trimmomatic_values.txt"),
	params:
		adapters=dirs_dict["ADAPTERS_DIR"] + "/" + config['adapters_file']
	message:
		"Trimming Illumina Adapters with Trimmomatic"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/trim_adapters_quality_illumina_PE/{sample}.tsv"
	threads: 8
	shell:
		"""
		echo leading {config[trimmomatic_leading]} trailing {config[trimmomatic_trailing]} winsize {config[trimmomatic_window_size]} \
			winqual {config[trimmomatic_window_quality]} minlnth {config[trimmomatic_minlen]} > {output.trimmomatic_values}
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
		sourmash scripts manysketch {output.manysketch_csv} -p k=31,abund -o {output.sketch} -c {threads}
		"""


#
# rule trim_adapters_quality_illumina_SE:
# 	input:
# 		forward_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq",
# 		qc_report=dirs_dict["QC_DIR"]+ "/pre_processing_multiqc_report.html"
# 	output:
# 		forward_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq"),
# 	params:
# 		adapters=dirs_dict["ADAPTERS_DIR"] + "/" + config['adapters_file']
# 	message:
# 		"Trimming Illumina Adapters with Trimmomatic"
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/trim_adapters_quality_illumina_SE/{sample}.tsv"
# 	threads: 2
# 	shell:
# 		"""
# 		trimmomatic SE -threads {threads} -phred33 {input.forward_file} {output.forward_unpaired} \
# 		ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:{config[trimmomatic_leading]} TRAILING:{config[trimmomatic_trailing]} \
# 		SLIDINGWINDOW:{config[trimmomatic_window_size]}:{config[trimmomatic_window_quality]} MINLEN:{config[trimmomatic_minlen]}
# 		"""

# rule contaminants_KRAKEN_pre:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
# 		kraken_db=(config['kraken_db']),
# 		kraken_tools=(config['kraken_tools']),
# 	output:
# 		kraken_output_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_pre.csv"),
# 		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_pre.csv"),
# 		kraken_domain=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_domains_pre.csv"),
# 	params:
# 		kraken_db=config['kraken_db'],
# 	message:
# 		"Assesing contamination with kraken2"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/kmer_rarefraction/{sample}_tot.tsv"
# 	threads: 8
# 	shell:
# 		"""
# 		kraken2 --db {params.kraken_db} --threads {threads} \
# 			--paired {input.forward_paired} {input.reverse_paired} \
# 			--output {output.kraken_output_paired} --report {output.kraken_report_paired}
# 		grep -P 'D\t' {output.kraken_report_paired} | sort -r > {output.kraken_domain}
# 		"""

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
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kraken/{sample}_preliminary.tsv"
	threads: 16
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

# rule contaminants_KRAKEN_microbial:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
# 		kraken_db=(config['kraken_db_nt']),
# 	output:
# 		kraken_output_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_microbial.tot.csv"),
# 		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_microbial.tot.csv"),
# 	params:
# 		kraken_db=config['kraken_db_nt'],
# 	message:
# 		"Assesing contamination with kraken2"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/kraken/{sample}_preliminary_microbial.tsv"
# 	threads: 32
# 	shell:
# 		"""
# 		kraken2 --db {params.kraken_db} --threads {threads} \
# 			--paired {input.forward_paired} {input.reverse_paired} \
# 			--output {output.kraken_output_paired} --report {output.kraken_report_paired} \
# 			--report-minimizer-data
		# """

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
		host_taxid=2759
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 4
	benchmark:
		dirs_dict["BENCHMARKS"] +"/remove_euk_PE/{sample}.tsv"
	resources:
		mem_gb=40
	shell:
		"""
		python {input.kraken_tools}/extract_kraken_reads.py -k {input.kraken_output_paired} \
			-s1 {input.forward_paired} -s2 {input.reverse_paired} \
			-o {output.forward_paired} -o2 {output.reverse_paired} \
			--exclude --taxid {params.host_taxid} --include-children -r {input.kraken_report_paired} --fastq-output
		python {input.kraken_tools}/extract_kraken_reads.py -k {input.kraken_output_unpaired} \
			-s {input.merged_unpaired} -o {output.unpaired} --exclude --taxid {params.host_taxid} --include-children \
			-r {input.kraken_report_unpaired} --fastq-output
		"""

rule remove_user_contaminants_PE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_noEuk.tot.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_noEuk.tot.fastq"),
		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_noEuk.tot.fastq"),
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
		mem_gb=40
	shell:
		"""
		cat {input.contaminants_fasta} > {output.phix_contaminants_fasta}
		#PE
		#PAIRED
		bbduk.sh -Xmx{resources.mem_gb}g in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
			ref={output.phix_contaminants_fasta} k=31 hdist=1 threads={threads} stats={output.stats}
		#UNPAIRED
		bbduk.sh -Xmx{resources.mem_gb}g in={input.unpaired} out={output.unpaired} ref={output.phix_contaminants_fasta} k=31 hdist=1 threads={threads}
		"""

# rule remove_human_PE:
# 	input:1
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_bbduk.tot.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_bbduk.tot.fastq"),
# 		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired_bbduk.tot.fastq"),
# 		kraken_db_human=(config['kraken_db_human']),
# 	output:
# 		forward_paired_temp=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken_paired_R_1.tot.fastq"),
# 		reverse_paired_temp=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken_paired_R_2.tot.fastq"),
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq",
# 		paired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.tot.txt",
# 		unpaired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.txt",
# 		kraken_output_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-kraken2-out_paired.txt",
# 		kraken_report_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-kraken2-report_paired.txt",
# 		kraken_output_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-kraken2-out_unpaired.txt",
# 		kraken_report_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-kraken2-report_unpaired.txt",
# 	message:
# 		"Removing human reads with Kraken"
# 	params:
# 		unclassified_name_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken_paired_R#.tot.fastq",
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	threads: 4
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/remove_human_PE/{sample}.tsv"
# 	resources:
# 		mem_gb=40
# 	shell:
# 		"""
# 		#PAIRED
# 		kraken2 --paired --db {input.kraken_db_human} --threads {threads} --output {output.kraken_output_paired} \
# 				--report {output.kraken_report_paired} --unclassified-out {params.unclassified_name_paired} \
# 				{input.forward_paired} {input.reverse_paired}
# 		cp {output.forward_paired_temp} {output.forward_paired}
# 		cp {output.reverse_paired_temp} {output.reverse_paired}
# 		grep -c "^@" {output.forward_paired} > {output.paired_size}
# 		#UNPAIRED
# 		kraken2 --db {input.kraken_db_human} --threads {threads} --output {output.kraken_output_unpaired} \
# 				--report {output.kraken_report_unpaired} --unclassified-out {output.unpaired} {input.unpaired}
# 		grep -c "^@" {output.unpaired} > {output.unpaired_size} ||  echo "0" > {output.unpaired_size}
# 		"""

#
# rule postQualityCheckIlluminaPE:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq",
# 	output:
# 		html_forward=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_forward_paired_clean.tot_fastqc.html"),
# 		zipped_forward=(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_forward_paired_clean.tot_fastqc.zip"),
# 		html_reverse=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_clean.tot_fastqc.html"),
# 		zipped_reverse=(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_clean.tot_fastqc.zip"),
# 		html_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_unpaired_clean.tot_fastqc.html"),
# 		zipped_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_unpaired_clean.tot_fastqc.zip"),
# 	params:
# 		postQC_dir=dirs_dict["CLEAN_DATA_DIR"] +"/postQC",
# 	message:
# 		"Performing fastqQC statistics"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/QC.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/postQualityCheckIlluminaPE/{sample}.tsv"
# 	shell:
# 		"""
# 		fastqc {input.forward_paired} -o {params.postQC_dir}
# 		fastqc {input.reverse_paired} -o {params.postQC_dir}
# 		fastqc {input.unpaired} -o {params.postQC_dir}
# 		"""

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
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kraken/{sample}_clean.tsv"
	priority: 1
	threads: 8
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output_paired} --report {output.kraken_report_paired} \
			--report-minimizer-data
		"""

# rule read_classification_BRACKEN_pre:
# 	input:
# 		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_tot.csv"),
# 		kraken_db=(config['kraken_db']),
# 		bracken_checkpoint=config['kraken_db'] + "../bracken_db_ckeckpoint.txt",
# 	output:
# 		bracken_report_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_bracken_{level}_report_paired_pre_tot.csv",
# 	message:
# 		"Creating taxonomic reports with Bracken"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/bracken/{sample}_{level}_tot.tsv"
# 	shell:
# 		"""
# 		bracken -d {input.kraken_db}  -i {input.kraken_report_paired}  -o {output.bracken_report_paired} -l {wildcards.level} -t 4000 || true
# 		touch {output.bracken_report_paired}
# 		"""
# rule read_classification_BRACKEN:
# 	input:
# 		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_clean_tot.csv"),
# 		kraken_db=(config['kraken_db']),
# 		bracken_checkpoint=config['kraken_db'] + "../bracken_db_ckeckpoint.txt",
# 		read_count=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot_read_count.txt"),
# 	output:
# 		bracken_report_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_bracken_{level}_report_paired_tot.csv",
# 	message:
# 		"Creating taxonomic reports with Bracken"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/bracken/{sample}_{level}_tot.tsv"
# 	shell:
# 		"""
# 		read_count=$(cat {input.read_count})
# 		threshold=$((read_count*500 / 1000000))
# 		echo threshold $threshold
# 		bracken -d {input.kraken_db}  -i {input.kraken_report_paired}  -o {output.bracken_report_paired} -l {wildcards.level} -t $threshold || true
# 		touch {output.bracken_report_paired}
# 		"""

# rule read_classification_BRACKEN_microbial:
# 	input:
# 		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_microbial.tot.csv"),
# 		kraken_db=(config['kraken_db_nt']),
# 		read_count=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_read_count.txt"),
# 	output:
# 		bracken_report_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nt_bracken_{level}_report_paired_tot.csv",
# 	message:
# 		"Creating taxonomic reports with Bracken"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/bracken/{sample}_{level}_tot.tsv"
# 	shell:
# 		"""
# 		read_count=$(cat {input.read_count})
# 		threshold=$((read_count*500 / 1000000))
# 		echo threshold $threshold
# 		est_abundance.py -k {input.kraken_db}/database150mers.kmer_distrib -i {input.kraken_report_paired}  -o {output.bracken_report_paired} -l {wildcards.level} -t $threshold || true
# 		touch {output.bracken_report_paired}
# 		"""

# rule krakenUnique:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
# 		krakenUniq_db=(config['krakenUniq_db']),
# 	output:
# 		krakenuniq_output_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_krakenuniq_output_paired_clean_tot.csv"),
# 		krakenuniq_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_krakenuniq_report_paired_clean_tot.csv"),
# 	message:
# 		"Assesing taxonomy with kraken2 on clean reads"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env6.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/krakenUniq/{sample}_clean.tsv"
# 	priority: 1
# 	threads: 8
# 	shell:
# 		"""
# 		krakenuniq --report-file {output.krakenuniq_report_paired} --db {input.krakenUniq_db} --threads {threads} --output {output.krakenuniq_output_paired} --paired {input.forward_paired} {input.reverse_paired}
# 		"""

# rule read_classification_BRACKENUniq:
# 	input:
# 		bracken_dir=config['bracken_dir'],
# 		krakenuniq_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_krakenuniq_report_paired_clean_tot.csv"),
# 		krakenuniq_db=(config['krakenUniq_db']),
# 		brackenuniq_checkpoint=config['krakenUniq_db'] + "_brackenuniq_db_ckeckpoint.txt",
# 	output:
# 		brackenuniq_report_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_brackenuniq_{level}_report_paired_tot.csv",
# 	message:
# 		"Creating taxonomic reports with Bracken"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/bracken/{sample}_{level}_tot.tsv"
# 	shell:
# 		"""
# 		bracken -d {input.krakenuniq_db}  -i {input.krakenuniq_report_paired}  -o {output.brackenuniq_report_paired} -l {wildcards.level} -t 4000 || true
# 		touch {output.brackenuniq_report_paired}
# 		"""

#
# rule postQualityCheckIlluminaSE:
# 	input:
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/post_{sample}_unpaired_clean.tot.fastq",
# 	output:
# 		html=temp(dirs_dict["CLEAN_DATA_DIR"] + "/post_{sample}_unpaired_clean.tot_fastqc.html"),
# 		zipped=temp(dirs_dict["CLEAN_DATA_DIR"] + "/post_{sample}_unpaired_clean.tot_fastqc.zip")
# 	message:
# 		"Performing fastqQC statistics"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/QC.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/postQualityCheckIlluminaSE/{sample}.tsv"
# 	shell:
# 		"""
# 		fastqc {input}
# 		"""

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
	shell:
		"""
		multiqc -f {input} -o {params.multiqc_dir} -n {params.html_name}
		"""


# rule subsampleReadsIllumina_PE:
# 	input:
# 		unpaired_sizes=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot_read_count.txt", sample=SAMPLES),
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz"
# 	output:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.sub.fastq.gz"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.sub.fastq.gz"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.sub.fastq.gz",
# 	message:
# 		"Subsampling Illumina reads with BBtools"
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/subsampleReadsIllumina_PE/{sample}.tsv"
# 	params:
# 		max_subsample=int(int(config['max_subsample'])/2),
# 		files_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/*_unpaired_clean.tot.txt",
# 		files_paired=dirs_dict["CLEAN_DATA_DIR"] + "/*_paired_clean.tot.txt"
# 	threads: 1
# 	resources:
# 		mem_mb=4000
# 	shell:
# 		"""
# 		#paired
# 		paired=$( cat {params.files_paired} | sort -n | head -1 )
# 		p=$([ $paired -le {params.max_subsample} ] && echo "$paired" || echo {params.max_subsample})
# 		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads=$p
# 		#unpaired
# 		unpaired_temp=$( cat {params.files_unpaired} | sort -n | head -1 )
# 		un=$([ $unpaired_temp -le {params.max_subsample} ] && echo "$unpaired_temp" || echo {params.max_subsample})
# 		reads_left=$(({params.max_subsample} - ($paired*2)))
# 		unpaired=$([ $un -le $reads_left ] && echo "$un" || echo $reads_left )
# 		reformat.sh in={input.unpaired} out={output.unpaired} reads=$unpaired
# 		"""
#
# rule subsampleReadsIllumina_SE:
# 	input:
# 		unpaired_sizes=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.txt", sample=SAMPLES),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq"
# 	output:
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.sub.fastq",
# 		unpaired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.sub.txt"
# 	message:
# 		"Subsampling Illumina reads with BBtools"
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/subsampleReadsIllumina_SE/{sample}.tsv"
# 	params:
# 		max_subsample=config['max_subsample'],
# 		files_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/*_unpaired_clean.tot.txt"
# 	threads: 1
# 	resources:
# 		mem_mb=4000
# 	shell:
# 		"""
# 		#unpaired
# 		unpaired_temp=$( cat {params.files_unpaired} | sort -n | head -1 )
# 		echo $unpaired_temp
# 		un=$([ $unpaired_temp -le {params.max_subsample} ] && echo "$unpaired_temp" || echo {params.max_subsample})
# 		echo $un
# 		reformat.sh in={input.unpaired} out={output.unpaired} reads=$un
# 		grep -c "^@" {output.unpaired} > {output.unpaired_size}
# 		"""

# rule errorCorrectReads_PE:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
# 	output:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_ecc.{sampling}.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_ecc.{sampling}.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_ecc.{sampling}.fastq",
# 	message:
# 		"Error correcting reads with BBtools"
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/errorCorrectReads_PE/{sample}_{sampling}.tsv"
# 	threads: 8
# 	resources:
# 		mem_mb=MEMORY_ECORR
# 	shell:
# 		"""
# 		#PE
# 		#paired
# 		tadpole.sh -Xmx{resources.mem_mb}m in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} threads={threads} mode=correct
# 		#unpaired
# 		tadpole.sh -Xmx{resources.mem_mb}m in={input.unpaired} out={output.unpaired} threads={threads} mode=correct
# 		"""

rule normalizeReads_PE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq.gz",
	output:
		forward_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq.gz"),
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
	resources:
		mem_mb=MEMORY_ECORR
	shell:
		"""
		#PE
		#paired
		bbnorm.sh -Xmx{resources.mem_mb}m in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
			target={params.max_depth} mindepth={params.min_depth} t={threads} khist={output.histogram_pre} peaks={output.peaks} khistout={output.histogram_post}
		#unpaired
		bbnorm.sh -Xmx{resources.mem_mb}m in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth} threads={threads}
		"""

# rule normalizeReads_PE_noecc:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
# 	output:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm_noecc.{sampling}.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm_noecc.{sampling}.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm_noecc.{sampling}.fastq",
# 	message:
# 		"Normalizing reads with BBtools"
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/normalizeReads_PE/{sample}_{sampling}.tsv"
# 	params:
# 		min_depth=config['min_norm'],
# 		max_depth=config['max_norm']
# 	threads: 8
# 	resources:
# 		mem_mb=MEMORY_ECORR
# 	shell:
# 		"""
# 		#PE
# 		#paired
# 		bbnorm.sh -Xmx{resources.mem_mb}m in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
# 			target={params.max_depth} mindepth={params.min_depth} t={threads}
# 		#unpaired
# 		bbnorm.sh -Xmx{resources.mem_mb}m in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth} threads={threads}
# 		"""

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
# rule normalizeReads_SE:
# 	input:
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq"
# 	output:
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq"
# 	message:
# 		"Normalizing and error correcting reads with BBtools"
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/normalizeReads_SE/{sample}_{sampling}.tsv"
# 	params:
# 		min_depth=config['min_norm'],
# 		max_depth=config['max_norm']
# 	threads: 4
# 	resources:
# 		mem_mb=8000
# 	shell:
# 		"""
# 		#SE
# 		#unpaired
# 		bbnorm.sh -Xmx{resources.mem_mb}m ecc in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth} t={threads}
# 		"""

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
		mem_mb=MEMORY_ECORR
	shell:
		"""
		bbcountunique.sh -Xmx{resources.mem_mb}m in1={input.forward_paired} in2={input.reverse_paired} out={output.histogram} interval={config[kmer_window]}
		"""

