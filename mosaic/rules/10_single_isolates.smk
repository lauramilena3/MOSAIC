rule estimateGenomeCompletnessIsolates:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		checkv_db=(config['checkv_db']),
	output:
		quality_summary=dirs_dict["vOUT_DIR"] + "/checkV_isolates_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["vOUT_DIR"] + "/checkV_isolates_{sampling}/completeness.tsv",
		contamination=dirs_dict["vOUT_DIR"] + "/checkV_isolates_{sampling}/contamination.tsv",
	params:
		checkv_outdir=dirs_dict["vOUT_DIR"] + "/checkV_isolates_{sampling}",
		checkv_db=dirs_dict["vOUT_DIR"] + "/checkV_isolates_{sampling}",
	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/estimateGenomeCompletness/isolates_{sampling}.tsv"
	threads: 4
	shell:
		"""
		checkv contamination {input.positive_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv completeness {input.positive_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv complete_genomes {input.positive_contigs} {params.checkv_outdir}
		checkv quality_summary {input.positive_contigs} {params.checkv_outdir}
		"""

rule filter_isolates:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		quality_summary=dirs_dict["vOUT_DIR"] + "/checkV_isolates_{sampling}/quality_summary.tsv",
	output:
		filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_list.{sampling}.txt",
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates.{sampling}.fasta",
	message:
		"Filtering vOTUs "
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/filter_vOTUs/{sampling}.tsv"
	threads: 4
	shell:
		"""
		grep "High-quality" {input.quality_summary} | cut -f1 > {output.filtered_list}
		seqtk subseq {input.positive_contigs} {output.filtered_list} > {output.filtered_positive_contigs}
		"""

rule makeblastdb_isolates:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta",
	output:
		orfs_blast_db=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta.pdb",
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_makeblast_{sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 8
	shell:
		"""
		makeblastdb -in {input.filtered_positive_contigs} -dbtype prot
		"""

rule blastall_isolates:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta",
		blast=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta.pdb",
	output:
		blast_output=(dirs_dict["ANNOTATION"] + "/isolates_blastall.{sampling}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_blastall_{sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 8
	shell:
		"""
		blastp -num_threads {threads} -db {input.filtered_positive_contigs} -query {input.filtered_positive_contigs} \
			-outfmt "6 qseqid sseqid qstart qend qlen slen qcovs qcovhsp length pident evalue positive gaps" > {output.blast_output}
		"""

def inputDatabase(wildcards):
	if wildcards.db=="RefSeqViral":
		return config["RefSeqViral_protein_db"]
	elif wildcards.db=="GenBankViral":
		return config["GenBankViral_protein_db"]
	elif wildcards.db=="IMGVR":
		return config["IMGVR_protein_db"]

rule blastp_Reference_db:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta",
		blast_database=inputDatabase
	output:
		blast_output=(dirs_dict["ANNOTATION"] + "/isolates_blast_{db}.{sampling}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_{db}_{sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 32
	shell:
		"""
		blastp -num_threads {threads} -db {input.blast_database} -query {input.filtered_positive_contigs} \
			-outfmt "6 qseqid sseqid qstart qend qlen slen qcovs qcovhsp length pident evalue positive gaps" > {output.blast_output}
		"""

# rule blastp_IMGVR:
# 	input:
# 		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta",
# 		blast="/home/lmf/db/IMG_VR_09_2020/ORFs/IMGVR_ORF.faa",
# 	output:
# 		blast_output=(dirs_dict["ANNOTATION"] + "/isolates_blast_IMGVR.{sampling}.csv"),
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/viga.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_IMGVR_{sampling}.tsv"
# 	message:
# 		"Annotating contigs with BLAST"
# 	threads: 8
# 	shell:
# 		"""
# 		blastp -num_threads {threads} -db {input.blast} -query {input.filtered_positive_contigs} \
# 			-outfmt "6 qseqid sseqid qstart qend qlen slen qcovs qcovhsp length pident evalue positive gaps" > {output.blast_output}
# 		"""

rule blastp_database_lengths:
	input:
		blast="{fasta}.faa",
	output:
		length_temp=temp("{fasta}_length_temp.txt"),
		names_temp=temp("{fasta}_names_temp.txt"),
		length="{fasta}_lengths_and_names.txt",
		cummulative_length="{fasta}_coding_lengths.txt",
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/lengths_{fasta}.tsv"
	message:
		"Getting database ORF lengths"
	threads: 1
	shell:
		"""
		cat {input.blast} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
			$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.length_temp}
		cat {output.length_temp} | cut -f1 | cut -f1 -d' ' | rev | cut -d_ -f2- | rev > {output.names_temp}
		paste {output.length_temp} {output.names_temp} > {output.length}
		cat {output.length} | cut -f3,2 | awk '{{arr[$2]+=$1}} END {{for (i in arr) {{print i,arr[i]}}}}' > {output.cummulative_length}
		"""

def inputDatabaseLen(wildcards):
	if wildcards.db=="RefSeqViral":
		return (config["RefSeqViral_protein_db"].split(".faa")[0] + "_lengths_and_names.txt")
	elif wildcards.db=="GenBankViral":
		return (config["GenBankViral_protein_db"].split(".faa")[0] + "_lengths_and_names.txt")
	elif wildcards.db=="IMGVR":
		return (config["IMGVR_protein_db"].split(".faa")[0] + "_lengths_and_names.txt")

def inputDatabaseCoding(wildcards):
	if wildcards.db=="RefSeqViral":
		return (config["RefSeqViral_protein_db"].split(".faa")[0] + "_coding_lengths.txt")
	elif wildcards.db=="GenBankViral":
		return (config["GenBankViral_protein_db"].split(".faa")[0] + "_coding_lengths.txt")
	elif wildcards.db=="IMGVR":
		return (config["IMGVR_protein_db"].split(".faa")[0] + "_coding_lengths.txt")

rule get_relatives_list:
	input:
		blast_output=dirs_dict["ANNOTATION"] + "/isolates_blast_{db}.{sampling}.csv",
		isolate_len_file=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs_length.{sampling}.txt",
		isolate_coding=dirs_dict["ANNOTATION"]+ "/filtered_isolates_ORFs_coding_lengths.{sampling}.txt",
		reference_len_file=inputDatabaseLen,
		reference_coding=inputDatabaseCoding,
	output:
		blast_relatives_phages=(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_phages_{db}.{sampling}.txt"),
		blast_relatives_proteins=(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_ORFs_{db}.{sampling}.txt"),
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_{db}_{sampling}.tsv",
	message:
		"Extracting relative sequences NAMES"
	params:
		similarity_cutoff=20,
		length_cutoff=20,
	threads: 1
	run:
		# ALL TOGETHER
		import pandas as pd
		from Bio import SeqIO
		import seaborn as sns
		import matplotlib.pyplot as plt
		import time
		#------------------------
		# DEFINE CONSTANTS
		#------------------------
		#input
		blast_out_file=input.blast_output
		isolate_len_file=input.isolate_len_file
		isolates_coding_len=input.isolate_coding
		reference_len_file=input.reference_len_file
		reference_coding_len=input.reference_coding
		#output
		positive_phages=output.blast_relatives_phages
		positive_orfs=output.blast_relatives_proteins

		similarity_cutoff=params.similarity_cutoff
		len_cutoff=params.length_cutoff
		start = time.time()

		#------------------------
		# READ CODING LENGTH
		#------------------------
		print("Calculate coding length")
		print(time.time() - start)

		isolates_len_df=pd.read_csv(isolate_len_file, sep="\t", names=["orf", "aa_length", "phage"]).dropna()
		isolates_coding_len_df=pd.read_csv(isolates_coding_len, sep=" ", names=["phage", "aa_length"]).dropna()

		Reference_len_df=pd.read_csv(reference_len_file, sep="\t", names=["orf", "aa_length", "phage"]).dropna()
		Reference_coding_len_df=pd.read_csv(reference_coding_len, sep=" ", names=["phage", "aa_length"]).dropna()

		isolates_coding_len_df

		#-------------------------
		# PARSE BLAST RESULTS
		#-------------------------
		print("PARSE BLAST RESULTS")
		print(time.time() - start)

		#BLASTALL
		blast_results_reference=pd.read_csv(blast_out_file, names=["qseqid", "sseqid", "qstart", "qend", "qlen", "slen", "qcovs", "qcovhsp", "length", "pident", "evalue", "positive", "gaps"], sep="\t")
		blast_results_reference["phage_A"]=blast_results_reference["qseqid"].str.rsplit("_",1).str[0]
		blast_results_reference["phage_B"]=blast_results_reference["sseqid"].str.rsplit("_",1).str[0]
		blast_results_reference["match_positions"]=blast_results_reference["pident"]*blast_results_reference["length"]/100

		#-------------------------
		# CALCULATE BLAST IDENTITY
		#-------------------------
		print("CALCULATE BLAST IDENTITY")
		print(time.time() - start)

		grouped = blast_results_reference.groupby(by=["phage_A", "phage_B"])
		phage_A=[]
		phage_B=[]
		positive=[]
		for name, group in grouped:
		    best_hit=group.groupby(by=["qseqid"]).first()
		    phage_A.append(name[0])
		    phage_B.append(name[1])
		    positive.append(best_hit["match_positions"].sum())

		similarity_df=pd.DataFrame()
		similarity_df["Phage_A"]=phage_A
		similarity_df["Phage_B"]=phage_B
		similarity_df["Positive"]=positive


		similarity_df1=similarity_df.merge(isolates_coding_len_df, left_on="Phage_A", right_on="phage").drop(columns=["phage"])
		similarity_df1 = similarity_df1.rename(columns={'aa_length': 'lenA'})
		similarity_df2=similarity_df1.merge(Reference_coding_len_df, left_on="Phage_B", right_on="phage").drop(columns=["phage"])
		similarity_df2 = similarity_df2.rename(columns={'aa_length': 'lenB'})
		similarity_df2["simAB"]=similarity_df2["Positive"]*100/similarity_df2["lenA"]
		similarity_df2["ratioAB"]=similarity_df2["lenA"]*100/similarity_df2["lenB"]

		similarity_df2
		#-------------------------
		#FILTER POSITIVE HITS
		#-------------------------
		print("FILTER POSITIVE HITS")
		print(time.time() - start)

		similarity_df2=similarity_df2[similarity_df2["simAB"]>similarity_cutoff]
		similarity_df2=similarity_df2[(similarity_df2["ratioAB"]>(100-len_cutoff)) & (similarity_df2["ratioAB"]<(100*100/(100-len_cutoff)))]

		positive_hits_nucleotide=list(set(similarity_df2["Phage_B"].to_list()))
		positive_hits_nucleotide
		positive_hits_protein=list(set(Reference_len_df[Reference_len_df["phage"].isin(positive_hits_nucleotide)]["orf"].to_list()))
		positive_hits_protein

		#-------------------------
		#PRINT TO FILE
		#-------------------------
		print("PRINTING RESULTS TO FILE")
		print(time.time() - start)

		with open(positive_phages, 'w') as f:
		    for item in positive_hits_nucleotide:
		        f.write("%s\n" % item)

		with open(positive_orfs, 'w') as f:
		    for item in positive_hits_protein:
		        f.write("%s\n" % item)

def inputDatabaseExtract(wildcards):
	if wildcards.type=="ORFs":
		if wildcards.db=="RefSeqViral":
			return config["RefSeqViral_protein_db"]
		elif wildcards.db=="GenBankViral":
			return config["GenBankViral_protein_db"]
		elif wildcards.db=="IMGVR":
			return config["IMGVR_protein_db"]
	elif wildcards.type=="phages":
		if wildcards.db=="RefSeqViral":
			return config["RefSeqViral_db"]
		elif wildcards.db=="GenBankViral":
			return config["GenBankViral_db"]
		elif wildcards.db=="IMGVR":
			return config["IMGVR_db"]

rule get_relatives_fasta:
	input:
		blast_relatives=(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_{type}_{db}.{sampling}.txt"),
		blast_database=inputDatabaseExtract
	output:
		blast_relatives=(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_{type}_{db}.{sampling}.fasta"),
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_blast_relatives_extract_{type}_{db}.{sampling}.tsv"
	message:
		"Extracting relative sequences FASTA"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		seqtk subseq {input.blast_database} {input.blast_relatives} > {output.blast_relatives}
		"""

rule makeblastdb_relatives_ORFs:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta",
		blast_relatives=expand(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_ORFs_{db}.{{sampling}}.fasta", db=REFERENCE_DATABASES),
	output:
		cat_isolates_relatives=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs.{sampling}.fasta",
		orfs_blast_db=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs.{sampling}.fasta.pdb",
		length_temp=temp(dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs_length_temp.{sampling}.txt"),
		names_temp=temp(dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs_names_temp.{sampling}.txt"),
		length=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs_lengths_and_names.{sampling}.txt",
		cummulative_length=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs_coding_lengths.{sampling}.txt",
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_relatives_makeblast_{sampling}.tsv"
	message:
		"Creating ORFs relatives blast database"
	threads: 8
	shell:
		"""
		cat {input.filtered_positive_contigs} {input.blast_relatives} > {output.cat_isolates_relatives}
		cat {output.cat_isolates_relatives} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
			$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.length_temp}
		cat {output.length_temp} | cut -f1 | cut -f1 -d' ' | rev | cut -d_ -f2- | rev > {output.names_temp}
		paste {output.length_temp} {output.names_temp} > {output.length}
		cat {output.length} | cut -f3,2 | awk '{{arr[$2]+=$1}} END {{for (i in arr) {{print i,arr[i]}}}}' > {output.cummulative_length}
		makeblastdb -in {output.cat_isolates_relatives} -dbtype prot
		"""

rule blastall_relatives_ORFs:
	input:
		cat_isolates_relatives=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs.{sampling}.fasta",
		orfs_blast_db=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs.{sampling}.fasta.pdb",
	output:
		isolates_relatives_blastall=(dirs_dict["ANNOTATION"] + "/isolates_relatives_blastall_ORFs.{sampling}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_blastall_{sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 32
	shell:
		"""
		blastp -num_threads {threads} -db {input.cat_isolates_relatives} -query {input.cat_isolates_relatives} \
			-outfmt "6 qseqid sseqid qstart qend qlen slen qcovs qcovhsp length pident evalue positive gaps" > {output.isolates_relatives_blastall}
		"""

rule cat_relatives_phages:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates.{sampling}.fasta",
		blast_relatives=expand(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_phages_{db}.{{sampling}}.fasta", db=REFERENCE_DATABASES),
	output:
		cat_isolates_relatives=dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages.{sampling}.fasta",
		length_temp=temp(dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages_length_temp.{sampling}.txt"),
		names_temp=temp(dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages_names_temp.{sampling}.txt"),
		length=dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages_lengths_and_names.{sampling}.txt",
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_relatives_cat_{sampling}.tsv"
	message:
		"Creating ORFs relatives blast database"
	threads: 8
	shell:
		"""
		cat {input.filtered_positive_contigs} {input.blast_relatives} > {output.cat_isolates_relatives}
		cat {output.cat_isolates_relatives} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
			$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.length_temp}
		cat {output.length_temp} | cut -f1 | cut -f1 -d' ' | rev | cut -d_ -f2- | rev > {output.names_temp}
		paste {output.length_temp} {output.names_temp} > {output.length}
		"""

# rule blastall_relatives_phages:
# 	input:
# 		cat_isolates_relatives=dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages.{sampling}.fasta",
# 		orfs_blast_db=dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages.{sampling}.fasta.ndb",
# 	output:
# 		isolates_relatives_blastall=(dirs_dict["ANNOTATION"] + "/isolates_relatives_blastall_phages.{sampling}.csv"),
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/viga.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_blastall_{sampling}.tsv"
# 	message:
# 		"Annotating contigs with BLAST"
# 	threads: 8
# 	shell:
# 		"""
# 		blastn -num_threads {threads} -db {input.cat_isolates_relatives} -query {input.cat_isolates_relatives} \
# 			-outfmt "6 qseqid sseqid qstart qend qlen slen qcovs qcovhsp length pident evalue positive gaps" > {output.isolates_relatives_blastall}
# 		"""

rule viridic_relatives_phages:
	input:
		cat_isolates_relatives=dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages.{sampling}.fasta",
		viridic_singularity_folder=config["viridic_folder"]
	output:
		viridic_out=directory(dirs_dict["ANNOTATION"] + "/VIRIDIC_isolates_relatives_phages.{sampling}/"),
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/isolates_VIRIDIC_{sampling}.tsv"
	message:
		"Finding simmilartiy contigs with VIRIDIC"
	threads: 8
	shell:
		"""
		{input.viridic_singularity_folder}/viridic.bash projdir={output.viridic_out} in={input.cat_isolates_relatives}
		"""

rule genomad_host:
	input:
		host_fasta = dirs_dict["HOST_DIR"] + "/{host}.fasta",
		genomad_db= config['genomad_db'],
	output:
		genomad_outdir=directory(dirs_dict["HOST_DIR"] + "/{host}_geNomad"),
		positive_contigs=dirs_dict["HOST_DIR"] + "/prophages/{host}_prophages.fasta",
	params:
		viral_fasta=dirs_dict["HOST_DIR"] + "/{host}_geNomad/{host}_find_proviruses/{host}_provirus.fna",
	message:
		"Identifying prophages in host genome {wildcards.host} with geNomad"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/geNomad_viralID/{host}_host.tsv"
	threads: 16
	shell:
		"""
		genomad end-to-end --cleanup -t {threads} {input.host_fasta} {output.genomad_outdir} {input.genomad_db} 
		cp {params.viral_fasta} {output.positive_contigs}
		"""

rule buildBowtieDB_host:
	input:
		host_fasta = dirs_dict["HOST_DIR"] + "/{host}.fasta",
	output:
		contigs_bt2_1=temp(dirs_dict["HOST_DIR"]+ "/{host}.1.bt2"),
		contigs_bt2_2=temp(dirs_dict["HOST_DIR"]+ "/{host}.2.bt2"),
		contigs_bt2_3=temp(dirs_dict["HOST_DIR"]+ "/{host}.3.bt2"),
		contigs_bt2_4=temp(dirs_dict["HOST_DIR"]+ "/{host}.4.bt2"),
	# wildcard_constraints:
		# host = r"[^/]" 
	params:
		prefix=dirs_dict["HOST_DIR"]+ "/{host}",
	message:
		"Creating contig DB with Bowtie2"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{host}_bowtie_host.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		bowtie2-build {input.host_fasta} {params.prefix} --threads {threads}
		"""

rule map_to_host:
	input:
		contigs_bt2=dirs_dict["HOST_DIR"]+ "/{host}.1.bt2",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_flagstats_{sample}_vs_{host}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_flagstats_filtered_{sample}_vs_{host}.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}_covstats.txt",
	params:
		prefix=dirs_dict["HOST_DIR"]+ "/{host}",
	message:
		"Mapping reads to contigs"
	wildcard_constraints:
		host = r"[^/]+$" 
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadsToContigsPE/{sample}_vs_{host}_host.tsv"
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
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		"""