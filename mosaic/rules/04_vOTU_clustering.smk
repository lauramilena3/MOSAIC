# ruleorder: vOUTclustering_references>vOUTclustering
# ruleorder: filter_vOTUs_references>filter_vOTUs


def input_vOTU_clustering(wildcards):
	input_list=expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta",sample=SAMPLES)
	if NANOPORE:
		#input_list.append(dirs_dict["VIRAL_DIR"]+ "/{sample}_"+ LONG_ASSEMBLER + "_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta")
		#input_list=expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta",sample=SAMPLES)
		input_list.extend(expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_"+ LONG_ASSEMBLER + "_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta", sample=NANOPORE_SAMPLES))
		# print(input_list)
	return input_list



if len(config['additional_reference_contigs'])==0:

	rule vOUTclustering:
		input:
			positive_contigs=input_vOTU_clustering
		output:
			combined_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
			clusters=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + ".{sampling}_95-85.clstr",
			blastout=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + ".{sampling}-blastout.csv",
			aniout=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + ".{sampling}-aniout.csv",
			representative_list=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.txt",
			representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
			representative_lengths=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_lengths.{sampling}.txt",
		message:
			"Creating vOUTs with CheckV aniclust"
		conda:
			dirs_dict["ENVS_DIR"] + "/vir.yaml"
		benchmark:
			dirs_dict["BENCHMARKS"] +"/vOUTclustering/{sampling}.tsv"
		threads: 64
		shell:
			"""
			cat {input.positive_contigs} > {output.combined_positive_contigs}
			makeblastdb -in {output.combined_positive_contigs} -dbtype nucl -out {output.combined_positive_contigs}
			blastn -query {output.combined_positive_contigs} -db {output.combined_positive_contigs} -outfmt '6 std qlen slen' \
				-max_target_seqs 10000 -out {output.blastout} -num_threads {threads}
			python scripts/anicalc_checkv.py  -i {output.blastout} -o {output.aniout}
			python scripts/aniclust_checkv.py --fna {output.combined_positive_contigs} --ani {output.aniout} --out {output.clusters} --min_ani 95 --min_tcov 85 --min_qcov 0
			cut {output.clusters} -f1 > {output.representative_list}
			seqtk subseq {output.combined_positive_contigs} {output.representative_list} > {output.representatives}
			cat {output.representatives} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
				$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.representative_lengths}
			"""
	rule filter_vOTUs:
		input:
			representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
			high_qualty_list=dirs_dict["vOUT_DIR"] + "/checkV_high_quality.{sampling}.txt",
			representative_lengths=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_lengths.{sampling}.txt",
		output:
			filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
			filtered_list_temp=temp(dirs_dict["vOUT_DIR"]+ "/filtered_temp_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt"),
			filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		message:
			"Filtering vOTUs "
		conda:
			dirs_dict["ENVS_DIR"] + "/env1.yaml"
		benchmark:
			dirs_dict["BENCHMARKS"] +"/filter_vOTUs/{sampling}.tsv"
		threads: 4
		shell:
			"""
			cat {input.representative_lengths} | awk '$2>=10000' | cut -f1 >  {output.filtered_list_temp}
			cat {input.high_qualty_list} {output.filtered_list_temp} | sort | uniq > {output.filtered_list}
			seqtk subseq {input.representatives} {output.filtered_list} > {output.filtered_representatives}
			"""


else:

	rule vOUTclustering_references:
		input:
			positive_contigs=input_vOTU_clustering,
			additional_reference_contigs=config['additional_reference_contigs'],
		output:
			combined_positive_contigs=(dirs_dict["vOUT_DIR"]+ "/viral_and_reference_contigs.{sampling}.fasta"),
			clusters=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs.{sampling}_95-85.clstr",
			blastout=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs.{sampling}-blastout.csv",
			aniout=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs.{sampling}-aniout.csv",
			representative_list=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.txt",
			representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
			representative_lengths=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs_lengths.{sampling}.txt",
		message:
			"Creating vOUTs with Chechv aniclust"
		params:
			additional_reference_contigs=config['additional_reference_contigs']
		conda:
			dirs_dict["ENVS_DIR"] + "/vir.yaml"
		benchmark:
			dirs_dict["BENCHMARKS"] +"/vOUTclustering_references/{sampling}.tsv"
		threads: 64
		shell:
			"""
			cat {input.positive_contigs} {input.additional_reference_contigs} > {output.combined_positive_contigs}
			makeblastdb -in {output.combined_positive_contigs} -dbtype nucl -out {output.combined_positive_contigs}
			blastn -query {output.combined_positive_contigs} -db {output.combined_positive_contigs} -outfmt '6 std qlen slen' \
				-max_target_seqs 10000 -out {output.blastout} -num_threads {threads}
			python scripts/anicalc_checkv.py -i {output.blastout} -o {output.aniout}
			python scripts/aniclust_checkv.py --fna {output.combined_positive_contigs} --ani {output.aniout} --out {output.clusters} --min_ani 95 --min_tcov 85 --min_qcov 0
			cut {output.clusters} -f1 > {output.representative_list}
			seqtk subseq {output.combined_positive_contigs} {output.representative_list} > {output.representatives}
			cat {output.representatives} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
				$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.representative_lengths}
			"""

	rule filter_vOTUs_references:
		input:
			representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
			high_qualty_list=dirs_dict["vOUT_DIR"] + "/checkV_high_quality.{sampling}.txt",
			representative_lengths=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs_lengths.tot.txt",
		output:
			filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
			filtered_list_temp=temp(dirs_dict["vOUT_DIR"]+ "/filtered_temp_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt"),
			filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		message:
			"Filtering vOTUs "
		conda:
			dirs_dict["ENVS_DIR"] + "/env1.yaml"
		benchmark:
			dirs_dict["BENCHMARKS"] +"/filter_vOTUs/{sampling}.tsv"
		threads: 1
		shell:
			"""
			cat {input.representative_lengths} | awk '$2>=10000' | cut -f1 >>  {output.filtered_list_temp}
			cat {input.high_qualty_list} {output.filtered_list_temp}| sort | uniq > {output.filtered_list}
			seqtk subseq {input.representatives} {output.filtered_list} > {output.filtered_representatives}
			"""

rule getHighquality:
	input:
		quality_summary=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_checkV/quality_summary.tsv",
	output:
		high_qualty_list=dirs_dict["vOUT_DIR"] + "/checkV_high_quality.{sampling}.txt",
	message:
		"Filtering vOTUs "
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/filter_vOTUs/{sampling}.tsv"
	threads: 1
	shell:
		"""
		cat {input.quality_summary} | grep "High-quality"  | cut -f1 > {output.high_qualty_list}
		"""
