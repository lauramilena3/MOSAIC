# ruleorder: vOUTclustering_references>vOUTclustering
# ruleorder: filter_vOTUs_references>filter_vOTUs

def input_vOTU_clustering(wildcards):
	input_list=[]
	input_list.extend(expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta",sample=SAMPLES))
	if NANOPORE:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_"+ LONG_ASSEMBLER + "_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta", sample=NANOPORE_SAMPLES))
	if CROSS_ASSEMBLY:
		input_list.append(dirs_dict["VIRAL_DIR"]+ "/ALL_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta")
	if SUBASSEMBLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_positive_" + VIRAL_ID_TOOL + ".{{sampling}}.fasta", sample=SAMPLES, subsample=subsample_test))
	if len(config['additional_reference_contigs'])>0:
		input_list.extend(config['additional_reference_contigs'])
	return input_list

# if len(config['additional_reference_contigs'])==0:

rule derreplicate_assembly:
	input:
		positive_contigs=input_vOTU_clustering
	output:
		combined_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		derreplicated_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.fasta",
		derreplicated_tmp=directory(dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_tmp"),
	params:
		rep_name="combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated",
		rep_name_full=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_rep_seq.fasta",
		rep_temp="combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_tmp",
		dir_votu=dirs_dict["vOUT_DIR"],
	message:
		"Derreplicating assembled contigs with mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/contig_derreplication/{sampling}.tsv"
	threads: 16
	shell:
		"""
		cat {input.positive_contigs} > {output.combined_positive_contigs}
		cd {params.dir_votu}
		mmseqs easy-cluster --createdb-mode 1 --min-seq-id 1 -c 1 --cov-mode 1 {output.combined_positive_contigs} {params.rep_name} {params.rep_temp}
		mv {params.rep_name_full} {output.derreplicated_positive_contigs}
		"""

rule vOUTclustering:
	input:
		derreplicated_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.fasta",
	output:
		clusters=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + ".{sampling}_95-85.clstr",
		blastout=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + ".{sampling}-blastout.csv",
		aniout=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + ".{sampling}-aniout.csv",
	message:
		"Creating vOUTs with CheckV aniclust"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/vOUTclustering/{sampling}.tsv"
	threads: 64
	shell:
		"""
		makeblastdb -in {input.derreplicated_positive_contigs} -dbtype nucl -out {input.derreplicated_positive_contigs}
		blastn -query {input.derreplicated_positive_contigs} -db {input.derreplicated_positive_contigs} -outfmt '6 std qlen slen' \
			-max_target_seqs 10000 -out {output.blastout} -num_threads {threads}
		python scripts/anicalc_checkv.py  -i {output.blastout} -o {output.aniout}
		python scripts/aniclust_checkv.py --fna {input.derreplicated_positive_contigs} --ani {output.aniout} --out {output.clusters} --min_ani 95 --min_tcov 85 --min_qcov 0
		"""

def input_getHighQuality(wildcards):
	input_list=[]
	input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{{sampling}}/quality_summary.tsv",sample=SAMPLES)),
	if NANOPORE:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/nanopore_{sample}_" + LONG_ASSEMBLER + "_checkV_{{sampling}}/quality_summary.tsv", sample=NANOPORE_SAMPLES)),
	if CROSS_ASSEMBLY:
		input_list.append(dirs_dict["vOUT_DIR"] + "/ALL_checkV_{{sampling}}/quality_summary.tsv"),
	if SUBASSEMBLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_" + VIRAL_ID_TOOL + "_checkV_{{sampling}}/quality_summary.tsv", sample=SAMPLES, subsample=subsample_test)),
	if len(config['additional_reference_contigs'])>0:
		input_list.append(dirs_dict["ANNOTATION"] + "/user_reference_contigs_checkV/quality_summary.tsv"),
	return input_list

rule getHighQuality:
	input:
		input_getHighQuality,
	output:
		quality_summary_concat=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
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
		awk 'FNR>1' {input} > {output.quality_summary_concat}
		grep "High-quality" {output.quality_summary_concat} | cut -f1 > {output.high_qualty_list}
		"""

rule select_vOTU_representative:
	input:
		merged_summary=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		cluster_file=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + ".{sampling}_95-85.clstr",
	output:
		representatives=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list.{sampling}.csv",
		checkv_categories=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list_checkv_per_category.{sampling}.csv",
	params:
		samples=SAMPLES,
		contig_dir=dirs_dict["ASSEMBLY_DIR"],
		viral_dir=dirs_dict['VIRAL_DIR'],
		subassembly=SUBASSEMBLY,
		cross_assembly=CROSS_ASSEMBLY,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/05_vOTU_representative.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/05_vOTU_representative.py.ipynb"

rule vOUTclustering_get_new_references:
	input:
		derreplicated_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.fasta",
		representative_list=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list.{sampling}.csv",
	output:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		representative_lengths=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_lengths.{sampling}.txt",
	message:
		"Selecting new representatives with seqtk"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/vOUTclustering/{sampling}.tsv"
	threads: 64
	shell:
		"""
		seqtk subseq {input.derreplicated_positive_contigs} {input.representative_list} > {output.representatives}
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
	params:
		min_votu_len=config['min_votu_length']
	message:
		"Filtering vOTUs "
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/filter_vOTUs/{sampling}.tsv"
	threads: 4
	shell:
		"""
		cat {input.representative_lengths} | awk '$2>={params.min_votu_len}' | cut -f1 >  {output.filtered_list_temp}
		cat {input.high_qualty_list} {output.filtered_list_temp} | sort | uniq > {output.filtered_list}
		seqtk subseq {input.representatives} {output.filtered_list} > {output.filtered_representatives}
		"""

# rule filter_vOTUs:
# 	input:
# 		merged_summary=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
# 		cluster_file=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + ".{sampling}_95-85.clstr",
# 		vibrant_circular=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_circular.{sampling}.csv",
# 		vibrant_positive=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_list.{sampling}.csv",
# 		vibrant_quality=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_quality.{sampling}.csv",
# 		vibrant_summary=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_summary_results.{sampling}.csv",
# 		table_virsorter=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/final-viral-score.tsv",
# 		positive_list=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/positive_VS_list_{sampling}.txt",	
# 	output:
# 		representatives=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list.{sampling}.csv",
# 		checkv_categories=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list_checkv_per_category.{sampling}.csv",
# 	params:
# 		samples=SAMPLES,
# 		contig_dir=dirs_dict["ASSEMBLY_DIR"],
# 		viral_dir=dirs_dict['VIRAL_DIR'],
# 		subassembly=SUBASSEMBLY,
# 		cross_assembly=CROSS_ASSEMBLY,
# 	log:
# 		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/05_vOTU_representative.{sampling}.ipynb"
# 	notebook:
# 		dirs_dict["RAW_NOTEBOOKS"] + "/05_vOTU_representative.py.ipynb"


# else:

# 	rule vOUTclustering_references:
# 		input:
# 			positive_contigs=input_vOTU_clustering,
# 			additional_reference_contigs=config['additional_reference_contigs'],
# 		output:
# 			combined_positive_contigs=(dirs_dict["vOUT_DIR"]+ "/viral_and_reference_contigs.{sampling}.fasta"),
# 			clusters=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs.{sampling}_95-85.clstr",
# 			blastout=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs.{sampling}-blastout.csv",
# 			aniout=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs.{sampling}-aniout.csv",
# 			representative_list=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.txt",
# 			representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
# 			representative_lengths=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs_lengths.{sampling}.txt",
# 		message:
# 			"Creating vOUTs with Chechv aniclust"
# 		params:
# 			additional_reference_contigs=config['additional_reference_contigs']
# 		conda:
# 			dirs_dict["ENVS_DIR"] + "/vir.yaml"
# 		benchmark:
# 			dirs_dict["BENCHMARKS"] +"/vOUTclustering_references/{sampling}.tsv"
# 		threads: 64
# 		shell:
# 			"""
# 			cat {input.positive_contigs} {input.additional_reference_contigs} > {output.combined_positive_contigs}
# 			makeblastdb -in {output.combined_positive_contigs} -dbtype nucl -out {output.combined_positive_contigs}
# 			blastn -query {output.combined_positive_contigs} -db {output.combined_positive_contigs} -outfmt '6 std qlen slen' \
# 				-max_target_seqs 10000 -out {output.blastout} -num_threads {threads}
# 			python scripts/anicalc_checkv.py -i {output.blastout} -o {output.aniout}
# 			python scripts/aniclust_checkv.py --fna {output.combined_positive_contigs} --ani {output.aniout} --out {output.clusters} --min_ani 95 --min_tcov 85 --min_qcov 0
# 			cut {output.clusters} -f1 > {output.representative_list}
# 			seqtk subseq {output.combined_positive_contigs} {output.representative_list} > {output.representatives}
# 			cat {output.representatives} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
# 				$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.representative_lengths}
# 			"""

# 	rule filter_vOTUs_references:
# 		input:
# 			representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
# 			high_qualty_list=dirs_dict["vOUT_DIR"] + "/checkV_high_quality.{sampling}.txt",
# 			representative_lengths=dirs_dict["vOUT_DIR"] + "/viral_and_reference_contigs_lengths.tot.txt",
# 		output:
# 			filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
# 			filtered_list_temp=temp(dirs_dict["vOUT_DIR"]+ "/filtered_temp_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt"),
# 			filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
# 		params:
# 			min_votu_len=config['min_votu_length']
# 		message:
# 			"Filtering vOTUs "
# 		conda:
# 			dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 		benchmark:
# 			dirs_dict["BENCHMARKS"] +"/filter_vOTUs/{sampling}.tsv"
# 		threads: 1
# 		shell:
# 			"""
# 			cat {input.representative_lengths} | awk '$2>={params.min_votu_len}' | cut -f1 >>  {output.filtered_list_temp}
# 			cat {input.high_qualty_list} {output.filtered_list_temp}| sort | uniq > {output.filtered_list}
# 			seqtk subseq {input.representatives} {output.filtered_list} > {output.filtered_representatives}
# 			"""
