# ruleorder: vOUTclustering_references>vOUTclustering
# ruleorder: filter_vOTUs_references>filter_vOTUs

def input_vOTU_clustering(wildcards):
	input_list=[]
	if not NANOPORE_ONLY:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta",sample=SAMPLES))
	if NANOPORE:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_"+ LONG_ASSEMBLER + "_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta", sample=NANOPORE_SAMPLES))
	if CROSS_ASSEMBLY:
		input_list.append(dirs_dict["VIRAL_DIR"]+ "/ALL_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta")
	if SUBASSEMBLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_positive_" + VIRAL_ID_TOOL + ".{{sampling}}.fasta", sample=SAMPLES, subsample=subsample_test))
	if len(config['additional_reference_contigs'])>0:
		input_list.append(config['additional_reference_contigs'])
	return input_list

# if len(config['additional_reference_contigs'])==0:

rule derreplicate_assembly:
	input:
		positive_contigs=input_vOTU_clustering
	output:
		combined_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		derreplicated_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.fasta",
		derreplicated_clusters=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_cluster.tsv",
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
		mmseqs easy-cluster --threads {threads} --createdb-mode 1 --min-seq-id 1 -c 1 --cov-mode 1 {output.combined_positive_contigs} {params.rep_name} {params.rep_temp} 
		mv {params.rep_name_full} {output.derreplicated_positive_contigs}
		"""

rule vOUTclustering:
	input:
		fasta="{basedir}/{sequence}.fasta",
	output:
		clusters="{basedir}/{sequence}_95-85.clstr",
		blastout="{basedir}/{sequence}-blastout.csv",
		aniout="{basedir}/{sequence}-aniout.csv",
	message:
		"Creating vOUTs with CheckV aniclust"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	#  benchmark:
	#		dirs_dict['BENCHMARKS']+ "/vOUTclustering/{sequence}.tsv",
	threads: 144
	wildcard_constraints:
		sequence="[^/]+"  # The 'sequence' wildcard cannot contain a slash
	shell:
		"""
		makeblastdb -in {input.fasta} -dbtype nucl -out {input.fasta}
		blastn -query {input.fasta} -db {input.fasta} -outfmt '6 std qlen slen' \
				-max_target_seqs 10000000 -out {output.blastout} -num_threads {threads}
		python scripts/anicalc_checkv.py  -i {output.blastout} -o {output.aniout}
		python scripts/aniclust_checkv.py --fna {input.fasta} --ani {output.aniout} --out {output.clusters} --min_ani 95 --min_tcov 85 --min_qcov 0
		"""

def input_getHighQuality(wildcards):
	input_list=[]
	if not NANOPORE_ONLY:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{{sampling}}/quality_summary.tsv",sample=SAMPLES)),
	if NANOPORE:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/nanopore_{sample}_" + LONG_ASSEMBLER + "_checkV_{{sampling}}/quality_summary.tsv", sample=NANOPORE_SAMPLES)),
	if CROSS_ASSEMBLY:
		input_list.append(dirs_dict["vOUT_DIR"] + "/ALL_checkV_{sampling}/quality_summary.tsv"),
	if SUBASSEMBLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_" + VIRAL_ID_TOOL + "_checkV_{{sampling}}/quality_summary.tsv", sample=SAMPLES, subsample=subsample_test)),
	if len(config['additional_reference_contigs'])>0:
		input_list.append(dirs_dict["vOUT_DIR"] + "/user_reference_contigs_checkV/quality_summary.tsv"),
	return input_list

rule getHighQuality:
	input:
		input_getHighQuality,
	output:
		quality_summary_concat=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		high_qualty_list=dirs_dict["vOUT_DIR"] + "/checkV_high_quality.{sampling}.txt",
	message:
		"Getting list high-quality vOTUs"
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

checkpoint getHighQuality_clusters_fasta:
	input:
		new_clusters = dirs_dict["vOUT_DIR"] + "/new_references_clusters.{sampling}.csv",
		high_quality_list = dirs_dict["vOUT_DIR"] + "/checkV_high_quality.{sampling}.txt",
		combined_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		complete_clusters = dirs_dict["vOUT_DIR"] + "/new_references_complete_clusters.{sampling}.csv",
		fasta_dir = directory(dirs_dict["vOUT_DIR"] + "/high_quality_fastas.{sampling}")
	message:
		"Filtering clusters and extracting high-quality vOTUs FASTA sequences"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/filter_vOTUs/{sampling}.tsv"
	threads: 1
	shell:
		"""
		mkdir -p {output.fasta_dir}
		awk 'NR==FNR {{a[$1]; next}} $1 in a' {input.high_quality_list} {input.new_clusters} > {output.complete_clusters}
		awk '{{print $2 >> "{output.fasta_dir}/" $1 ".list"}}' {output.complete_clusters}

		for listfile in {output.fasta_dir}/*.list; do
			rep=$(basename "$listfile" .list)
			n_lines=$(wc -l < "$listfile" | tr -d '[:space:]')  # make sure to strip spaces
			if [ "$n_lines" -gt 1 ]; then
				seqtk subseq {input.combined_positive_contigs} "$listfile" > {output.fasta_dir}/"$rep".fasta
			fi
			rm "$listfile"
			 {input.results_dir_taxmyphage}/Results_per_genome/${{rep}}/query.fasta
		done
		"""

rule combine_with_taxmyphage:
    input:
        ref_fasta = "{contigs}.fasta",
        tax_fasta = lambda wildcards: dirs_dict["ANNOTATION"] + f"/taxmyphage_tot/Results_per_genome/{os.path.basename(wildcards.contigs)}/query.fasta"
    output:
        combined = "{contigs}_with_references.fasta"
    message:
        "Combining {input.ref_fasta} with taxmyphage result: {input.tax_fasta}"
    shell:
        """
        cat {input.ref_fasta} {input.tax_fasta} > {output.combined}
        """


rule select_vOTU_representative:
	input:
		merged_summary=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		cluster_file=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}_95-85.clstr",
		derreplicated_clusters=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_cluster.tsv",
	output:
		representatives=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list.{sampling}.csv",
		checkv_categories=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list_checkv_per_category.{sampling}.csv",
		new_clusters=dirs_dict["vOUT_DIR"]+ "/new_references_clusters.{sampling}.csv"
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
		combined_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
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
	threads: 1
	shell:
		"""
		seqtk subseq {input.combined_positive_contigs} {input.representative_list} > {output.representatives}
		cat {output.representatives} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
			$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.representative_lengths}
		"""
		
rule get_list_filtered_vOTUs:
	input:
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.{sampling}.csv",
		vOTUs_prefiltered=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		merged_summary=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		vibrant_circular=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_circular.{sampling}.csv",
		vibrant_positive=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_list.{sampling}.csv",
		vibrant_quality=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_quality.{sampling}.csv",
		vibrant_summary=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_summary_results.{sampling}.csv",
		virsorter_table=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/final-viral-score.tsv",
		virsorter_positive_list=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/positive_VS_list_{sampling}.txt",	
		genomad_virus_summary=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_virus_summary.tsv",
		genomad_plasmid_summary=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_plasmid_summary.tsv",
		genomad_viral_fasta=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/formatted_viral_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",										
		genomad_viral_fasta_conservative=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/formatted_viral_" + REPRESENTATIVE_CONTIGS_BASE + "_conservative.{sampling}.fasta",
		map_unfiltered=expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_unfiltered_contigs.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT),

	output:
		summary=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_summary.{sampling}.csv",
		filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
	params:
		samples=SAMPLES,
		contig_dir=dirs_dict["ASSEMBLY_DIR"],
		viral_dir=dirs_dict['VIRAL_DIR'],
		mapping_dir=dirs_dict['MAPPING_DIR'],
		subassembly=SUBASSEMBLY,
		cross_assembly=CROSS_ASSEMBLY,
		min_votu_len=config['min_votu_length'],
		key_samples=SAMPLES_key,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/05_vOTU_filtering.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/05_vOTU_filtering.py.ipynb"

rule filter_vOTUs:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
	output:
		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		min_votu_len=config['min_votu_length']
	message:
		"Filtering vOTUs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/filter_vOTUs/{sampling}.tsv"
	threads: 2
	shell:
		"""
		seqtk subseq {input.representatives} {input.filtered_list} > {output.filtered_representatives}
		"""

rule clustered_with_filter_vOTUs:
	input:
		derreplicated_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.fasta",
		new_clusters=dirs_dict["vOUT_DIR"]+ "/new_references_clusters.{sampling}.csv",
		filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
	output:
		cluster_filtered_representatives_list=dirs_dict["vOUT_DIR"]+ "/viral_contigs_clustered_with_filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
		cluster_filtered_representatives_fasta=dirs_dict["vOUT_DIR"]+ "/viral_contigs_clustered_with_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		min_votu_len=config['min_votu_length']
	message:
		"Filtering vOTUs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/filter_vOTUs/cluster_filtered_representatives_{sampling}.tsv"
	threads: 2
	shell:
		"""
		grep -f {input.filtered_list} {input.new_clusters} | cut -f2 > {output.cluster_filtered_representatives_list}
		seqtk subseq {input.derreplicated_positive_contigs} {output.cluster_filtered_representatives_list} > {output.cluster_filtered_representatives_fasta}
		"""


