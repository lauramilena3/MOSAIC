# ruleorder: gbk_to_faa>getORFs_assembly
def input_representative(wildcards):
	if METAGENOME:
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


rule getORFs_prodigal_gv:
	input:
		nuc_fasta="{fasta}.{sampling}.fasta",
	output:
		coords="{fasta}_ORFs.{sampling}.coords",
		aa="{fasta}_ORFs.{sampling}.faa",
	message:
		"Calling ORFs with prodigal-gv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		python ./scripts/parallel-prodigal-gv.py -q -i {input.nuc_fasta} -o {output.coords} -a {output.aa} -t {threads}
		"""

rule getORFs_coding_length:
	input:
		aa="{fasta}_ORFs.{sampling}.faa",
	output:
		length_temp=temp("{fasta}_ORFs_length_temp.{sampling}.txt"),
		names_temp=temp("{fasta}_ORFs_names_temp.{sampling}.txt"),
		length="{fasta}_ORFs_length.{sampling}.txt",
		cummulative_length="{fasta}_ORFs_coding_lengths.{sampling}.txt",
	message:
		"Calculating ORFs length"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		awk '$0 ~ ">" {{if (NR > 1) print c; c = 0; printf substr($0,2,100) "\t"}} $0 !~ ">" {{c += length($0)}} END {{print c}}' {input.aa} > {output.length_temp}
		cat {output.length_temp} | cut -f1 | cut -f1 -d' ' | rev | cut -d_ -f2- | rev > {output.names_temp}
		paste {output.length_temp} {output.names_temp} > {output.length}
		cat {output.length} | cut -f3,2 | awk '{{arr[$2]+=$1}} END {{for (i in arr) {{print i,arr[i]}}}}' > {output.cummulative_length}
		"""

rule clusterTaxonomy:
	input:
		aa=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs.{sampling}.faa",
		clusterONE_dir=config["clusterONE_dir"],
		gene2genome_millard=("db/vcontact2/6Nov2024_vConTACT2_gene_to_genome.csv"),
		vcontact_aa_millard=("db/vcontact2/6Nov2024_vConTACT2_proteins.faa"),
	output:
		gene2genome=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/gene2genome.csv",
		merged_gene2genome=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/gene2genome_merged.csv",
		merged_ORFs=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs_merged.{sampling}.fasta",
		genome_file=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/genome_by_genome_overview.csv",
		viral_cluster_overview=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/viral_cluster_overview.csv",
	params:
		out_dir=directory(dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}"),
		reference_genomes=config["reference_genomes_vcontact"],
		#reference_genomes='ProkaryoticViralRefSeq94-Merged'
	message:
		"Clustering viral genomes with vContact2"
	conda:
		dirs_dict["ENVS_DIR"] + "/vcontact.yaml"
	# conda:
	# 	dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/clusterTaxonomy/{sampling}.tsv"
	threads: 64
	shell:
		"""
		rm -rf {params.out_dir}
		mkdir {params.out_dir}
		vcontact2_gene2genome -p {input.aa} -s Prodigal-FAA -o {output.gene2genome}
		cat {output.gene2genome} {input.gene2genome_millard}  > {output.merged_gene2genome}
		cat {input.aa} {input.vcontact_aa_millard} > {output.merged_ORFs}
		dos2unix {output.merged_gene2genome}
		vcontact2 --raw-proteins {output.merged_ORFs} --rel-mode 'Diamond' --proteins-fp {output.merged_gene2genome} \
		--db {params.reference_genomes} --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {input.clusterONE_dir}/cluster_one-1.0.jar \
		--output-dir {params.out_dir} --threads {threads} -f
		#|| true
		"""

rule parseVcontact:
	input:
		viral_cluster_overview=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/viral_cluster_overview.csv",
		formatting_taxonomy_affiliations=config["taxonomy_file"],
	output:
		taxonomy_results=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vcontact2_taxonomy.{sampling}.csv",
	message:
		"Assigning viral taxonomy with vContact2 results"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/parseVcontact/{sampling}.tsv"
	threads: 1
	run:
		import pandas as pd
		def is_unique(s):
			a = s.to_numpy()
			return (a[0] == a).all()

		taxonomy_df=pd.read_csv(input.formatting_taxonomy_affiliations)
		df=pd.read_csv(input.viral_cluster_overview, index_col=0)

		df[['cluster','subcluster']] = df.VC.str.rsplit('_', 1, expand=True)
		grouped_df=df.groupby('cluster')
		grouped_results_df=pd.DataFrame()
		members=[]
		vc=[]
		for name, group in grouped_df:
			members.append(','.join([str(elem) for elem in group.Members.tolist()]))
			vc.append(name)
		grouped_results_df["Members"]=members
		grouped_results_df["VC"]=vc


		df=grouped_results_df[grouped_results_df['Members'].str.contains("NODE|tig0")]
		df["Members"]=df["Members"].str.split(",")
		accessions=[]
		nodes=[]
		taxonomies=[]

		with open(output.taxonomy_results, 'w') as f:
			for index, row in df.iterrows():
				accession1=([x for x in row['Members'] if not 'NODE' in x])
				accession=[item for item in accession1 if not item.startswith("tig00")]
				accession=([x for x in accession if not '~' in x])
				accessions.append(accession)
				#node=[x for x in row['Members'] if 'NODE' in x]
				node1=([x for x in row['Members'] if 'NODE' in x])
				node2=([x for x in row['Members'] if 'tig00' in x])
				node=node1+node2
				print(node)
				nodes.append(node)
				taxonomy=[]
				for acc in accession:
					print(acc)
					print(taxonomy_df[taxonomy_df["acc"]==acc]["lineage"].values[0].split(";"))
					taxonomy.append(taxonomy_df[taxonomy_df["acc"]==acc]["lineage"].values[0].split(";"))
				if taxonomy:
					tax_df=pd.DataFrame(taxonomy)
					tax_df.columns=["kindom", "phylum", "class", "order", "family", "genus", "species"]
					#print(tax_df)
					tax_df=tax_df.drop(columns="species")
					consensus_tax=""
					for (columnName, columnData) in tax_df.iteritems():
						#print('Colunm Name : ', columnName)
						if (is_unique(tax_df[columnName])):
							if not tax_df[columnName].to_numpy()[0] =="__":
								consensus_tax=(columnName, tax_df[columnName].to_numpy()[0])
					#print(tax_df)
					for n in node:
						print(n + "\t" + consensus_tax[1] +" [" + consensus_tax[0] + "]", file=f)

# rule mmseqsTaxonomy:
# 	input:
# 		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
# 		mmseqs_dir=(os.path.join(workflow.basedir, config['mmseqs_dir'])),
# 		refseq=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna")),
# 		refseq_taxid=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna.taxidmapping")),
# 	output:
# 		mmseqsdir=directory(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/"),
# 		html=(dirs_dict["vOUT_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.html"),
# 		tsv=(dirs_dict["vOUT_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tsv"),
# 		table=(dirs_dict["vOUT_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tbl"),
# 	message:
# 		"Taxonomy Assignment with MMseqs2"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env4.yaml"
# 	params:
# 		tmp=(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/tmp"),
# 		positive_contigsDB=(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/positive_contigsDB"),
# 		taxonomyResultDB=(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/taxonomyResultDB"),
# 		taxdump=(os.path.join(workflow.basedir,"db/ncbi-taxdump/")),
# 		refDB=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fnaDB")),
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env4.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/mmseqsTaxonomy/{sampling}.tsv"
# 	threads: 8
# 	shell:
# 		"""
# 		#analyse
# 		mkdir -p {output.mmseqsdir}
# 		mmseqs createdb {input.filtered_representatives} {params.positive_contigsDB}
# 		mmseqs taxonomy --threads {threads} {params.positive_contigsDB} {params.refDB} \
# 			{params.taxonomyResultDB} {params.tmp} --search-type 3 --lca-mode 2 -c 0.3 --cov-mode 2
# 		#results
# 		mmseqs createtsv {params.positive_contigsDB} {params.taxonomyResultDB} {output.tsv}
# 		mmseqs taxonomyreport {params.refDB} {params.taxonomyResultDB} {output.table}
# 		mmseqs taxonomyreport {params.refDB} {params.taxonomyResultDB} {output.html} --report-mode 1
# 	 	"""

rule PhaGCNTaxonomy:
	input:
		PhaGCN_newICTV_dir=config['PhaGCN_newICTV_dir'],
		fasta=dirs_dict["vOUT_DIR"] + "/{sequence}.fasta",
	output:
		taxonomy_table=dirs_dict["ANNOTATION"] + "/PhaGCN_taxonomy_report_{sequence}.csv",
	message:
		"Taxonomy Assignment with PhaGCN"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	params:
		taxonomy_table_temp=("final_prediction.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/PhaGCN_Taxonomy/{sequence}.tsv"
	threads: 32
	wildcard_constraints:
		  sequence="[^/]+"  # The 'sequence' wildcard cannot contain a slash
	shell:
		"""
		#analyse
		cd {input.PhaGCN_newICTV_dir}
		python run_Speed_up.py --contigs {input.fasta} --len 2000 --threads {threads}
		mv {params.taxonomy_table_temp} {output.taxonomy_table}
	 	"""

rule hostID_iphop:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/{sequence}.fasta",
		iphop_db=(config['iphop_db']),
	output:
		results_dir=directory(dirs_dict["ANNOTATION"] + "/iphop_hostID_{sequence}_resultsDir"),
	message:
		"Host finding with iphop"
	conda:
		dirs_dict["ENVS_DIR"] + "/env2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/iphop/{sequence}.tsv"
	threads: 64
	wildcard_constraints:
		  sequence="[^/]+"  # The 'sequence' wildcard cannot contain a slash
	shell:
		"""
		iphop predict --fa_file {input.fasta} --db_dir {input.iphop_db}/Aug_2023_pub_rw --out_dir {output.results_dir} --num_threads {threads}
		rm -rf {output.results_dir}/Wdir
		"""

rule single_fasta_filtered:
	input:
		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		filtered_representatives_checkpoint=temp(dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/filtered_representatives_dir_checkpoint.txt"),
	params:
		filtered_representatives_dir=((dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}")),
	message:
		"formating filtered vOTUs into single fasta"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	threads: 1
	shell:
		"""
		seqkit split --quiet -i {input.filtered_representatives} --out-dir {params.filtered_representatives_dir}
		touch {output.filtered_representatives_checkpoint}
	 	"""

rule match_spacers:
	input:
		spacers=config['microbial_spacers'],
		filtered_representatives_checkpoint=(dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/filtered_representatives_dir_checkpoint.txt"),
	output:
		spacer_match=dirs_dict["ANNOTATION"] + "/spacepharer_minced_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tsv",
	message:
		"Matching microbial spacers"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	params:
		filtered_representatives_dir=((dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}")),
		viralTargetDB=temp(directory(dirs_dict["ANNOTATION"] + "/viralTargetDB.{sampling}")),
		viralTargetDB_rev=temp(directory(dirs_dict["ANNOTATION"] + "/viralTargetDB_rev.{sampling}")),
		spacers_mincedSetDB=temp(directory(dirs_dict["ANNOTATION"] + "/spacers_mincedSetDB.{sampling}")),
		tmpFolder=temp(directory(dirs_dict["ANNOTATION"] + "/tmpFolder.{sampling}")),	
	benchmark:
		dirs_dict["BENCHMARKS"] +"/spacepharer/{sampling}.tsv"
	threads: 1
	shell:
		"""
		ulimit -S -s unlimited
		rm -rf {params.spacers_mincedSetDB}* {params.viralTargetDB}* {params.viralTargetDB_rev}* {params.tmpFolder}*
		spacepharer createsetdb {params.filtered_representatives_dir}/*fasta {params.viralTargetDB} {params.tmpFolder}
		spacepharer createsetdb {params.filtered_representatives_dir}/*fasta {params.viralTargetDB_rev} {params.tmpFolder} --reverse-fragments 1
		spacepharer createsetdb {input.spacers} {params.spacers_mincedSetDB} {params.tmpFolder} --extractorf-spacer 1
		spacepharer predictmatch {params.spacers_mincedSetDB} {params.viralTargetDB} {params.viralTargetDB_rev} {output.spacer_match} {params.tmpFolder} -s 7.5 
		rm -rf {params.spacers_mincedSetDB}* {params.viralTargetDB}* {params.viralTargetDB_rev}* {params.tmpFolder}*
	 	"""

rule match_spacers_dion:
	input:
		spacers_dion_db=(os.path.join(workflow.basedir, config['dion_db'])),
		filtered_representatives_dir=dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	output:
		spacer_match=dirs_dict["ANNOTATION"] + "/spacepharer_dion_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tsv",
		spacer_match_header=dirs_dict["ANNOTATION"] + "/spacepharer_dion_" + REPRESENTATIVE_CONTIGS_BASE + "_header.{sampling}.tsv",
	message:
		"Matching microbial spacers with the DION database"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	params:
		viralTargetDB=temp(directory(dirs_dict["ANNOTATION"] + "/viralTargetDB_dion.{sampling}")),
		viralTargetDB_rev=temp(directory(dirs_dict["ANNOTATION"] + "/viralTargetDB_dion_rev.{sampling}")),
		# spacers_dionSetDB=temp(directory(dirs_dict["ANNOTATION"] + "/spacers_dionSetDB.{sampling}")),
		tmpFolder=temp(directory(dirs_dict["ANNOTATION"] + "/tmpFolder_dion.{sampling}")),	
	benchmark:
		dirs_dict["BENCHMARKS"] +"/spacepharer/{sampling}_dion.tsv"
	threads: 1
	shell:
		"""
		rm -rf {params.viralTargetDB}* {params.viralTargetDB_rev}* {params.tmpFolder}*	 	
		spacepharer createsetdb {input.filtered_representatives_dir}/*fasta {params.viralTargetDB} {params.tmpFolder}
		spacepharer createsetdb {input.filtered_representatives_dir}/*fasta {params.viralTargetDB_rev} {params.tmpFolder} --reverse-fragments 1
		spacepharer predictmatch {input.spacers_dion_db}/dionSetDB {params.viralTargetDB} {params.viralTargetDB_rev} {output.spacer_match} {params.tmpFolder} -s 7.5 
		grep "^#" {output.spacer_match} > {output.spacer_match_header}
		rm -rf {params.viralTargetDB}* {params.viralTargetDB_rev}* {params.tmpFolder}*	 	
		"""


rule taxmyphage:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/{sequence}.fasta",
		taxmyphage_db=(config['taxmyphage_db']),
	output:
		results_dir=directory(dirs_dict["ANNOTATION"] + "/taxmyphage_{sequence}"),
	params:
		results_dir=directory(dirs_dict["ANNOTATION"] + "/taxmyphage_results"),
		annotation_dir=directory(dirs_dict["ANNOTATION"]),
	message:
		"Assigning taxonomy with taxmyphage"
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/taxmyphage/{sequence}.tsv"
	threads: 32
	shell:
		"""
		taxmyphage run -i {input.fasta} -t {threads} -db {input.taxmyphage_db} -o {params.results_dir}
		mv {params.results_dir} {output.results_dir}
		"""