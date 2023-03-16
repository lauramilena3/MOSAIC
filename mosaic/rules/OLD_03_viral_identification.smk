rule virSorter:
	input:
		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
		virSorter_dir=config['virSorter_dir'],
		virSorter_db=config['virSorter_db']
	output:
		out=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/VIRSorter_global-phage-signal.csv",
		results=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/virSorterCategories.txt"
	params:
		out_folder=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}"
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 32
	shell:
		"""
		WD=$(pwd)
		{config[virSorter_dir]}/wrapper_phage_contigs_sorter_iPlant.pl -f {input.merged_assembly} \
			--db 2 \
			--wdir {params.out_folder} \
			--ncpu {threads} \
			--data-dir $WD/{input.virSorter_db} \
			--virome
		grep ">" {params.out_folder}/Predicted_viral_sequences/VIRSorter_*fasta > {output.results}
		"""

# rule virFinder:
# 	input:
# 		representatives=dirs_dict["vOUT_DIR"] + "/representative_contigs.{sampling}.fasta",
# 		virFinder_dir=config['virFinder_dir'],
# 	output:
# 		pvalues=dirs_dict["VIRAL_DIR"] + "/virFinder_pvalues.{sampling}.txt"
# 	params:
# 		virFinder_script="scripts/virfinder_wrapper.R"
# 	message:
# 		"Scoring virus VirFinder"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/vir.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		Rscript {params.virFinder_script} {input.representatives} {output.pvalues}
# 		"""

rule viralID_VIBRANT:
	input:
		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
		VIBRANT_dir=os.path.join(workflow.basedir, config['vibrant_dir']),
	output:
		plus5000_list=dirs_dict["VIRAL_DIR"] + "/merged_scaffolds_over5000.{sampling}.txt",
		plus5000_contigs=dirs_dict["VIRAL_DIR"] + "/merged_scaffolds_over5000.{sampling}.fasta",
		vibrant=directory(dirs_dict["VIRAL_DIR"] + "/VIBRANT_merged_scaffolds_over5000.{sampling}"),
		vibrant_circular=dirs_dict["VIRAL_DIR"] + "/VIBRANT_merged_scaffolds_circular.{sampling}.txt",
		vibrant_positive=dirs_dict["VIRAL_DIR"] + "/VIBRANT_merged_scaffolds_positive_list.{sampling}.txt",
		vibrant_quality=dirs_dict["VIRAL_DIR"] + "/VIBRANT_merged_scaffolds_positive_quality.{sampling}.txt",
	params:
		viral_dir=directory(dirs_dict["VIRAL_DIR"]),
		minlen=5000,
		name_circular=dirs_dict["VIRAL_DIR"] + "/VIBRANT_merged_scaffolds_over5000.{sampling}/VIBRANT_results*/VIBRANT_complete_circular*.{sampling}.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	message:
		"Identifying viral contigs with VIBRANT"
	threads: 8
	shell:
		"""
		cd {params.viral_dir}
		cat {input.merged_assembly} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' | \
			awk  '$2 > {params.minlen}' | cut -f1 > {output.plus5000_list}
		head {output.plus5000_list}
		seqtk subseq {input.merged_assembly} {output.plus5000_list} > {output.plus5000_contigs}
		{input.VIBRANT_dir}/VIBRANT_run.py -i {output.plus5000_contigs} -t {threads}
		cut -f1 {params.name_circular} > {output.vibrant_circular}
		touch {output.vibrant_circular}
		cp {output.vibrant}/VIBRANT_phages_*/*phages_combined.txt {output.vibrant_positive}
		cp {output.vibrant}/VIBRANT_results*/VIBRANT_genome_quality*.tsv {output.vibrant_quality}
		#vibrant_figures=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_figures_" +VIRAL_CONTIGS_BASE + ".tot"),
		#vibrant_tables_parsed=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_parsed_" +VIRAL_CONTIGS_BASE + ".tot"),
		#vibrant_tables_unformated=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_unformatted_" +VIRAL_CONTIGS_BASE + ".tot"),
		#vibrant_phages=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_unformatted_" +VIRAL_CONTIGS_BASE + ".tot"),
		#vibrant_results=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_unformatted_" +VIRAL_CONTIGS_BASE + ".tot"),
		"""

rule parseViralTable:
	input:
		categories=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/virSorterCategories.txt",
		vibrant=directory(dirs_dict["VIRAL_DIR"] + "/VIBRANT_merged_scaffolds_over5000.{sampling}"),
		vibrant_circular=dirs_dict["VIRAL_DIR"] + "/VIBRANT_merged_scaffolds_circular.{sampling}.txt",
		vibrant_positive=dirs_dict["VIRAL_DIR"] + "/VIBRANT_merged_scaffolds_positive_list.{sampling}.txt",
		vibrant_quality=dirs_dict["VIRAL_DIR"] + "/VIBRANT_merged_scaffolds_positive_quality.{sampling}.txt",
		merged_assembly_len=dirs_dict["VIRAL_DIR"] + "/merged_scaffolds_lengths.{sampling}.txt",
	output:
		circular_unk=dirs_dict["VIRAL_DIR"]+ "/unknown_circular_list.{sampling}.txt",
		table=dirs_dict["VIRAL_DIR"]+ "/viral_table.{sampling}.csv",
		positive_list=dirs_dict["VIRAL_DIR"]+ "/positive_VS_VB_list.{sampling}.txt",
	message:
		"Parsing VirSorter and VIBRANT results"
	run:
		import pandas as pd
		VS=input.categories
		VB_circular = input.vibrant_circular
		lenghts=input.merged_assembly_len
		above5000=input.vibrant_positive
		qualityAbove=input.vibrant_quality

		circular_unk=output.circular_unk
		positive=output.positive_list
		viral_table=output.table
		#VirSorter
		nameVS=[]
		catVS=[]
		circularVS=[]
		typeVS=[]

		with open(VS) as fp:
			line = fp.readline()
			cnt = 1
			while line:
				circular="N"
				contig=line.split(">")[1]
				contigName=contig.split("-")[0].split("VIRSorter_")[1].replace(".", "_")
				category=float(contig.split("-")[-1].split("_")[1])
				if category<4:
					contig_type="Complete phage contigs"
				else:
					contig_type="Prophage"
					category=category-3
					contigName=contigName.split("_gene")[0]
				if "-circular" in contig:
					circular="Y"
				nameVS.append(contigName)
				catVS.append(category)
				circularVS.append(circular)
				typeVS.append(contig_type)
				line = fp.readline()
				cnt += 1

		resultsVS=pd.DataFrame()
		resultsVS["name"]=nameVS
		resultsVS["catVS"]=catVS
		resultsVS["circularVS"]=circularVS
		resultsVS["typeVS"]=typeVS
		resultsVS=resultsVS.drop_duplicates(subset=['name'], keep="first")

		#VIBRANT
		nameVB=[]
		circular_VB=[]

		with open(VB_circular) as fp:
			line = fp.readline()
			cnt = 1
			while line:
				contigName=line.replace(".", "_").rstrip("\n\r")
				nameVB.append(contigName)
				circular_VB.append("Y")
				line = fp.readline()
				cnt += 1

		circularVB=pd.DataFrame()
		circularVB["nameVB"]=nameVB
		circularVB["circularVB"]=circular_VB

		nameVIBRANT=[]
		VB_positive=[]

		with open(above5000) as fp:
			line = fp.readline()
			cnt = 1
			while line:
				if cnt != 1:
					nameVIBRANT.append(line.rstrip("\n\r").replace(".", "_"))
					VB_positive.append(1)
				line = fp.readline()
				cnt += 1

		resultsVIBRANT=pd.DataFrame()

		resultsVIBRANT["nameVIBRANT"]=nameVIBRANT
		resultsVIBRANT["VB_positive"]=VB_positive
		resultsVIBRANT


		nameVIBRANT=[]
		VB_quality=[]

		with open(qualityAbove) as fp:
			line = fp.readline()
			cnt = 1
			while line:
				if cnt != 1:
					quality=line.rstrip("\n\r").replace(".", "_").split("\t")[2]
					if not quality == "complete circular":
						nameVIBRANT.append(line.rstrip("\n\r").replace(".", "_").split("\t")[0].split("_fragment")[0])
						VB_quality.append(line.rstrip("\n\r").replace(".", "_").split("\t")[2].split(" ")[0])
				line = fp.readline()
				cnt += 1


		qualityVIBRANT=pd.DataFrame()

		qualityVIBRANT["nameVIBRANT"]=nameVIBRANT
		qualityVIBRANT["VB_quality"]=VB_quality
		qualityVIBRANT

		#LENGTH
		nameLEN=[]
		LEN=[]

		with open(lenghts) as fp:
			line = fp.readline()
			cnt = 1
			while line:
				if cnt != 1:
					nameLEN.append(line.split("\t")[0].rstrip("\n\r").replace(".", "_"))
					LEN.append(int(line.split("\t")[1].rstrip("\n\r")))
				line = fp.readline()
				cnt += 1

		resultsLEN=pd.DataFrame()
		resultsLEN["nameLEN"]=nameLEN
		resultsLEN["LEN"]=LEN
		resultsLEN

		resultsVS=resultsVS.drop_duplicates(subset=['name'], keep="first")
		qualityVIBRANT=qualityVIBRANT.drop_duplicates(subset=['nameVIBRANT'], keep="first")

		new_df = pd.merge(resultsLEN, circularVB,  how='left', left_on=['nameLEN'], right_on = ['nameVB'])
		new_df2 = pd.merge(new_df, qualityVIBRANT,  how='left', left_on=['nameLEN'], right_on = ['nameVIBRANT'])
		new_df3 = pd.merge(new_df2, resultsVS,  how='left', left_on=['nameLEN'], right_on = ['name'])
		final = pd.merge(new_df3, resultsVIBRANT,  how='left', left_on=['nameLEN'], right_on = ['nameVIBRANT'])


		final['catVS'] = final['catVS'].fillna(value=int(4))
		final["catVS"] = final["catVS"].astype('float')
		final['circularVS']=final['circularVS'].fillna(value="NA")
		final['circularVB']=final['circularVB'].fillna(value="NA")
		final['LEN'] = final['LEN'].fillna(value=float(0))
		final["LEN"] = final["LEN"].astype('float')
		final['typeVS']=final['typeVS'].fillna(value="NA")
		final['VB_positive']=final['VB_positive'].fillna(value=float(0))

		final.loc[final['catVS'] <= float(2), 'VS_positive'] = float(1)
		final['VS_positive'] = final['VS_positive'].fillna(value=float(0))
		final.loc[(final['VS_positive'] == float(1)) & (final['VB_positive'] == float(1)), 'VS_VB'] = float(1)
		final.loc[((final['VS_positive'] == float(1)) | (final['VB_positive'] == float(1))), 'POSITIVE' ] = "Y"
		final.loc[((final['circularVB'] == "Y" )| (final['circularVS'] == "Y")), 'CIRCULAR' ] = "Y"
		final.loc[((final['POSITIVE'] =="Y") & (final['LEN'] >= float(5000))), 'VIRAL' ] = "Y"
		final.loc[((final["CIRCULAR"]=="Y") | (final["VB_quality"]=="high")| (final["VB_quality"]=="medium")) &  (final["LEN"]<float(5000))  &  (final["POSITIVE"]=="Y"), 'VIRAL' ] = "Y"

		viral=final[final["VIRAL"]=="Y"].nameLEN.tolist()
		corr=[]
		for s in viral:
			k = s.rfind("_")
			new_string = s[:k] + "." + s[k+1:]
			corr.append(new_string)
		f=open(positive, 'w')
		f.write("\n".join(corr))
		f.close()

		rep_check=final[(final["CIRCULAR"]=="Y") & (final["LEN"]<=float(7000)) & (final["VIRAL"]!="Y")].name.tolist()
		f=open(circular_unk, 'w')
		f.write("\n".join(rep_check))
		f.close()

		final.to_csv(viral_table, sep="\t")

rule hmmCircularContigs:
	input:
		circular_unk=dirs_dict["VIRAL_DIR"]+ "/unknown_circular_list.{sampling}.txt",
		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
	output:
		edited_fasta=(dirs_dict["vOUT_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
		circular_unk_fasta=dirs_dict["VIRAL_DIR"]+ "/unknown_circular.{sampling}.fna",
		coords=dirs_dict["VIRAL_DIR"] + "/unknown_circular.{sampling}.coords",
		aa=dirs_dict["VIRAL_DIR"] + "/unknown_circular.{sampling}.faa",
		hmm_results=dirs_dict["VIRAL_DIR"]+ "/hmm_parsed.{sampling}.out",
		hmm_out=dirs_dict["VIRAL_DIR"]+ "/hmmsearch.{sampling}.out",
		hmm_list=dirs_dict["VIRAL_DIR"]+ "/positive_rep_list.{sampling}.txt",
	params:
		hmm="db/hmm/ssDNA.hmm",
		min_score=50,
		min_eval=0.001,
	message:
		"Selecting Viral Circular Contigs with hmmsearch"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	shell:
		"""
		sed 's/\./_/g' {input.merged_assembly} > {output.edited_fasta}
		seqtk subseq {output.edited_fasta} {input.circular_unk} > {output.circular_unk_fasta}
		if [ -s {output.circular_unk_fasta} ]
		then
			prodigal -i {output.circular_unk_fasta} -o {output.coords} -a {output.aa} -p meta
			hmmsearch --tblout {output.hmm_out} -E {params.min_eval} {params.hmm} {output.aa}
			cat {output.hmm_out} | grep -v '^#' | awk '{{ if ( $6 > {params.min_score} ) {{print $1,$3,$5,$6}} }}' > {output.hmm_results} || true
			cut -d' ' -f1 {output.hmm_results} | sort | uniq > {output.hmm_list}
		else
			touch {output.hmm_out}
			touch {output.hmm_results}
			touch {output.hmm_list}
			touch {output.coords}
			touch {output.aa}
		fi
		"""
rule extractViralContigs:
	input:
		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
		positive_rep_list=dirs_dict["VIRAL_DIR"]+ "/positive_rep_list.{sampling}.txt",
		positive_VS_VB_list=dirs_dict["VIRAL_DIR"]+ "/positive_VS_VB_list.{sampling}.txt",
	output:
		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/formatted_merged_scaffolds.{sampling}.fasta"),
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		positive_rep_contigs=dirs_dict["VIRAL_DIR"]+ "/positive_rep_contigs.{sampling}.fasta",
		positive_VS_VB_contigs=dirs_dict["VIRAL_DIR"]+ "/positive_VS_VB_contigs.{sampling}.fasta",
	message:
		"Selecting Viral Contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	shell:
		"""
		sed 's/\./_/g' {input.merged_assembly} > {output.merged_assembly}
		sed -i 's/\./_/g' {input.positive_VS_VB_list}
		seqtk subseq {output.merged_assembly} {input.positive_rep_list} > {output.positive_rep_contigs}
		seqtk subseq {output.merged_assembly} {input.positive_VS_VB_list} > {output.positive_VS_VB_contigs}
		cat {output.positive_rep_contigs} {output.positive_VS_VB_contigs} > {output.positive_contigs}
		"""
