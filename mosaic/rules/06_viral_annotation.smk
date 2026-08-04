rule lifestyle_bacphlip:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/{sequence}.fasta",
	output:
		results_bacphlip_final=(dirs_dict["ANNOTATION"] + "/{sequence}_bacphlip.csv"),
	params:
		results_bacphlip=(dirs_dict["vOUT_DIR"] + "/{sequence}.fasta.bacphlip"),
		results_dir=((dirs_dict["vOUT_DIR"] + "/{sequence}.fasta.BACPHLIP_DIR")),
	message:
		"Predicting lifecycle with BACPHLIP"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/bacphlip/{sequence}.tsv"
	threads: 1
	wildcard_constraints:
		  sequence="[^/]+"  # The 'sequence' wildcard cannot contain a slash
	shell:
		"""
		rm -rf {params.results_dir}
		mkdir {params.results_dir}
		cd {params.results_dir}
		bacphlip -i {input.fasta} --multi_fasta -f
		mv {params.results_bacphlip} {output.results_bacphlip_final}
		rm -rf {params.results_dir}
		"""


rule estimateGenomeCompletness_long:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/{sample}_"+ LONG_ASSEMBLER + "_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		checkv_db=(config['checkv_db']),
	output:
		quality_summary=dirs_dict["vOUT_DIR"] + "/nanopore_{sample}_" + LONG_ASSEMBLER + "_checkV_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["vOUT_DIR"] + "/nanopore_{sample}_" + LONG_ASSEMBLER + "_checkV_{sampling}/completeness.tsv",
		contamination=dirs_dict["vOUT_DIR"] + "/nanopore_{sample}_" + LONG_ASSEMBLER + "_checkV_{sampling}/contamination.tsv",
	params:
		checkv_outdir=dirs_dict["vOUT_DIR"] + "/nanopore_{sample}_" + LONG_ASSEMBLER + "_checkV_{sampling}",
		tmp=((dirs_dict["vOUT_DIR"] + "/nanopore_{sample}_" + LONG_ASSEMBLER + "_checkV_{sampling}/tmp")),
	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/estimateGenomeCompletness/nanopore_{sample}_{sampling}.tsv"
	threads: 4
	shell:
		"""
		if [ -s {input.positive_contigs} ]; then
		    		            	checkv contamination {input.positive_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		    		            	checkv completeness {input.positive_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		    		            	checkv complete_genomes {input.positive_contigs} {params.checkv_outdir}
		    		            	checkv quality_summary {input.positive_contigs} {params.checkv_outdir}
										rm -rf {params.tmp}
		else
		    		            	echo "The FASTA file {input.positive_contigs} is empty"
		    		            	mkdir -p {params.checkv_outdir}
		    		            	touch {output.quality_summary}
		    		            	touch {output.completeness}
		    		            	touch {output.contamination}
		fi

		"""

rule estimateGenomeCompletness:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		checkv_db=(config['checkv_db']),
	output:
		quality_summary=dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{sampling}/completeness.tsv",
		contamination=dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{sampling}/contamination.tsv",
	params:
		checkv_outdir=dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{sampling}",
		tmp=((dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{sampling}/tmp")),
	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/estimateGenomeCompletness/{sample}_{sampling}.tsv"
	threads: 4
	shell:
		"""
		if [ -s {input.positive_contigs} ]; then
		    		            	checkv contamination {input.positive_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		    		            	checkv completeness {input.positive_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		    		            	checkv complete_genomes {input.positive_contigs} {params.checkv_outdir}
		    		            	checkv quality_summary {input.positive_contigs} {params.checkv_outdir}
										rm -rf {params.tmp}

		else
		    		            	echo "The FASTA file {input.positive_contigs} is empty"
		    		            	mkdir -p {params.checkv_outdir}
		    		            	touch {output.quality_summary}
		    		            	touch {output.completeness}
		    		            	touch {output.contamination}
		fi
		"""

rule estimateGenomeCompletness_reference:
	input:
		reference_contigs=config['additional_reference_contigs'],
		checkv_db=(config['checkv_db']),
	output:
		quality_summary=dirs_dict["vOUT_DIR"] + "/user_reference_contigs_checkV/quality_summary.tsv",
		completeness=dirs_dict["vOUT_DIR"] + "/user_reference_contigs_checkV/completeness.tsv",
		contamination=dirs_dict["vOUT_DIR"] + "/user_reference_contigs_checkV/contamination.tsv",
	params:
		checkv_outdir=dirs_dict["vOUT_DIR"] + "/user_reference_contigs_checkV",
		tmp=dirs_dict["vOUT_DIR"] + "/user_reference_contigs_checkV/tmp",
	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/estimateGenomeCompletness/user_reference_contigs_checkV.tsv"
	threads: 32
	shell:
		"""
		rm -rf {params.checkv_outdir} || true
		if [ -s {input.reference_contigs} ]; then
		    		            	checkv contamination {input.reference_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		    		            	checkv completeness {input.reference_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		    		            	checkv complete_genomes {input.reference_contigs} {params.checkv_outdir}
		    		            	checkv quality_summary {input.reference_contigs} {params.checkv_outdir}
										rm -rf {params.tmp}
		else
		    		            	echo "The FASTA file {input.reference_contigs} is empty"
		    		            	mkdir -p {params.checkv_outdir}
		    		            	touch {output.quality_summary}
		    		            	touch {output.completeness}
		    		            	touch {output.contamination}
		fi
		"""

rule virSorter2_DRAM:
	input:
		# cluster_filtered_representatives_fasta=dirs_dict["vOUT_DIR"]+ "/viral_contigs_clustered_with_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		virSorter_db=config['virSorter_db'],
	output:
		positive_fasta=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/final-viral-combined.fa",
		table_virsorter=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/final-viral-score.tsv",
		positive_list=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/positive_VS_list_{sampling}.txt",
		DRAM_tab=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/for-dramv/viral-affi-contigs-for-dramv.tab",
		DRAM_fasta=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/for-dramv/final-viral-combined-for-dramv.fa",
		iteration=directory(dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/iter-0"),
	params:
		out_folder=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}"
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/VirSorter2_DRAM/{sampling}_illumina.tsv"
	threads: 64
	shell:
		"""
		virsorter run -w {params.out_folder} -i {input.representatives} -j {threads} --db-dir {input.virSorter_db} \
				--include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --seqname-suffix-off  --provirus-off --min-length 0 \
				--viral-gene-enrich-off --prep-for-dramv --keep-original-seq --min-score 0
		grep ">" {output.positive_fasta} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list} || true
		"""

rule DRAMv_annotation:
	input:
		DRAM_tab=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/for-dramv/viral-affi-contigs-for-dramv.tab",
		DRAM_fasta=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/for-dramv/final-viral-combined-for-dramv.fa",		
		DRAM_db=config['DRAM_db'],
	output:
		DRAM_output=directory(dirs_dict["ANNOTATION"]+ "/vDRAM_annotate_results_{sampling}"),
	params:
		DRAM_annotations=dirs_dict["ANNOTATION"]+ "/vDRAM_annotate_results_{sampling}/annotations.tsv",
		# trna=directory(dirs_dict["vOUT_DIR"]+ "/DRAM_combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq_{sampling}/trnas.tsv"),
		# rrna=directory(dirs_dict["vOUT_DIR"]+ "/DRAM_combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq_{sampling}/rrnas.tsv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/DRAM/{sampling}.tsv"
	message:
		"Annotate contigs with DRAM"
	threads: 64
	shell:
		"""
		DRAM-v.py annotate -i {input.DRAM_fasta} -v {input.DRAM_tab} -o {output.DRAM_output} --threads 64
		"""

rule DRAMv_distill:
	input:
		DRAM_tab=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/for-dramv/viral-affi-contigs-for-dramv.tab",
		DRAM_fasta=dirs_dict["ANNOTATION"] + "/VirSorter2_DRAM_{sampling}/for-dramv/final-viral-combined-for-dramv.fa",		
		DRAM_db=config['DRAM_db'],
		DRAM_output=dirs_dict["ANNOTATION"]+ "/vDRAM_annotate_results_{sampling}",
	output:
		DRAM_summary=directory(dirs_dict["ANNOTATION"]+ "/vDRAM_distill_results_{sampling}"),
	params:
		DRAM_annotations=dirs_dict["ANNOTATION"]+ "/vDRAM_annotate_results_{sampling}/annotations.tsv",
		# trna=directory(dirs_dict["vOUT_DIR"]+ "/DRAM_combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq_{sampling}/trnas.tsv"),
		# rrna=directory(dirs_dict["vOUT_DIR"]+ "/DRAM_combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq_{sampling}/rrnas.tsv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/DRAM/{sampling}.tsv"
	message:
		"Annotate contigs with DRAM"
	threads: 64
	shell:
		"""
		DRAM-v.py distill -i {params.DRAM_annotations} -o {output.DRAM_summary} 
		"""

rule DRAMv_extract_genes:
	input:
		DRAM_output=dirs_dict["ANNOTATION"]+ "/vDRAM_annotate_results_{sampling}",
	output:
		mmseqs_all=temp((dirs_dict["ANNOTATION"]+ "/predicted_genes_95_85_{sampling}_all_seqs.fasta")),
		genes_fna=dirs_dict["ANNOTATION"]+ "/predicted_genes_{sampling}.fna",
		NR_fna=dirs_dict["ANNOTATION"]+ "/predicted_genes_NR_95_85_{sampling}.fna",
		NR_clusters=dirs_dict["ANNOTATION"]+ "/predicted_genes_NR_95_85_{sampling}_clusters.tsv",
		NR_fna_150=dirs_dict["ANNOTATION"]+ "/predicted_genes_NR_95_85_150bp_{sampling}.fna",
		NR_fna_150_list=dirs_dict["ANNOTATION"]+ "/predicted_genes_NR_95_85_150bp_list_{sampling}.txt",
		NR_fna_150_temp_gff=temp(dirs_dict["ANNOTATION"]+ "/predicted_genes_NR_95_85_150bp_{sampling}.gff_temp"),
		NR_fna_150_gff=dirs_dict["ANNOTATION"]+ "/predicted_genes_NR_95_85_150bp_{sampling}.gff",
	params:
		DRAM_fna=dirs_dict["ANNOTATION"]+ "/vDRAM_annotate_results_{sampling}/genes.fna",
		DRAM_gff=dirs_dict["ANNOTATION"]+ "/vDRAM_annotate_results_{sampling}/genes.gff",
		mmseqs_clusters=(dirs_dict["ANNOTATION"]+ "/predicted_genes_95_85_{sampling}_cluster.tsv"),
		mmseqs_temp=((dirs_dict["ANNOTATION"]+ "/temp_NR_mmseqs_{sampling}")),
		mmseqs_name="predicted_genes_95_85_{sampling}",
		annotation_dir=dirs_dict["ANNOTATION"],
	conda:
		dirs_dict["ENVS_DIR"] + "/vcontact.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/AMG/{sampling}.tsv"
	message:
		"Derreplicate genes with mmseqs"
	threads: 64
	shell:
		"""
		cd {params.annotation_dir}
		mmseqs easy-cluster --threads {threads} --createdb-mode 1 --min-seq-id 0.95 -c 0.85 --cov-mode 1 \
			{params.DRAM_fna} {params.mmseqs_name} {params.mmseqs_temp}
		rm -rf {params.mmseqs_temp}
		cp {params.DRAM_fna} {output.genes_fna}
		mv {params.mmseqs_clusters} {output.NR_clusters}
		mv {params.mmseqs_name}_rep_seq.fasta {output.NR_fna}
		seqtk seq -L 150 {output.NR_fna} > {output.NR_fna_150}
		grep ">" {output.NR_fna_150} | sed "s/>//g" > {output.NR_fna_150_list}
		awk -F'\t' '
		BEGIN {{ OFS = "\t" }}
		$3 == "CDS" {{
			split($9, a, ";")
			for (i in a) {{
				if (a[i] ~ /^ID=/) {{
						id = substr(a[i], 4)
						$1 = id                    # Set seqname = ID
						len = $5 - $4 + 1          # Calculate feature length
						$4 = 1                     # Set start to 1
						$5 = len                   # Set end to length
						$9 = "ID=" id              # Rewrite attributes to just ID
						print
				}}
			}}
		}}
		' {params.DRAM_gff} > {output.NR_fna_150_temp_gff}
    	awk 'NR==FNR {{ ids[$1]; next }} $1 in ids' {output.NR_fna_150_list} {output.NR_fna_150_temp_gff} > {output.NR_fna_150_gff}
		"""

rule pharokka_annotation_genbank:
	input:
		DRAM_output=dirs_dict["ANNOTATION"]+ "/vDRAM_annotate_results_{sampling}",
		pharokka_db = config["pharokka_db"]
	output:
		pharokka_output=directory(dirs_dict["ANNOTATION"]+ "/vDRAM_results_{sampling}_pharokka"),
	params:
		DRAM_gbk=dirs_dict["ANNOTATION"]+ "/vDRAM_annotate_results_{sampling}/genbank/final-viral-combined-for-dramv.gbk",
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	message:
		"Annotate contigs with pharokka"
	threads: 144
	shell:
		"""
		pharokka.py -i {params.DRAM_gbk} -o {output.pharokka_output} -d {input.pharokka_db} -t {threads} -m -f --genbank
		"""

rule pharokka_annotation:
	input:
		fasta="{contigs}.fasta",
		pharokka_db = config["pharokka_db"]
	output:
		pharokka_output = directory("{contigs}_pharokka")
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	message:
		"Annotate contigs with pharokka"
	threads: 16
	shell:
		"""
		pharokka.py -i {input.fasta} -o {output.pharokka_output} -d {input.pharokka_db} -t {threads} -m -f --dnaapler
		"""

rule annotate_VIGA:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/{contigs}.tot.fasta",
		VIGA_dir=os.path.join(workflow.basedir, config['viga_dir']),
		piler_dir=os.path.join(workflow.basedir, (config['piler_dir'])),
		trf_dir=os.path.join(workflow.basedir, (config['trf_dir'])),
	output:
		modifiers=temp(dirs_dict["ANNOTATION"] + "/modifiers_{contigs}.txt"),
		temp_symlink=temp(dirs_dict["ANNOTATION"] + "/{contigs}.tot.fasta_symlink"),
		temp_viga_dir=temp(directory(dirs_dict["ANNOTATION"] + "/{contigs}_tempVIGA")),
		GenBank_file=dirs_dict["ANNOTATION"] + "/{contigs}.tot_VIGA_annotated.gbk",
		GenBank_table_temp1=temp(dirs_dict["ANNOTATION"] + "/{contigs}.tot_annotated.tbl"),
		GenBank_table_temp2=temp(dirs_dict["ANNOTATION"] + "/{contigs}.tot_annotated.tbl2"),
		GenBank_table=dirs_dict["ANNOTATION"] + "/{contigs}.tot_VIGA_annotated.tbl",
		GenBank_fasta=dirs_dict["ANNOTATION"] + "/{contigs}.tot_VIGA__annotated.fasta",
		csv=dirs_dict["ANNOTATION"] + "/{contigs}.tot_VIGA_annotated.csv",
		viga_names=temp(dirs_dict["ANNOTATION"] + "/viga_names_{contigs}.tot.txt"),
		viga_topology_temp=temp(dirs_dict["ANNOTATION"] + "/viga_topology_temp_{contigs}_tot.txt"),
		viga_topology=(dirs_dict["ANNOTATION"] + "/viga_topology_{contigs}_tot.txt"),
	params:
		viga_log=dirs_dict["ANNOTATION"] + "/viga_log_{contigs}.tot.txt",
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/{contigs}",
		results_name=dirs_dict["MMSEQS"] + "/{contigs}_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
		VIGA_dir=directory("../" + config['viga_dir']),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_VIGA_{contigs}/tot.tsv"
	message:
		"Annotating contigs with VIGA"
	threads: 8
	shell:
		"""
		PILER={input.piler_dir}
		PATH=$PILER:$PATH
		TRF={input.trf_dir}
		PATH=$TRF:$PATH
		ln -sfn {input.representatives} {output.temp_symlink}
		mkdir -p {output.temp_viga_dir}
		cd {output.temp_viga_dir}
		touch {output.modifiers}
		echo "viga.1"
		{input.VIGA_dir}/VIGA.py --input {output.temp_symlink} --diamonddb {input.VIGA_dir}/databases/RefSeq_Viral_DIAMOND/refseq_viral_proteins.dmnd \
		--blastdb {input.VIGA_dir}/databases/RefSeq_Viral_BLAST/refseq_viral_proteins --hmmerdb {input.VIGA_dir}/databases/pvogs/pvogs.hmm \
		--rfamdb {input.VIGA_dir}/databases/rfam/Rfam.cm --modifiers {output.modifiers} --threads {threads} &> {params.viga_log}
		echo "viga.2"
		cat {params.viga_log} | grep "was renamed as" > {output.viga_names}
		echo "viga.3"
		cat {params.viga_log} | grep "according to LASTZ" > {output.viga_topology_temp}
		echo "viga.4"
		cat {output.viga_names} | while read line
		do
		    		            stringarray=($line)
		    		            new=${{stringarray[-1]}}
		    		            old=${{stringarray[1]}}
		    		            sed -i -e "s/${{new}}\t/${{old}}\t/g" -e "s/${{new}}_/${{old}}_/g" {output.csv}
		    		            sed -i -e "s/${{new}}$/${{old}}/g" -e "s/${{new}} /${{old}} /g" -e "s/${{new}}_/${{old}}_/g" {output.GenBank_file}
		    		            sed -i -e "s/${{new}}$/${{old}}/g" -e "s/${{new}} /${{old}} /g" -e "s/${{new}}_/${{old}}_/g" {output.GenBank_table_temp1}
		    		            sed -i "s/>${{new}} $/>${{old}}/g" {output.GenBank_fasta}
		    		            sed -i -e "s/${{new}} /${{old}} /g" {output.viga_topology_temp}
		done
		echo "viga.5"
		awk  '{{print $1 "\t" $6}}'  {output.viga_topology_temp} > {output.viga_topology}
		echo "viga.6"
		grep -v "gene$" {output.GenBank_table_temp1} > {output.GenBank_table_temp2}
		echo "viga.7"
		grep -n "CDS$" {output.GenBank_table_temp2} | cut -d : -f 1 | awk '{{$1+=-1}}1' | sed 's%$%d%' | sed -f - {output.GenBank_table_temp2} > {output.GenBank_table}
		echo "viga.8"
		sed -i "s/tRNA-?(Asp|Gly)(atcc)/tRNA-Xxx/g" {output.GenBank_table}
		echo "viga.9"
		"""

rule annotate_BLAST:
	input:
		faa=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_viga_ORFs.tot.faa",
		blast=(config['blast_db']),
	output:
		blast_output=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + "_blast_viralRefSeq.{sampling}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/{sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 32
	shell:
		"""
		blastp -num_threads {threads} -db {input.blast} -query {input.faa} \
		-outfmt "6 qseqid sseqid stitle qstart qend qlen slen qcovs evalue length"  > {output.blast_output}
		"""

rule cluster_proteins_viga:
	input:
		faa=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_viga_ORFs.tot.faa",
	output:
		mmseqs_out=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + "_viga_cluster.tsv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	# benchmark:
	# 	dirs_dict["BENCHMARKS"] +"/annotate_BLAST/{sampling}.tsv"
	message:
		"Clustering proteins with mmseqs"
	params:
		rep_name=REPRESENTATIVE_CONTIGS_BASE + "_viga",
		dir_mmseqs=dirs_dict["ANNOTATION"],
	threads: 16
	shell:
		"""
		cd {params.dir_mmseqs}
		mmseqs easy-cluster {input.faa} {params.rep_name} tmp --threads {threads} 
		"""

rule cluster_proteins:
	input:
		faa=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs.tot.faa",
	output:
		mmseqs_out=(dirs_dict["vOUT_DIR"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + "_cluster.tsv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	# benchmark:
	# 	dirs_dict["BENCHMARKS"] +"/annotate_BLAST/{sampling}.tsv"
	message:
		"Clustering proteins with mmseqs"
	params:
		rep_name=REPRESENTATIVE_CONTIGS_BASE,
		dir_mmseqs=dirs_dict["vOUT_DIR"],
	threads: 16
	shell:
		"""
		cd {params.dir_mmseqs}
		mmseqs easy-cluster {input.faa} {params.rep_name} tmp --threads {threads}
		"""

rule phynteny_annotation:
	input:
		pharokka_output="{contigs}_pharokka",
		phynteny_db=config["phynteny_db"]
	output:
		phynteny_output=directory("{contigs}_phynteny")
	params:
		pharokka_gbk="{contigs}_pharokka/pharokka.gbk"
	conda:
		dirs_dict["ENVS_DIR"] + "/phynteny.yaml"
	message:
		"Annotating Pharokka GenBank with Phynteny"
	threads: 1
	shell:
		"""
		rm -rf {output.phynteny_output}
		if [ -d {input.phynteny_db}/models ] && [ -s {input.phynteny_db}/confidence_dict.pkl ]; then
			phynteny {params.pharokka_gbk} -o {output.phynteny_output} -m {input.phynteny_db}/models -t {input.phynteny_db}/confidence_dict.pkl
		elif [ -s {input.phynteny_db}/confidence_dict.pkl ]; then
			phynteny {params.pharokka_gbk} -o {output.phynteny_output} -m {input.phynteny_db} -t {input.phynteny_db}/confidence_dict.pkl
		else
			phynteny {params.pharokka_gbk} -o {output.phynteny_output}
		fi
		if [ "$(find {output.phynteny_output} -type f \( -name '*.gbk' -o -name '*.gb' \) | wc -l)" -eq 0 ]; then
			cp {params.pharokka_gbk} {output.phynteny_output}/pharokka.gbk
		fi
		"""

rule clinker_figure:
	input:
		phynteny_output=("{contigs}_phynteny"),
	output:
		clinker=("{contigs}_clinker.html"),
	conda:
		dirs_dict["ENVS_DIR"] + "/clinker.yaml"
	message:
		"Creating genome visualization with clinker"
	params:
		clinker_dir="{contigs}_clinker",
		gb="{contigs}_genbank*.gbk"
	threads: 16
	shell:
		"""
		rm -rf {params.clinker_dir} {wildcards.contigs}_genbank*
		mkdir -p {params.clinker_dir}
		cd {params.clinker_dir}
		find {input.phynteny_output} -type f \( -name "*.gbk" -o -name "*.gb" \) -exec cp {{}} . \;
		gbk_files=$(find . -maxdepth 1 -type f \( -name "*.gbk" -o -name "*.gb" \) | sort)
		if [ -z "$gbk_files" ]; then
			echo "No GenBank files found in {input.phynteny_output}" >&2
			exit 1
		fi

		awk '{{f="tmp_record_" NR; print $0 "//" > f}}' RS='//' $gbk_files
		rm -f $gbk_files
		find . -type f -size -10c -delete

		for f in tmp_record_*; do
			locus=$(awk '/^LOCUS/ {{print $2; exit}}' "$f")
			if [ -n "$locus" ]; then
				mv "$f" "${{locus}}.gbk"
			else
				echo "Warning: could not extract LOCUS from $f" >&2
				rm "$f"
			fi
		done
		clinker *.gbk -p {output.clinker} -j {threads}

		rm -rf {params.clinker_dir}
		"""

rule lovis4u_figure:
	input:
		phynteny_output=("{contigs}_phynteny"),
	output:
		lovis4u=directory("{contigs}_lovis4u"),
	conda:
		dirs_dict["ENVS_DIR"] + "/lovis4u.yaml"
	message:
		"Creating genome visualization with LoVis4u"
	params:
		output_dir=lambda wildcards, output: os.path.abspath(str(output.lovis4u)),
		work_dir=lambda wildcards, output: os.path.abspath(str(output.lovis4u) + "_work")
	threads: 16
	shell:
		"""
		out_dir="{params.output_dir}"
		work_dir="{params.work_dir}"
		gb_dir="$work_dir/genbank"

		rm -rf "$out_dir" "$work_dir"
		mkdir -p "$gb_dir"
		find {input.phynteny_output} -type f \( -name "*.gbk" -o -name "*.gb" \) -exec cp {{}} "$gb_dir"/ \;
		if [ "$(find "$gb_dir" -maxdepth 1 -type f \( -name "*.gbk" -o -name "*.gb" \) | wc -l)" -eq 0 ]; then
			echo "No GenBank files found in {input.phynteny_output}" >&2
			exit 1
		fi

		mkdir -p "$(dirname "$out_dir")"
		cd "$work_dir"
		lovis4u --data
		lovis4u --linux
		lovis4u -gb "$gb_dir" -hl -smp mmseqs -c A4L -o "$out_dir"
		rm -rf "$work_dir"
		"""

checkpoint filtered_vOTU_visualization_decision:
	input:
		fasta="{contigs}.fasta",
	output:
		summary="{contigs}_visualization.tsv",
	params:
		max_contigs=VISUALIZATION_MAX_CONTIGS,
		tool=VISUALIZATION_TOOL,
		lovis4u_output="{contigs}_lovis4u/lovis4u.pdf",
		clinker_output="{contigs}_clinker.html"
	message:
		"Checking whether filtered vOTUs are small enough to visualize"
	shell:
		"""
		seq_count=$(grep -c '^>' {input.fasta} || true)
		status="render"
		visualization_output="{params.lovis4u_output}"
		if [ "{params.tool}" = "clinker" ]; then
			visualization_output="{params.clinker_output}"
		fi
		if [ "$seq_count" -eq 0 ]; then
			status="skipped_empty_fasta"
			visualization_output=""
		elif [ "$seq_count" -ge {params.max_contigs} ]; then
			status="skipped_too_many_contigs"
			visualization_output=""
		fi
		printf "tool\tmax_contigs\tseq_count\tstatus\toutput\n" > {output.summary}
		printf "{params.tool}\t{params.max_contigs}\t%s\t%s\t%s\n" "$seq_count" "$status" "$visualization_output" >> {output.summary}
		"""

def filtered_vOTU_visualization_inputs(wildcards):
	ckpt = checkpoints.filtered_vOTU_visualization_decision.get(contigs=wildcards.contigs)
	summary = ckpt.output["summary"]
	inputs = [summary, wildcards.contigs + "_phynteny"]

	with open(summary) as handle:
		header = handle.readline().strip().split("\t")
		values = handle.readline().strip().split("\t")
	row = dict(zip(header, values))

	if row.get("status") == "render":
		if VISUALIZATION_TOOL == "clinker":
			inputs.append(wildcards.contigs + "_clinker.html")
		else:
			inputs.append(wildcards.contigs + "_lovis4u")
	return inputs

rule filtered_vOTU_visualization:
	input:
		filtered_vOTU_visualization_inputs
	output:
		done="{contigs}_visualization.done",
	message:
		"Finalizing filtered vOTU annotation and visualization"
	shell:
		"""
		touch {output.done}
		"""

rule blasToRefSeq:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/{sequence}.fasta",
		refseq_db=(config['RefSeqViral_db']),
	output:
		blast_output=(dirs_dict["ANNOTATION"] + "/blast_output_ViralRefSeq_{sequence}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/blasToRefSeq/{sequence}.tsv"
	message:
		"Blast contigs agaist RefSeq database"
	threads: 32
	shell:
		"""
		blastn -num_threads {threads} -db {input.refseq_db} -query {input.fasta} \
		-outfmt "6 qseqid sseqid salltitles qstart qend qlen slen qcovs evalue length pident" > {output.blast_output}
		"""

rule blastToIMGVR:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/{sequence}.fasta",
		img_vr_db=(config['IMGVR_db']),
	output:
		blast_output=(dirs_dict["ANNOTATION"] + "/blast_output_IMGVR_{sequence}.csv"),
	params:
			img_vr_db=(config['IMGVR_db'] + "IMGVR_all_nucleotides"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/blasToIMGVR/{sequence}.tsv"
	message:
		"Blast contigs agaist IMG/VR database"
	threads: 32
	shell:
		"""
		blastn -num_threads {threads} -db {params.img_vr_db} -query {input.fasta} \
		-outfmt "6 qseqid sseqid salltitles qstart qend qlen slen qcovs evalue length pident" > {output.blast_output}
		"""

rule create_dbs_mmseqs2:
	input:
		MMseqs2_dir=(config['mmseqs_dir']),
		representatives=dirs_dict["vOUT_DIR"] + "/merged_scaffolds.tot_95-80.fna",
		reference=REPRESENTATIVE_CONTIGS
	output:
		index_representatives=dirs_dict["MMSEQS"] + "/representatives.index",
		index_reference=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".index",
		idx_reference=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".idx",
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REPRESENTATIVE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
	message:
		"Creating databases for reference and assembly mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/create_dbs_mmseqs2/tot.tsv"
	threads: 4
	shell:
		"""
		mmseqs createdb {input.representatives} {params.representatives_name}
		mmseqs createdb {input.reference} {params.reference_name}
		mmseqs createindex {params.reference_name} tmp --search-type 3
		"""
		
rule search_contigs_mmseqs2:
	input:
		index_representatives=dirs_dict["MMSEQS"] + "/representatives.index",
		index_reference=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".index",
		idx_reference=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".idx",
	output:
		results_index=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_search_results.index",
		results_table=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_best_search_results.txt",
		temp_dir=temp(directory(dirs_dict["MMSEQS"] + "/tmp")),
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REPRESENTATIVE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
	message:
		"Comparing reference and assembly mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/search_contigs_mmseqs2/tot.tsv"
	threads: 16
	shell:
		"""
		mkdir {output.temp_dir}
		mmseqs search {params.representatives_name} {params.reference_name} {params.results_name} {output.temp_dir} \
		--start-sens 1 --sens-steps 3 -s 7 --search-type 3 --threads {threads}
		mmseqs convertalis {params.representatives_name} {params.reference_name} {params.results_name} {output.results_table}
		"""

rule mapReadstoContigsPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{contigs}.fasta",
	output:
		sam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired.{sampling}_to_{contigs}.sam",
		bam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired.{sampling}_to_{contigs}.bam",
		sorted_bam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{sampling}_to_{contigs}.bam",
		sorted_bam_paired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{sampling}_to_{contigs}.bam.bai",
		sam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{sampling}_to_{contigs}.sam",
		bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{sampling}_to_{contigs}.bam",
		sorted_bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{sampling}_to_{contigs}.bam",
		sorted_bam_unpaired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{sampling}_to_{contigs}.bam.bai",
	params:
		db_name=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_{sampling}_to_{contigs}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/mapReadstoContigsPE/{sample}_{sampling}_{contigs}.tsv"
	threads: 8
	shell:
		"""
		bowtie2-build -f {input.scaffolds} {params.db_name} --threads {threads}
		#paired
		bowtie2 -x {params.db_name} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam_paired} --threads {threads}
		samtools view -b -S {output.sam_paired} > {output.bam_paired}
		samtools sort {output.bam_paired} -o {output.sorted_bam_paired}
		samtools index {output.sorted_bam_paired}
		#unpaired
		bowtie2 -x {params.db_name} -U {input.unpaired} -S {output.sam_unpaired} --threads {threads}
		samtools view -b -S {output.sam_unpaired} > {output.bam_unpaired}
		samtools sort {output.bam_unpaired} -o {output.sorted_bam_unpaired}
		samtools index {output.sorted_bam_unpaired}
		"""

rule detectNucleotideModifications:
	input:
		fastq_file=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore.fastq",
		fast5_dir=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore_single_fast5",
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
	output:
		plus_wig=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + ".fraction_modified_reads.plus.wig"),
		minus_wig=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + ".fraction_modified_reads.minus.wig"),
	params:
		representative_basename=REPRESENTATIVE_CONTIGS_BASE,
	message:
		"Detecting nucleotide modifications with tombo"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_tombo.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/detectNucleotideModifications/tot.tsv"
	threads: 16
	shell:
		"""
		tombo preprocess annotate_raw_with_fastqs --fast5-basedir {input.fast5_dir} --fastq-filenames {input.fastq_file} --overwrite --processes {theads}
		tombo resquiggle {input.fast5_dir} {input.representatives} --processes {threads}
		tombo detect_modifications de_novo --fast5-basedirs {input.fast5_dir} --statistics-file-basename {params.representative_basename}.de_novo --processes {threads}
		tombo text_output browser_files --fast5-basedirs {input.fast5_dir} --statistics-filename {params.representative_basename}.de_novo.tombo.stats \
		--genome-fasta {input.representatives} --browser-file-basename {params.representative_basename} --file-types fraction
		"""
rule gbk_to_faa:
	input:
		genbank=dirs_dict["ANNOTATION"] + "/{contigs}.tot_VIGA_annotated.gbk",
	output:
		faa=dirs_dict["ANNOTATION"] + "/{contigs}_viga_ORFs.tot.faa",
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	message:
		"Formating genbank proteins as amino acid fasta"
	shell:
		"""
		python2 ./scripts/gbk_to_faa.py {input.genbank} {output.faa}
		sed -i 's/\s.*$//' {output.faa}
		"""


checkpoint split_multi_fasta:
	input:
		mmseqs_out=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + "_viga_cluster.tsv"),
	output:
		faa_dir=directory(dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_hhpred/"),
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	message:
		"Splitting amino acid fasta into protein cluster fasta"
	params:
		base=REPRESENTATIVE_CONTIGS_BASE
	shell:
		"""
		rm -rf {output.faa_dir} || true
		mkdir {output.faa_dir}
		cd {output.faa_dir}
		awk '{{data=$2; print data >> "cluster_"$1"_list.txt"}}' {input.mmseqs_out}
		"""

rule fasta_to_a2m:
	input:
		faa=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_viga_ORFs.tot.faa",
		cluster=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_hhpred/" + "cluster_{protein}_list.txt",
	output:
		faa=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_hhpred/" + "cluster_{protein}.faa",
		aln=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_hhpred/" + "cluster_{protein}.aln",
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	message:
		"Splitting amino acid fasta into protein cluster fasta"
	params:
		base=REPRESENTATIVE_CONTIGS_BASE
	shell:
		"""
		seqtk subseq {input.faa} {input.cluster} > {output.faa}
		muscle -in {output.faa} -out {output.aln}
		"""
def aggregate_input_annotation(wildcards):
	checkpoint_output = checkpoints.split_multi_fasta.get(**wildcards).output[0]
	return expand(dirs_dict["ANNOTATION"] + "/results_" + REPRESENTATIVE_CONTIGS_BASE  + "_hhpred/" + "cluster_{protein}.txt",
		protein=glob_wildcards(os.path.join(checkpoint_output, "cluster_{protein}_list.txt")).protein)


rule merge_annotations:
	input:		
		blast_output=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + "_blast_viralRefSeq.tot.csv"),
		hhr=aggregate_input_annotation,
		csv=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot_VIGA_annotated.csv",
		mmseqs_out=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + "_viga_cluster.tsv"),
	output:
		annotation_table=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot" + "_annotation_table_merged.csv",
		annotation_table_hhpred=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot" + "_annotation_table_hhpred.csv",
	params:
		hhpred_dir=dirs_dict["ANNOTATION"] + "/results_" + REPRESENTATIVE_CONTIGS_BASE + "_hhpred/"
	message:
		"Merging annotation results"
	threads: 1
	run:
		import pandas as pd
		import os
		import numpy as np
		import sys
		from collections import Counter
		###################
		name=REPRESENTATIVE_CONTIGS_BASE

		qseqid =[]
		sseqid =[]
		sname =[]
		qstart =[]
		qend =[]
		qlen =[]
		slen =[]
		qcovs =[]
		pident =[]
		evalue =[]
		length=[]

		with open(input.blast_output) as f:
		    content = f.readlines()

		for line in content:
		    qseqid.append((line.strip().split("\t")[0]))
		    sseqid.append((line.strip().split("\t")[1].split(":")[-1]))
		    sname.append((line.strip().split("\t")[2]))
		    qstart.append(int(line.strip().split("\t")[3]))
		    qend.append(int(line.strip().split("\t")[4]))
		    qlen.append(int(line.strip().split("\t")[5]))
		    slen.append(int(line.strip().split("\t")[6]))
		    qcovs.append(int(line.strip().split("\t")[7]))
		    evalue.append(float(line.strip().split("\t")[8]))
		    length.append(int(line.strip().split("\t")[9]))

		blast_df=pd.DataFrame()
		blast_df["qseqid"]=qseqid
		blast_df["sseqid"]=sseqid
		blast_df["sname"]=sname
		blast_df["qstart"]=qstart
		blast_df["qend"]=qend
		blast_df["qlen"]=qlen
		blast_df["slen"]=slen
		blast_df["qcovs"]=qcovs
		blast_df["evalue"]=evalue
		blast_df["length"]=length

		blast_df=blast_df[~blast_df["sname"].str.contains("hypothetical")]
		blast_df=blast_df[~blast_df["sname"].str.contains("Hypothetical")]
		blast_df=blast_df[~blast_df["sname"].str.contains('^putative \[')]
		blast_df=blast_df[~blast_df["sname"].str.contains('^phage protein \[')]
		blast_df=blast_df[~blast_df["sname"].str.contains('\|phage protein \[')]
		blast_df=blast_df[~blast_df["sname"].str.contains('^Phage protein \[')]
		blast_df=blast_df[~blast_df["sname"].str.contains('^unnamed protein product')]
		blast_df=blast_df[~blast_df["sname"].str.contains('^uncharacterized protein')]
		blast_df=blast_df[~blast_df["sname"].str.contains('^PHIKZ')]

		blast_gr_df=blast_df.groupby('qseqid').first().sort_values(by=['qcovs'])
		blast_gr_df=blast_gr_df[(blast_gr_df["qcovs"]>95)& (blast_gr_df["evalue"]<0.0000000005)]
		blast_gr_df["Method"]='BLAST'
		blast_gr_df=blast_gr_df.add_prefix('BLAST_')

		blast_gr_df["BLAST_Description_short"] = blast_gr_df["BLAST_sname"].str.split("|", 1).str[-1].str.split("[", 1).str[0]

		#####
		#VIGA
		#####

		tbl_file=input.csv

		viga_df = pd.read_csv(tbl_file, sep="\t")
		viga_df=viga_df[(viga_df["Description"]!="Hypothetical protein") &
		               (~viga_df["Description"].str.contains('^PHIKZ')) &
		               (~viga_df["Description"].str.contains('kDa protein'))]
		viga_df.set_index('Protein ID', inplace=True)
		viga_df["Method"]='VIGA'
		viga_df=viga_df.add_prefix('VIGA_')
		viga_df["VIGA_Description_short"] = viga_df["VIGA_Description"].str.split("[").str[0]
		viga_df

		protein=[]
		length=[]
		description=[]
		prob=[]
		evalue=[]
		score=[]
		cols=[]
		identities=[]
		similarity=[]
		hh_df=pd.DataFrame()
		print(params.hhpred_dir)

		for filename in os.listdir(params.hhpred_dir):
		     if filename.endswith("hhr"):
		        print(filename)
		        with open(params.hhpred_dir + filename,"r", encoding='latin-1') as fi:
		            for ln in fi:
		                if ln.startswith("Match_columns"):
		                    len_prot=ln.split()[1]
		                if ln.startswith(">"):
		                    descr=ln.strip()
		                if ln.startswith("Probab"):
		                    description.append(descr)
		#                    print(filename.rsplit(".",1)[0])
		                    protein.append(filename.rsplit(".",1)[0])
		                    length.append(int(len_prot))
		                    prob.append(float(ln.strip().split(" ")[0].split("=")[1]))
		                    evalue.append(float(ln.strip().split(" ")[2].split("=")[1]))
		                    score.append(float(ln.strip().split(" ")[4].split("=")[1]))
		                    cols.append(int(ln.strip().split(" ")[6].split("=")[1]))
		                    identities.append(float(ln.strip().split(" ")[8].split("=")[1][:-1]))
		                    similarity.append(float(ln.strip().split(" ")[10].split("=")[1]))
		hh_df["protein"]=protein
		hh_df["length"]=length
		hh_df["description"]=description
		hh_df["prob"]=prob
		hh_df["score"]=score
		hh_df["evalue"]=evalue
		hh_df["score"]=score
		hh_df["cols"]=cols
		hh_df["identities"]=identities
		hh_df["similarity"]=similarity
		hh_df["percentage"]=hh_df["cols"]/hh_df["length"]


		hh_df=hh_df[hh_df["prob"]>95]
		print(len(hh_df))
		# hh_df["protein"]=hh_df["protein"].str.split("^cluster_").str[-1]
		cluster_file=input.mmseqs_out
		cluster_df=pd.read_csv(cluster_file,sep="\t", names=["rep", "mem"])
		hh_df=hh_df.merge(cluster_df, left_on="protein", right_on="rep", how="left").drop(columns=["rep", "protein"])
		hh_df=hh_df.rename(columns = {'mem':'protein'})


		hh_df.to_csv(output.annotation_table_hhpred, sep="\t")
		#CONFIDENCE 1
		hh_df1=hh_df[(hh_df["description"].str.contains("phage")) |
		      (hh_df["description"].str.contains("Phage")) |
		      (hh_df["description"].str.contains("virus")) |
		      (hh_df["description"].str.contains("Virus"))]

		#hh_df1=hh_df1[(~hh_df1["description"].str.contains("DUF"))]


		hh_df1=hh_df1[~hh_df1["description"].str.contains('Putative')]
		hh_df1=hh_df1[~hh_df1["description"].str.contains('putative')]
		hh_df1=hh_df1[~hh_df1["description"].str.contains('hypothetical')]
		hh_df1=hh_df1[~hh_df1["description"].str.contains('Hypothetical')]

		hh_df1["confidence"] = "1 Viral protein"
		#CONFIDENCE 2

		# hh_df2=hh_df[~(hh_df["description"].str.contains("phage")) &
		#       (~hh_df["description"].str.contains("Phage")) &
		#       (~hh_df["description"].str.contains("virus")) &
		#       (~hh_df["description"].str.contains("Virus")) &
		#       (hh_df["description"].str.contains("DUF")) ]
		hh_df2=hh_df[(hh_df["description"].str.contains("DUF"))]
		hh_df2["confidence"] = "2 DUF"

		#CONFIDENCE 3
		hh_df3=hh_df
		hh_df3["confidence"] = "3 Other proteins"

		hh_phred_combined=pd.concat([hh_df1, hh_df2, hh_df3])
		hh_grouped_df=hh_phred_combined.groupby('protein').first().sort_values(by=['confidence'])
		hh_grouped_df["Method"]='HHPHRED'
		hh_grouped_df=hh_grouped_df.add_prefix('HHPHRED_')
		hh_grouped_df

		df1=viga_df.merge(blast_gr_df, how='outer',left_index=True, right_index=True)
		df2=df1.merge(hh_grouped_df, how='outer',left_index=True, right_index=True)
		df2["annotation"]=df2["BLAST_Description_short"]
		df2["Method"]=df2["BLAST_Method"]
		df2.annotation.fillna(df2.VIGA_Description_short, inplace=True)
		df2.Method.fillna(df2.VIGA_Method, inplace=True)
		df2.annotation.fillna(df2.HHPHRED_description, inplace=True)
		df2.Method.fillna(df2.HHPHRED_Method, inplace=True)

		df2['protein'] = df2.index
		df2[['phage', 'protein_number']] = df2["protein"].str.rsplit("_",1, expand=True)
		df2['protein_number'] = pd.to_numeric(df2["protein_number"])

		df2=df2.sort_values(by=['phage', 'protein_number'])

		df2.to_csv(output.annotation_table, sep="\t")
rule makeblastdb:
	input:
		aa=dirs_dict["vOUT_DIR"]+ "/{fasta_name}_ORFs.{sampling}.faa",
	output:
		orfs_blast_db=dirs_dict["vOUT_DIR"]+ "/{fasta_name}_ORFs.{sampling}.faa.pdb",
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/{fasta_name}_makeblast_{sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 8
	shell:
		"""
		makeblastdb -in {input.aa} -dbtype prot
		"""

rule diamond:
	input:
		aa=dirs_dict["vOUT_DIR"]+ "/{fasta_name}_ORFs.{sampling}.faa",
		# orfs_blast_db=dirs_dict["vOUT_DIR"]+ "/{fasta_name}_ORFs.{sampling}.faa.pdb",
	output:
		diamond_output=(dirs_dict["ANNOTATION"] + "/{fasta_name}_ORFs_diamond_all.{sampling}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_BLAST/{fasta_name}_diamond_{sampling}.tsv"
	message:
		"Performing protein diamond for ORFs"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env4.yaml",
	threads: 64
	shell:
		"""
		diamond makedb --in {input.aa} --db {input.aa}
		diamond blastp --threads {threads} --db {input.aa} --query {input.aa} --max-target-seqs 0 \
		  --outfmt 6 qseqid sseqid length pident \
		  --out {output.diamond_output}
		"""

awk_command=r"""
{
	# Record keys and store the data
	if (!(($1,$2) in data)) {
		keys1[$1] = 1;
		keys2[$2] = 1;
	}
	data[$1, $2] = $3;
}
END {
	n = asorti(keys1, sorted_keys1);
	m = asorti(keys2, sorted_keys2);

	# Print header row
	printf "\t";
	for (j = 1; j <= m; j++) {
		printf "%s", sorted_keys2[j];
		if (j < m) printf "\t";
	}
	printf "\n";

	# Print data rows
	for (i = 1; i <= n; i++) {
		printf "%s", sorted_keys1[i];
		for (j = 1; j <= m; j++) {
			key = sorted_keys1[i] SUBSEP sorted_keys2[j];
			printf "\t%s", ((key in data) ? data[key] : "");
		}
		printf "\n";
	}
}"""

rule parse_diamond:
	input:
		diamond=(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs_diamond_all.tot.csv"),
		cummulative_length=(dirs_dict["vOUT_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs_coding_lengths.tot.txt")
	output:
		parsed_diamond=temp(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_parsed_diamond.txt"),
		parsed_diamond_first=(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_parsed_diamond_first.txt"),
		similarity=temp(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_similarity.txt"),
		similarity_dup=temp(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_similarity_dup.txt"),
		similarity_dup2=temp(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_similarity_dup2.txt"),
		similarity_dup3=temp(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_similarity_dup3.txt"),
		distance=temp(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_distance.txt"),
		distance_short=temp(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_distance_short.txt"),
		distance_short_full=temp(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_distance_short_full.txt"),
		pivot=(dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_distance_matrix_AAI.txt"),
	benchmark:
		dirs_dict["BENCHMARKS"] +"/BLAST_viridic/diamond_parsing.tsv"
	message:
		"Parsing blast results to AAI distance matrix"
	threads: 1
	shell:
		"""
		time awk 'BEGIN{{OFS="\t"}} {{split($1, a, "_"); split($2, b, "_"); $5 = substr($1, 1, length($1) - length(a[length(a)]) - 1); $6 = substr($2, 1, length($2) - length(b[length(b)]) - 1); $17 = ($3 * $4) / 100; print}}' {input.diamond} > {output.parsed_diamond}
		time awk -F'\t' '!seen[$1,$5,$6]++ {{print}}' {output.parsed_diamond} > {output.parsed_diamond_first}
		time awk 'BEGIN{{OFS="\t"}} {{key=$5 "\t" $6; sum[key]+=$7; count[key]++}} END{{for (key in sum) print key, sum[key], count[key]}}' {output.parsed_diamond_first} > {output.similarity}
		time awk 'NR==FNR{{a[$1,$2]=$3 FS $4; next}} {{if(($2,$1) in a) print $0, a[$2,$1]}}' {output.similarity} {output.similarity} > {output.similarity_dup}
		time awk 'NR==FNR{{a[$1]=$2; next}} {{if($1 in a) print $0, a[$1]}}' {input.cummulative_length} {output.similarity_dup} > {output.similarity_dup2}
		time awk 'NR==FNR{{a[$1]=$2; next}} {{if($2 in a) print $0, a[$2]}}' {input.cummulative_length} {output.similarity_dup2} > {output.similarity_dup3}
		time awk '{{ col9 = ( ($3 + $5) * 100 ) / ($7 + $8); col10 = 100 - col9; print $0, col9, col10 }}' {output.similarity_dup3}  | sed 's/\t/ /g' > {output.distance}
		time cut -d' ' -f1,2,10 {output.distance} > {output.distance_short}
		time awk 'BEGIN {{OFS=" "}} {{print}} {{matrix[$1][$2]=$3; contigs[$1]; contigs[$2]}} END {{for (i in contigs) {{for (j in contigs) {{if (!(i in matrix) || !(j in matrix[i])) {{print i, j, 100}}}}}}}}' {output.distance_short} > {output.distance_short_full}
		time awk {awk_command:q} {output.distance_short_full} > {output.pivot}
		"""

rule parse_diamond_isolates:
	input:
		diamond=(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_ORFs_diamond_all.tot.csv"),
		cummulative_length=(dirs_dict["vOUT_DIR"] + "/combined_positive_viral_contigs_ORFs_coding_lengths.tot.txt")
	output:
		parsed_diamond=temp(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_parsed_diamond.txt"),
		parsed_diamond_first=(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_parsed_diamond_first.txt"),
		similarity=temp(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_similarity.txt"),
		similarity_dup=temp(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_similarity_dup.txt"),
		similarity_dup2=temp(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_similarity_dup2.txt"),
		similarity_dup3=temp(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_similarity_dup3.txt"),
		distance=temp(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_distance.txt"),
		distance_short=temp(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_distance_short.txt"),
		distance_short_full=temp(dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_distance_short_full.txt"),
		pivot=dirs_dict["ANNOTATION"] + "/combined_positive_viral_contigs_distance_matrix_AAI.txt",
	benchmark:
		dirs_dict["BENCHMARKS"] +"/BLAST_viridic/diamond_parsing_combined.tsv"
	message:
		"Parsing blast results to AAI distance matrix"
	threads: 1
	shell:
		"""
		time awk 'BEGIN{{OFS="\t"}} {{split($1, a, "_"); split($2, b, "_"); $5 = substr($1, 1, length($1) - length(a[length(a)]) - 1); $6 = substr($2, 1, length($2) - length(b[length(b)]) - 1); $17 = ($3 * $4) / 100; print}}' {input.diamond} > {output.parsed_diamond}
		time awk -F'\t' '!seen[$1,$5,$6]++ {{print}}' {output.parsed_diamond} > {output.parsed_diamond_first}
		time awk 'BEGIN{{OFS="\t"}} {{key=$5 "\t" $6; sum[key]+=$7; count[key]++}} END{{for (key in sum) print key, sum[key], count[key]}}' {output.parsed_diamond_first} > {output.similarity}
		time awk 'NR==FNR{{a[$1,$2]=$3 FS $4; next}} {{if(($2,$1) in a) print $0, a[$2,$1]}}' {output.similarity} {output.similarity} > {output.similarity_dup}
		time awk 'NR==FNR{{a[$1]=$2; next}} {{if($1 in a) print $0, a[$1]}}' {input.cummulative_length} {output.similarity_dup} > {output.similarity_dup2}
		time awk 'NR==FNR{{a[$1]=$2; next}} {{if($2 in a) print $0, a[$2]}}' {input.cummulative_length} {output.similarity_dup2} > {output.similarity_dup3}
		time awk '{{ col9 = ( ($3 + $5) * 100 ) / ($7 + $8); col10 = 100 - col9; print $0, col9, col10 }}' {output.similarity_dup3}  | sed 's/\t/ /g' > {output.distance}
		time cut -d' ' -f1,2,10 {output.distance} > {output.distance_short}
		time awk 'BEGIN {{OFS=" "}} {{print}} {{matrix[$1][$2]=$3; contigs[$1]; contigs[$2]}} END {{for (i in contigs) {{for (j in contigs) {{if (!(i in matrix) || !(j in matrix[i])) {{print i, j, 100}}}}}}}}' {output.distance_short} > {output.distance_short_full}
		time awk {awk_command:q} {output.distance_short_full} > {output.pivot}
		"""

rule get_composition:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/{sequence}.fasta",
	output:
		composition=dirs_dict["ANNOTATION"]+ "/nucleotide_content_{sequence}.tsv",
	message:
		"Getting vOTUs nucleotide composition"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	threads: 1
	shell:
		"""
		seqtk comp {input.fasta} > {output.composition}
		"""
