# if VIRSORTER:
# 	rule virSorter2:
# 		input:
# 			scaffolds_spades=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
# 			virSorter_db=config['virSorter_db'],
# 		output:
# 			positive_fasta=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter_{sampling}/final-viral-combined.fa",
# 			table_virsorter=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter_{sampling}/final-viral-score.tsv",
# 			positive_list=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter_{sampling}/positive_VS_list_{sampling}.txt",
# 			viral_boundary=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter_{sampling}/final-viral-boundary.tsv",
# 			iter=directory(dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter_{sampling}/iter-0"),
# 		params:
# 			out_folder=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter_{sampling}"
# 		message:
# 			"Classifing contigs with VirSorter"
# 		conda:
# 			dirs_dict["ENVS_DIR"] + "/vir2.yaml"
# 		benchmark:
# 			dirs_dict["BENCHMARKS"] +"/virSorter2/{sample}_{sampling}_illumina.tsv"
# 		threads: 8
# 		shell:
# 			"""
# 			virsorter run -w {params.out_folder} -i {input.scaffolds_spades} -j {threads} --db-dir {input.virSorter_db}
# 			grep ">" {output.positive_fasta} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list} || true
# 			"""
# 	rule extractViralContigs_VS:
# 		input:
# 			scaffolds_spades=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
# 			positive_list=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter_{sampling}/positive_VS_list_{sampling}.txt",
# 		output:
# 			positive_list_derr=temp(dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter_derr_{sampling}/positive_VS_list_{sampling}.txt"),
# 			positive_contigs=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
# 		message:
# 			"Selecting Viral Contigs"
# 		conda:
# 			dirs_dict["ENVS_DIR"] + "/vir.yaml"
# 		benchmark:
# 			dirs_dict["BENCHMARKS"] +"/extractViralContigs/{sample}_{sampling}.tsv"
# 		threads: 1
# 		shell:
# 			"""
# 			cat {input.positive_list} | sort | uniq > {output.positive_list_derr}
# 			seqtk subseq {input.scaffolds_spades} {output.positive_list_derr} > {output.positive_contigs}
# 			"""

# 	#FIX FOR ONLY LONG READS
# 	rule virSorter2_nanopore:
# 		input:
# 			scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{sampling}.fasta"),
# 			virSorter_db=config['virSorter_db'],
# 		output:
# 			positive_fasta=dirs_dict["VIRAL_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_virSorter_{sampling}/final-viral-combined.fa",
# 			table_virsorter=dirs_dict["VIRAL_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_virSorter_{sampling}/final-viral-score.tsv",
# 			positive_list=dirs_dict["VIRAL_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_virSorter_{sampling}/positive_VS_list_{sampling}.txt",
# 			viral_boundary=dirs_dict["VIRAL_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_virSorter_{sampling}/final-viral-boundary.tsv",
# 			iter=directory(dirs_dict["VIRAL_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_virSorter_{sampling}/iter-0"),
# 		params:
# 			out_folder=dirs_dict["VIRAL_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_virSorter_{sampling}"
# 		message:
# 			"Classifing contigs with VirSorter"
# 		conda:
# 			dirs_dict["ENVS_DIR"] + "/vir2.yaml"
# 		benchmark:
# 			dirs_dict["BENCHMARKS"] +"/virSorter2/{sample}_{sampling}_nanopore.tsv"
# 		threads: 8
# 		shell:
# 			"""
# 			virsorter run -w {params.out_folder} -i {input.scaffolds} -j {threads} --db-dir {input.virSorter_db}
# 			grep ">" {output.positive_fasta} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list} || true
# 			"""

# 	rule extractViralContigs_nanopore:
# 		input:
# 			scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{sampling}.fasta"),
# 			positive_list=dirs_dict["VIRAL_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_virSorter_{sampling}/positive_VS_list_{sampling}.txt",
# 		output:
# 			positive_contigs=dirs_dict["VIRAL_DIR"]+ "/{sample}_"+ LONG_ASSEMBLER + "_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
# 		message:
# 			"Selecting Viral Contigs"
# 		conda:
# 			dirs_dict["ENVS_DIR"] + "/vir.yaml"
# 		benchmark:
# 			dirs_dict["BENCHMARKS"] +"/extractViralContigs/{sample}_{sampling}.tsv"
# 		threads: 1
# 		shell:
# 			"""
# 			seqtk subseq {input.scaffolds} {input.positive_list} > {output.positive_contigs}
# 			"""
# else:

rule satellite_finder:
	input:
		genomad_outdir=dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}/",
	output:
		satellite_finder_outdir=directory(dirs_dict["VIRAL_DIR"] + "/{sample}_{sampling}_satellite_finder_{model}/"),
		faa_temp=dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}/{sample}_spades_filtered_scaffolds.{sampling}_annotate/{sample}_spades_filtered_scaffolds.{sampling}_proteins.faa_fixed",
	params:
		faa=dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}/{sample}_spades_filtered_scaffolds.{sampling}_annotate/{sample}_spades_filtered_scaffolds.{sampling}_proteins.faa",
		model="{model}"
	message:
		"Identifying viral satellites with satellite_finder"
	conda:
		dirs_dict["ENVS_DIR"] + "/satellite_finder.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/satellite_finder/{sample}_{sampling}_{model}.tsv"
	threads: 8
	shell:
		"""
		awk '/^>/ {print $1} !/^>/ {print}' {params.faa} | sed -e '/^>/ {{s/_/-/g; s/-\([^_-]*\)$/_\1/}}' > {output.faa_temp}
		apptainer run -H ${{HOME}} docker://gempasteur/satellite_finder:0.9.1 --db-type gembase --models {params.model} --sequence-db {output.faa_temp} -w {threads} -o {output.satellite_finder_outdir} --mute
		"""

rule genomad_viral_id:
	input:
		scaffolds_spades=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
		genomad_db=(config['genomad_db']),
	output:
		genomad_outdir=directory(dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}/"),
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		viral_fasta=dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}/{sample}_spades_filtered_scaffolds.{sampling}_summary/{sample}_spades_filtered_scaffolds.{sampling}_virus.fna",
	message:
		"Identifying viral contigs with geNomad"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/geNomad_viralID/{sample}_{sampling}_illumina.tsv"
	threads: 8
	shell:
		"""
		genomad end-to-end --cleanup --splits 8 -t {threads} {input.scaffolds_spades} {output.genomad_outdir} {input.genomad_db} --relaxed
		cat {params.viral_fasta} | sed "s/|/_/g" > {output.positive_contigs}
		"""

rule genomad_viral_id_nanopore:
	input:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{sampling}.fasta"),
		genomad_db=(config['genomad_db']),
	output:
		genomad_outdir=directory(dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}_"+ LONG_ASSEMBLER + "/"),
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/{sample}_"+ LONG_ASSEMBLER + "_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		viral_fasta=dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}_"+ LONG_ASSEMBLER + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{sampling}_summary/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{sampling}_virus.fna",
	message:
		"Identifying viral contigs with geNomad"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/geNomad_viralID/{sample}_{sampling}_nanopore.tsv"
	threads: 8
	shell:
		"""
		genomad end-to-end --cleanup --splits 8 -t {threads} {input.scaffolds} {output.genomad_outdir} {input.genomad_db} --relaxed
		cat {params.viral_fasta} | sed "s/|/_/g" > {output.positive_contigs}
		"""

# VIRAL FILTERING vOTUS
rule virSorter2:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		virSorter_db=config['virSorter_db'],
	output:
		positive_fasta=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/final-viral-combined.fa",
		table_virsorter=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/final-viral-score.tsv",
		positive_list=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/positive_VS_list_{sampling}.txt",
		# DRAM_tab=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/for-dramv/viral-affi-contigs-for-dramv.tab",
		# DRAM_fasta=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/for-dramv/final-viral-combined-for-dramv.fa",
		iter=directory(dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/iter-0"),
	params:
		out_folder=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}"
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/VirSorter2_/" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}_illumina.tsv"
	threads: 64
	shell:
		"""
		virsorter run -w {params.out_folder} -i {input.representatives} -j {threads} --db-dir {input.virSorter_db} \
				--include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --seqname-suffix-off  --provirus-off --min-length 0
		grep ">" {output.positive_fasta} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list} || true
		"""


rule genomad_vOTUs:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		genomad_db=(config['genomad_db']),
	output:
		virus_summary=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_virus_summary.tsv",
		plasmid_summary=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_plasmid_summary.tsv",
		viral_fasta=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_virus.fna",																																																
		positive_contigs=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/formatted_viral_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		positive_contigs_conservative=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/formatted_viral_" + REPRESENTATIVE_CONTIGS_BASE + "_conservative.{sampling}.fasta",
	params:
		genomad_outdir=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/",
	message:
		"Identifying viral contigs with geNomad"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/geNomad_viralID_filtering/{sampling}_nanopore.tsv"
	threads: 32
	shell:
		"""
		genomad end-to-end --cleanup --splits 8 -t {threads} {input.representatives} {params.genomad_outdir} {input.genomad_db} --conservative  
		cat {output.viral_fasta} | sed "s/|/_/g" > {output.positive_contigs_conservative}
		genomad end-to-end --cleanup --splits 8 -t {threads} {input.representatives} {params.genomad_outdir} {input.genomad_db}  
		cat {output.viral_fasta} | sed "s/|/_/g" > {output.positive_contigs}
		# leave genomad folder as --relaxed to have all contigs in the summary file
		genomad end-to-end --cleanup --splits 8 -t {threads} {input.representatives} {params.genomad_outdir} {input.genomad_db}  --relaxed
		"""

rule annotate_VIBRANT:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		VIBRANT_dir=os.path.join(workflow.basedir, config['vibrant_dir']),
	output:
		vibrant_circular=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_circular.{sampling}.csv",
		vibrant_positive=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_list.{sampling}.csv",
		vibrant_quality=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_quality.{sampling}.csv",
		vibrant_summary=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_summary_results.{sampling}.csv",
	params:
		vibrant_outdir=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/annotate_VIBRANT/{sampling}.tsv"
	message:
		"Annotating viral contigs with VIBRANT"
	threads: 64
	shell:
		"""
		rm -rf {params.vibrant_outdir} || true
		mkdir {params.vibrant_outdir} ; cd {params.vibrant_outdir}
		{input.VIBRANT_dir}/VIBRANT_run.py -i {input.representatives} -t {threads} -virome
		cut -f1 {params.vibrant_outdir}/VIBRANT_*/VIBRANT_results*/*complete_circular*tsv > {output.vibrant_circular}
		cp {params.vibrant_outdir}/VIBRANT_*/VIBRANT_phages_*/*phages_combined.txt {output.vibrant_positive}
		cp {params.vibrant_outdir}/VIBRANT_*/VIBRANT_results*/VIBRANT_genome_quality*.tsv {output.vibrant_quality}
		cp {params.vibrant_outdir}/VIBRANT_*/VIBRANT_results*/VIBRANT_summary_results*.tsv {output.vibrant_summary}
		"""

# rule what_the_phage_vOTUs:
# 	input:
# 		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
# 		WTP_dir=directory(os.path.join(workflow.basedir, config['WTP_dir'])),
# 	output:
# 		out_folder=directory(dirs_dict["VIRAL_DIR"] + "/WhatThePhage_{sampling}"),
# 	# params:
# 	# 	out_folder=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}"
# 	message:
# 		"Classifing contigs with VirSorter"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/whatThePhage/{sampling}.tsv"
# 	threads: 32
# 	shell:
# 		"""
# 		nextflow run replikation/What_the_Phage -r v1.2.0 --cores {threads} --output {output.out_folder} -profile local,docker --fasta {input.representatives}
# 		"""
