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
		positive_fasta=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_virSorter_{sampling}/final-viral-combined.fa",
		table_virsorter=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_virSorter_{sampling}/final-viral-score.tsv",
		positive_list=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_virSorter_{sampling}/positive_VS_list_{sampling}.txt",
		viral_boundary=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_virSorter_{sampling}/final-viral-boundary.tsv",
		iter=directory(dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_virSorter_{sampling}/iter-0"),
	params:
		out_folder=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_virSorter_{sampling}"
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/virSorter2/" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}_illumina.tsv"
	threads: 32
	shell:
		"""
		virsorter run -w {params.out_folder} -i {input.representatives} -j {threads} --db-dir {input.virSorter_db}
		grep ">" {output.positive_fasta} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list} || true
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
