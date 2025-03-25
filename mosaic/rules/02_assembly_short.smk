#ruleorder: shortReadAsemblySpadesPE > shortReadAsemblySpadesSE
def input_error_correction(wildcards):
	params_ecc=""
	if wildcards.sample=="ALL":
		params_ecc="--only-assembler"
	return params_ecc

def input_threads_assembler(wildcards):
	use_threads=12
	if wildcards.sample=="ALL":
		use_threads=64
	return use_threads

rule shortReadAsemblySpadesPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq.gz",
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta"),
		assembly_graph=dirs_dict["ASSEMBLY_DIR"] +"/{sample}_assembly_graph_spades.{sampling}.fastg",
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/scaffolds.fasta",
		assembly_graph=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/assembly_graph.fastg",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}"),
		metagenomic_flag=METAGENOME_FLAG,
		error_correction=input_error_correction,
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/filtered_list.txt"),
	message:
		"Assembling PE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/shortReadAsemblySpadesPE/{sample}_{sampling}.tsv"
	threads: input_threads_assembler
	resources:
		mem_gb=450
	priority: 1
	shell:
		"""
		rm -rf {params.assembly_dir}
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		{params.metagenomic_flag} -t {threads} --memory {resources.mem_gb} {params.error_correction}
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {params.filtered_list}
		seqtk subseq {params.raw_scaffolds} {params.filtered_list} > {output.scaffolds}
		cp {params.assembly_graph} {output.assembly_graph}
		sed "s/>/>{wildcards.sample}_/g" -i {output.scaffolds}
		rm -rf {params.assembly_dir}
		"""
# rule shortReadAsemblySpadesSE:
# 	input:
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq"
# 	output:
# 		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta"),
# 		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/filtered_list.txt")
# 	params:
# 		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/scaffolds.fasta",
# 		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}")
# 	message:
# 		"Assembling SE reads with metaSpades"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	threads: 16
# 	shell:
# 		"""
# 		spades.py -s {input.unpaired} -o {params.assembly_dir} \
# 		--sc -t {threads} --only-assembler
# 		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
# 		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
# 		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
# 		"""

# def input_Quast(wildcards):
# 	input_list=[]
# 	input_list.extend(expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta",sample=SAMPLES))
# 	if NANOPORE:
# 		input_list.extend(expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_"+ LONG_ASSEMBLER + "_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta", sample=NANOPORE_SAMPLES))
# 	if CROSS_ASSEMBLY:
# 		input_list.append(dirs_dict["VIRAL_DIR"]+ "/ALL_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta")
# 	return input_list

def input_Quast(wildcards):
	input_list=[]
	input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta",sample=SAMPLES))
	if NANOPORE:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{{sampling}}.fasta", sample=NANOPORE_SAMPLES))
	if CROSS_ASSEMBLY:
		input_list.append(dirs_dict["ASSEMBLY_DIR"] + "/ALL_spades_filtered_scaffolds.{sampling}.fasta")
	return input_list

rule assemblyStats:
	input:
		scaffolds=input_Quast,
		quast_dir=(config["quast_dir"]),
	output:
		quast_report_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}"),
		quast_txt=dirs_dict["ASSEMBLY_DIR"] + "/assembly_quast_report.{sampling}.txt",
	message:
		"Creating assembly stats with quast"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/assemblyStats/{sampling}.tsv"
	threads: 1
	shell:
		"""
		{input.quast_dir}/quast.py {input.scaffolds} -o {output.quast_report_dir}
		cp {output.quast_report_dir}/report.txt {output.quast_txt}
		"""


# rule mergeAssembliesSHORT:
# 	input:
# 		scaffolds_spades=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta",sample=SAMPLES)
# 	output:
# 		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
# 		merged_assembly_len=dirs_dict["VIRAL_DIR"] + "/merged_scaffolds_lengths.{sampling}.txt",
# 	message:
# 		"Merging assembled contigs"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/mergeAssembliesSHORT/{sampling}.tsv"
# 	threads: 1
# 	shell:
# 		"""
# 		cat {input.scaffolds_spades} > {output.merged_assembly}
# 		cat {output.merged_assembly} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
# 		$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.merged_assembly_len}
# 		"""