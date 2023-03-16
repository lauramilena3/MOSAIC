rule subsampleReadsIllumina_PE_test_depth:
	input:
		paired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}_read_count.txt",),
		unpaired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}_read_count.txt",),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq.gz"),

	output:
		forward_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_unpaired_clean.{sampling}.fastq.gz"),
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/subsampleReadsIllumina_PE_test_depth/{sample}_{subsample}_{sampling}.tsv"
	threads: 1
	resources:
		mem_mb=4000
	shell:
		"""
		#paired
		paired=$( cat {input.paired_sizes} )
		p=$(echo "$paired"*{wildcards.subsample}/100 | bc)
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads=$p sampleseed=1
		#unpaired
		unpaired=$( cat {input.unpaired_sizes} )
		up=$(echo "$unpaired"*{wildcards.subsample}/100 | bc)
		reformat.sh in={input.unpaired} out={output.unpaired} reads=$up sampleseed=1
		"""

rule normalizeReads_test_depth:
	input:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_unpaired_clean.{sampling}.fastq.gz"),
	output:
		forward_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_unpaired_norm.{sampling}.fastq.gz"),
	message:
		"Normalizing reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/normalizeReads_test_depth/{sample}_{subsample}_{sampling}.tsv"
	params:
		min_depth=config['min_norm'],
		max_depth=config['max_norm']
	threads: 4
	resources:
		mem_mb=6000
	shell:
		"""
		#PE
		#paired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
		target={params.max_depth} mindepth={params.min_depth} t={threads}
		#unpaired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth}
		"""

rule metaspadesPE_test_depth:
	input:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_unpaired_norm.{sampling}.fastq.gz"),
	output:
		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_filtered_scaffolds.{sampling}.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_{sampling}/filtered_list.txt"),
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_{sampling}/scaffolds.fasta",
		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_{sampling}"),
	message:
		"Assembling PE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/metaspadesPE_test_depth/{sample}_{subsample}_{sampling}.tsv",
	threads: 16
	shell:
		"""
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		--meta -t {threads} --only-assembler --memory 350
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		"""

# rule spadesPE_test_depth:
# 	input:
# 		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_norm.{sampling}.fastq.gz"),
# 		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_norm.{sampling}.fastq.gz"),
# 	output:
# 		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_filtered_scaffolds.{sampling}.fasta"),
# 		filtered_list=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}/filtered_list.txt"),
# 	params:
# 		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}/scaffolds.fasta",
# 		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}"),
# 	message:
# 		"Assembling PE reads with metaSpades"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/spadesPE_test_depth/{sample}_{subsample}_{sampling}.tsv"
# 	threads: 4
# 	shell:
# 		"""
# 		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired} -o {params.assembly_dir} \
# 		-t {threads} --only-assembler --memory 350
# 		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
# 		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
# 		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
# 		"""

# rule IDBA_UD_test_depth:
# 	input:
# 		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_norm.{sampling}.fastq"),
# 		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_norm.{sampling}.fastq"),
# 	output:
# 		interleaved=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_interleaved_paired_norm.{sampling}.fastq"),
# 		#scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{iteration}_spades_filtered_scaffolds.{sampling}.fasta"),
# 		scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_idbaud_{sampling}/contig.fa",
# 		filtered_scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_idbaud_filtered_scaffolds.{sampling}.fasta"),
# 	params:
# 		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}/scaffolds.fasta",
# 		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_idbaud_{sampling}"),
# 	message:
# 		"Assembling PE reads with metaSpades"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/IDBA_UD_test_depth/{sample}_{subsample}_{sampling}.tsv"
# 	threads: 4
# 	shell:
# 		"""
# 		./tools/idba/idba/bin/fq2fa --merge {input.forward_paired} {input.forward_paired} {output.interleaved}
# 		./tools/idba/idba/bin/idba_ud -r {output.interleaved} --num_threads {threads} -o {params.assembly_dir}
# 		cp {output.scaffolds} {output.filtered_scaffolds}
# 		sed -i "s/>/>{wildcards.sample}_{wildcards.subsample}_/g" {output.filtered_scaffolds}
# 		"""

###ITERATIVE ASSEMBLY

# rule shortReadAsemblySpadesPE__test_depth:
# 	input:
# 		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_iteration_{iteration}_forward_paired_norm.tot.fastq"),
# 		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_iteration_{iteration}}_reverse_paired_norm.tot.fastq"),
# 	output:
# 		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_filtered_scaffolds.{sampling}.fasta"),
# 		filtered_list=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}/filtered_list.txt"),
# 	params:
# 		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}/scaffolds.fasta",
# 		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}"),
# 	message:
# 		"Assembling PE reads with metaSpades"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	threads: 4
# 	shell:
# 		"""
# 		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired} -o {params.assembly_dir} \
# 		--meta -t {threads} --only-assembler --memory 350
# 		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
# 		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
# 		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
# 		"""

rule virSorter_test_depth:
	input:
		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_filtered_scaffolds.{sampling}.fasta"),
		virSorter_db=config['virSorter_db'],
	output:
		positive_fasta=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/final-viral-combined.fa",
		table_virsorter=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/final-viral-score.tsv",
		viral_boundary=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/final-viral-boundary.tsv",
		positive_list=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/{sample}_{subsample}_positive_VS_list_{sampling}.txt",
	params:
		out_folder=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}",
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/virSorter_test_depth/{sample}_{subsample}_{sampling}.tsv"
	threads: 4
	shell:
		"""
		virsorter run -w {params.out_folder} -i {input.scaffolds} -j {threads} --db-dir {input.virSorter_db}
		grep ">" {output.positive_fasta} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list}
		"""

rule extractViralContigs_test_depth:
	input:
		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_filtered_scaffolds.{sampling}.fasta"),
		positive_list=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/{sample}_{subsample}_positive_VS_list_{sampling}.txt",
	output:
		positive_list_derr=temp(dirs_dict["VIRAL_DIR"] + "/{sample}_{subsample}_virSorter_derr_{sampling}/positive_VS_list_{sampling}.txt"),
		final_viral_contigs=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_positive_virsorter.{sampling}.fasta",
	message:
		"Selecting Viral Contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/extractViralContigs_test_depth/{sample}_{subsample}_{sampling}.tsv"
	threads: 1
	shell:
		"""
		cat {input.positive_list} | sort | uniq > {output.positive_list_derr}
		seqtk subseq {input.scaffolds} {output.positive_list_derr} > {output.final_viral_contigs}
		sed "s/>/>{wildcards.sample}_{wildcards.subsample}_/g" -i {output.final_viral_contigs}

		"""

rule estimateGenomeCompletness_test_depth:
	input:
		final_viral_contigs=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_positive_virsorter.{sampling}.fasta",
	output:
		quality_summary=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virsorter_checkV_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virsorter_checkV_{sampling}/completeness.tsv",
		contamination=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virsorter_checkV_{sampling}/contamination.tsv",
		tmp=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virsorter_checkV_{sampling}/tmp"),
	params:
		checkv_outdir=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virsorter_checkV_{sampling}",
#		checkv_db=dirs_dict["ASSEMBLY_TEST"] + "/95-80_merged_positive_virsorter_checkV_{sampling}",
	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/estimateGenomeCompletness_test_depth/{sample}_{subsample}_{sampling}.tsv"
	threads: 4
	shell:
		"""
		rm -rf {params.checkv_outdir} || true
		checkv contamination {input.final_viral_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv completeness {input.final_viral_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv complete_genomes {input.final_viral_contigs} {params.checkv_outdir}
		checkv quality_summary {input.final_viral_contigs} {params.checkv_outdir}
		"""

rule vOUTclustering_test_depth:
	input:
		final_viral_contigs=dirs_dict["ASSEMBLY_TEST"] + "/merged_positive_virsorter.{sampling}.fasta",
	output:
		clusters=dirs_dict["ASSEMBLY_TEST"] + "/merged_positive_virsorter.{sampling}_95-80.clstr",
		representatives_temp=temp(dirs_dict["ASSEMBLY_TEST"]+ "/merged_positive_virsorter.{sampling}_95-80.fna"),
		representatives=dirs_dict["ASSEMBLY_TEST"]+ "/95-80_merged_positive_virsorter.{sampling}.fasta",
		representative_lengths=dirs_dict["ASSEMBLY_TEST"] + "/95-80_merged_positive_virsorter_lengths.{sampling}.txt",
	message:
		"Creating vOUTs with stampede"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/vOUTclustering_test_depth/{sampling}.tsv"
	threads: 1
	shell:
		"""
		./scripts/stampede-Cluster_genomes.pl -f {input.final_viral_contigs} -c 80 -i 95
		cat {output.representatives_temp} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
		$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.representative_lengths}
		cp {output.representatives_temp} {output.representatives}
		"""

rule viralStatsILLUMINA_test_depth:
	input:
		quast_dir=(config["quast_dir"]),
		representatives=dirs_dict["ASSEMBLY_TEST"]+ "/95-80_merged_positive_virsorter.{sampling}.fasta",
	output:
		quast_report_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/viral_representatives_statistics_quast_{sampling}"),
		quast_txt=dirs_dict["ASSEMBLY_TEST"] + "/viral_representatives_quast_report.{sampling}.txt",
	message:
		"Creating viral stats with quast"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/viralStatsILLUMINA_test_depth/{sampling}.tsv"
	threads: 4
	shell:
		"""
		{input.quast_dir}/quast.py {input.representatives} -o {output.quast_report_dir} --threads {threads}
		cp {output.quast_report_dir}/report.txt {output.quast_txt}
		"""

rule assemblyStatsILLUMINA_test_depth:
	input:
		quast_dir=(config["quast_dir"]),
		scaffolds_assembly=expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_{assemblers}_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES, subsample=subsample_test, assemblers=["metaspades"]),
	output:
		quast_report_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/assembly_statistics_quast_{sampling}"),
		quast_txt=dirs_dict["ASSEMBLY_TEST"] + "/assembly_quast_report.{sampling}.txt",
	message:
		"Creating viral stats with quast"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/assemblyStatsILLUMINA_test_depth/{sampling}.tsv"
	threads: 4
	shell:
		"""
		{input.quast_dir}/quast.py {input.scaffolds_assembly} -o {output.quast_report_dir} --threads {threads}
		cp {output.quast_report_dir}/report.txt {output.quast_txt}
		"""


# ###TEST ASSEMBLY TYPES
# rule spades_no_ecc_no_norm:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
# 	output:
# 		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_filtered_scaffolds_no_ecc_no_norm.{sampling}.fasta"),
# 		filtered_list=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_no_ecc_no_norm_{sampling}/filtered_list.txt"),
# 		assembly_graph=dirs_dict["ASSEMBLY_TEST"] +"/{sample}_assembly_graph_spades_no_ecc_no_norm.{sampling}.fastg",
# 	params:
# 		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_no_ecc_no_norm_{sampling}/scaffolds.fasta",
# 		assembly_graph=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_no_ecc_no_norm_{sampling}/assembly_graph.fastg",
# 		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_no_ecc_no_norm_{sampling}"),
# 		metagenomic_flag=METAGENOME_FLAG
# 	message:
# 		"Assembling PE reads with metaSpades"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/shortReadAsemblySpadesPE/{sample}_{sampling}.tsv"
# 	threads: 8
# 	shell:
# 		"""
# 		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
# 		{params.metagenomic_flag} -t {threads} --memory 350
#
# 		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
# 		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
# 		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
# 		cp {params.assembly_graph} {output.assembly_graph}
# 		sed "s/>/>{wildcards.sample}_/g" -i {output.scaffolds}
# 		"""
#
# rule spades_no_ecc_yes_norm:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm_noecc.{sampling}.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm_noecc.{sampling}.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm_noecc.{sampling}.fastq",
# 	output:
# 		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_filtered_scaffolds_no_ecc_yes_norm.{sampling}.fasta"),
# 		filtered_list=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_no_ecc_yes_norm_{sampling}/filtered_list.txt"),
# 		assembly_graph=dirs_dict["ASSEMBLY_TEST"] +"/{sample}_assembly_graph_spades_no_ecc_yes_norm.{sampling}.fastg",
# 	params:
# 		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_no_ecc_yes_norm_{sampling}/scaffolds.fasta",
# 		assembly_graph=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_no_ecc_yes_norm_{sampling}/assembly_graph.fastg",
# 		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_no_ecc_yes_norm_{sampling}"),
# 		metagenomic_flag=METAGENOME_FLAG
# 	message:
# 		"Assembling PE reads with metaSpades"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/shortReadAsemblySpadesPE/{sample}_{sampling}.tsv"
# 	threads: 8
# 	shell:
# 		"""
# 		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
# 		{params.metagenomic_flag} -t {threads} --memory 350
# 		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
# 		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
# 		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
# 		cp {params.assembly_graph} {output.assembly_graph}
# 		sed "s/>/>{wildcards.sample}_/g" -i {output.scaffolds}
# 		"""
#
# rule spades_yes_ecc_no_norm:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_ecc.{sampling}.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_ecc.{sampling}.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_ecc.{sampling}.fastq",
# 	output:
# 		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_filtered_scaffolds_yes_ecc_no_norm.{sampling}.fasta"),
# 		filtered_list=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_yes_ecc_no_norm_{sampling}/filtered_list.txt"),
# 		assembly_graph=dirs_dict["ASSEMBLY_TEST"] +"/{sample}_assembly_graph_spades_yes_ecc_no_norm.{sampling}.fastg",
# 	params:
# 		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_yes_ecc_no_norm_{sampling}/scaffolds.fasta",
# 		assembly_graph=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_yes_ecc_no_norm_{sampling}/assembly_graph.fastg",
# 		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_yes_ecc_no_norm_{sampling}"),
# 		metagenomic_flag=METAGENOME_FLAG
# 	message:
# 		"Assembling PE reads with metaSpades"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/shortReadAsemblySpadesPE/{sample}_{sampling}.tsv"
# 	threads: 8
# 	shell:
# 		"""
# 		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
# 		{params.metagenomic_flag} -t {threads} --memory 350 --only-assembler
# 		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
# 		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
# 		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
# 		cp {params.assembly_graph} {output.assembly_graph}
# 		sed "s/>/>{wildcards.sample}_/g" -i {output.scaffolds}
# 		"""
#
# rule spades_yes_ecc_yes_norm:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq",
# 	output:
# 		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_filtered_scaffolds_yes_ecc_yes_norm.{sampling}.fasta"),
# 		filtered_list=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_yes_ecc_yes_norm_{sampling}/filtered_list.txt"),
# 		assembly_graph=dirs_dict["ASSEMBLY_TEST"] +"/{sample}_assembly_graph_spades_yes_ecc_yes_norm.{sampling}.fastg",
# 	params:
# 		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_yes_ecc_yes_norm_{sampling}/scaffolds.fasta",
# 		assembly_graph=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_yes_ecc_yes_norm_{sampling}/assembly_graph.fastg",
# 		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_spades_yes_ecc_yes_norm_{sampling}"),
# 		metagenomic_flag=METAGENOME_FLAG
# 	message:
# 		"Assembling PE reads with metaSpades"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/shortReadAsemblySpadesPE/{sample}_{sampling}.tsv"
# 	threads: 8
# 	shell:
# 		"""
# 		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
# 		{params.metagenomic_flag} -t {threads} --memory 350 --only-assembler
# 		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
# 		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
# 		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
# 		cp {params.assembly_graph} {output.assembly_graph}
# 		sed "s/>/>{wildcards.sample}_/g" -i {output.scaffolds}
# 		"""
