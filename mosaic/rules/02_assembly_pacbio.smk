rule hifiasmPacbio:
	input:
		pacbio=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_pacbio}_pacbio_clean.{sampling}.fastq.gz"
	output:
		gfa=dirs_dict["ASSEMBLY_DIR"] + "/{sample_pacbio}_hifiasm.{sampling}.bp.p_ctg.gfa",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample_pacbio}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta"
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/{sample_pacbio}_hifiasm.{sampling}"
	message:
		"Assembling PacBio HiFi reads with hifiasm"
	conda:
		dirs_dict["ENVS_DIR"] + "/pacbio.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/hifiasmPacbio/{sample_pacbio}_{sampling}.tsv"
	threads: 16
	shell:
		"""
		hifiasm -o {params.prefix} -t {threads} {input.pacbio}
		awk '/^S/{{print ">"$2"\\n"$3}}' {output.gfa} > {output.scaffolds}
		sed "s/>/>{wildcards.sample_pacbio}_/g" -i {output.scaffolds}
		"""

rule errorCorrectPolypolishPacbioPE:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_pacbio}_forward_paired_clean.{sampling}.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_pacbio}_reverse_paired_clean.{sampling}.fastq.gz",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample_pacbio}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta"
	output:
		sam1=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample_pacbio}_polypolish_1.{sampling}.sam"),
		sam2=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample_pacbio}_polypolish_2.{sampling}.sam"),
		scaffolds_polypolish=dirs_dict["ASSEMBLY_DIR"] + "/polypolish_{sample_pacbio}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta"
	message:
		"Correcting PacBio HiFi assembly with Illumina reads using Polypolish"
	conda:
		dirs_dict["ENVS_DIR"] + "/pacbio.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/errorCorrectPolypolishPacbioPE/{sample_pacbio}_{sampling}.tsv"
	threads: 8
	shell:
		"""
		bwa index {input.scaffolds}
		bwa mem -t {threads} -a {input.scaffolds} {input.forward_paired} > {output.sam1}
		bwa mem -t {threads} -a {input.scaffolds} {input.reverse_paired} > {output.sam2}
		polypolish {input.scaffolds} {output.sam1} {output.sam2} > {output.scaffolds_polypolish}
		sed "s/>/>polypolish_/g" -i {output.scaffolds_polypolish}
		"""
	
rule estimateFungalGenomeCompletnessBUSCO:
	input:
		assembly=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta"
	output:
		summary=dirs_dict["vOUT_DIR"] + "/{sample}_BUSCO_{sampling}/short_summary.specific."+ str(config['busco_lineage']) +".{sample}_BUSCO_{sampling}.txt"
	params:
		out="{sample}_BUSCO_{sampling}",
		out_path=dirs_dict["vOUT_DIR"],
		lineage=config['busco_lineage']
	message:
		"Estimating fungal genome completeness with BUSCO"
	conda:
		dirs_dict["ENVS_DIR"] + "/pacbio.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/estimateFungalGenomeCompletnessBUSCO/{sample}_{sampling}.tsv"
	threads: 8
	shell:
		"""
		busco -i {input.assembly} -o {params.out} --out_path {params.out_path} -l {params.lineage} -m genome -c {threads}
		"""
