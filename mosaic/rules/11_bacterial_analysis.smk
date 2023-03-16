rule estimateBacterialGenomeCompletness:
    input:
        corrected2_racon=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
        checkm_db=(config['checkm_db']),
    output:
        checkMoutdir_temp=temp(directory(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}_temp")),
        checkMoutdir=directory(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}"),
    params:
        checkv_db=dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{sampling}",
    log:
        checkMoutdir=(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}_log"),
    message:
        "Estimating genome completeness with CheckM "
    conda:
        dirs_dict["ENVS_DIR"] + "/env5.yaml"
    benchmark:
        dirs_dict["BENCHMARKS"] +"/estimateGenomeCompletness/{sample}_{sampling}_checkm.tsv"
    threads: 4
    shell:
        """
        mkdir {output.checkMoutdir_temp}
        cp {input.corrected2_racon} {output.checkMoutdir_temp}
        cd {output.checkMoutdir_temp}
        checkm lineage_wf -t {threads} -x fasta {output.checkMoutdir_temp} {output.checkMoutdir} 1> {log}
        """
