use_conda = True

rule all:
    input:
        '/lab/solexa_sun/shared_Data/scRNAseq/yang/yang_scRNAseq_formatted.h5ad'

rule create_h5ad:
    input:
        '/lab/solexa_sun/shared_Data/scRNAseq/yang/sample_mtxs'
    output:
        '/lab/solexa_sun/shared_Data/scRNAseq/yang/yang_scRNAseq.h5ad'
    conda:
        "envs/scanpy_env.yaml"
    shell:
        """
        python create_countsmtx.py \
            --input {input} \
            --output {output}
        """

rule initial_qc:
    input:
        '/lab/solexa_sun/shared_Data/scRNAseq/yang/yang_scRNAseq.h5ad'
    output:
        '/lab/solexa_sun/shared_Data/scRNAseq/yang/yang_scRNAseq_formatted.h5ad'
    conda:
        'envs/scanpy_env.yaml'
    shell:
        """
        python initial_qc.py
        """