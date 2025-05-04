use_conda = True
configfile: "config.yaml"
DATA=config["sample_data_folders"]
DATA_DIR=config["data_dir"]
METADATA=config["metadata_path"]

rule all:
    input:
        f"{DATA_DIR}/yang_scRNAseq_formatted.h5ad", 'results/qc_report.txt'

rule create_h5ad:
    input:
       data=DATA, metadata=METADATA
    output:
        f"{DATA_DIR}/yang_scRNAseq.h5ad"
    conda:
        "envs/scanpy_env.yaml"
    shell:
        """
        python create_countsmtx.py --data_dir {input.data} --meta {input.metadata} --output {output}
        """

rule initial_qc:
    input:
        f"{DATA_DIR}/yang_scRNAseq.h5ad"
    output:
        formatted_h5ad=f"{DATA_DIR}/yang_scRNAseq_formatted.h5ad",
        qc_report="results/qc_report.txt"
    conda:
        'envs/scanpy_env.yaml'
    shell:
        """
        python initial_qc.py --input {input} --output {output.formatted_h5ad} --report {output.qc_report}
        """