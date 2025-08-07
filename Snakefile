use_conda = True
configfile: "config.yaml"
DATA=config["sample_data_folders"]
DATA_DIR=config["data_dir"]
METADATA=config["metadata_path"]
MOUSE_REF=config["mouse_brain_reference"]
FIGURES_PATH=config['figures_path']
CELLBENDER_PARAMS=config['cellbender_params']
SAMPLES = config['samples']

rule all:
    input:
        f"{DATA_DIR}/yang_scRNAseq_mapped.h5ad", 'results/qc_report.txt'

rule remove_ambient_rna:
    input:
        matrix_dir = f"{DATA_DIR}/sample_mtxs/{{sample}}_filtered_feature_bc_matrix"
    output:
        h5 = "results/cellbender/{sample}_cellbender.h5"
    params:
        total_droplets = lambda wildcards: CELLBENDER_PARAMS[wildcards.sample]["total_droplets"],
        expected_cells = lambda wildcards: CELLBENDER_PARAMS[wildcards.sample]["expected_cells"]
    threads: 8
    resources:
        mem_mb=100000
    conda:
        "envs/cellbender_env.yaml"
    shell:
        """
        cellbender remove-background --input {input.matrix_dir} --output {output.h5} --expected-cells {params.expected_cells} \
            --total-droplets-included {params.total_droplets} --cuda --fpr 0.01
        """

rule create_h5ad:
    input:
       data=expand("results/cellbender/{sample}_cellbender_filtered.h5", sample=SAMPLES), metadata=METADATA
    output:
        f"{DATA_DIR}/yang_scRNAseq_cellbender.h5ad"
    conda:
        "envs/scanpy_env.yaml"
    shell:
        """
        python create_countsmtx.py --data_dir "results/cellbender" --meta {input.metadata} --output {output}
        """

rule initial_qc:
    input:
        f"{DATA_DIR}/yang_scRNAseq_cellbender.h5ad"
    output:
        formatted_h5ad=f"{DATA_DIR}/yang_scRNAseq_formatted.h5ad", qc_report="results/qc_report.txt"
    conda:
        'envs/scanpy_env.yaml'
    shell:
        """
        python initial_qc.py --input {input} --output {output.formatted_h5ad} --report {output.qc_report}
        """

rule map_annotations:
    input:
        reference=MOUSE_REF, formatted_h5ad=f"{DATA_DIR}/yang_scRNAseq_formatted.h5ad"
    output:
        annotated=f"{DATA_DIR}/yang_scRNAseq_mapped.h5ad"
    conda:
        'envs/tacco.yaml'
    shell:
        """
        python cell_annotations.py --reference {input.reference} --formatted {input.formatted_h5ad} --output {output.annotated}
        """