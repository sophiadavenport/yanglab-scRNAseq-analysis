use_conda = True
configfile: "config.yaml"
DATA=config["sample_data_folders"]
DATA_DIR=config["data_dir"]
METADATA=config["metadata_path"]
MOUSE_REF=config["mouse_brain_reference"]
FIGURES_PATH=config['figures_path']

rule all:
    input:
        'results2/qc_report.txt', f"{FIGURES_PATH}/results2/annotation_vis/celltype_counts_batch.png", f"results2/subtype_counts.xlsx"

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
        qc_report="results2/qc_report.txt"
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
        annotated=f"{DATA_DIR}/yang_scRNAseq_formatted_annotated.h5ad"
    conda:
        'envs/tacco.yaml'
    shell:
        """
        python cell_annotations.py --reference {input.reference} --formatted {input.formatted_h5ad} --output {output.annotated}
        """
    
rule visualize_anno_qc:
    input:
        adata=f"{DATA_DIR}/yang_scRNAseq_formatted_annotated.h5ad"
    output:
        celltype_counts=f"{FIGURES_PATH}/results2/annotation_vis/celltype_counts_batch.png", counts_excel=f"results2/subtype_counts.xlsx"
    conda:
        'envs/scanpy_env.yaml'
    params:
        save_folder=f"{FIGURES_PATH}/results2/annotation_vis/"
    shell:
        """
        python initial_visualizations.py --adata {input.adata} --save_folder {params.save_folder} --celltype_counts {output.celltype_counts} --excel {output.counts_excel}
        """