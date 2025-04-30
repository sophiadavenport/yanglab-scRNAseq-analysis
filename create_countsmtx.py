import scanpy as sc
import os
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Path to directory with 10x sample folders")
parser.add_argument("--output", required=True, help="Path to output .h5ad file")
args = parser.parse_args()

data_dir = Path(args.input)
sample_dirs = [f for f in data_dir.glob("*_filtered_feature_bc_matrix") if f.is_dir()]

adatas = []

for sample_dir in sample_dirs:
    sample_id = sample_dir.name.replace("_filtered_feature_bc_matrix", "")
    adata = sc.read_10x_mtx(
        sample_dir,  #path to matrix folder
        var_names='gene_ids',
        make_unique=True
    )

    adata.obs['sample_id'] = sample_id
    adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names] #Add sample ID as prefix to barcodes
    
    adatas.append(adata)

adata_combined = adatas[0].concatenate(*adatas[1:], batch_key="sample_batch", batch_categories=[a.obs['sample_id'][0] for a in adatas])

adata_combined.write(args.output)
