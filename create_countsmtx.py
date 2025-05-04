import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--data_dir", required=True, help="Path to directory with 10x sample folders")
parser.add_argument("--meta", required=True, help="Path to metadata file (should be in excel format)")
parser.add_argument("--output", required=True, help="Path to output .h5ad file")
args = parser.parse_args()

data_dir = Path(args.data_dir)
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

adata_combined = adatas[0].concatenate(*adatas[1:], batch_key="sample_id", batch_categories=[a.obs['sample_id'][0] for a in adatas])

metadata=pd.read_excel(args.meta)
metadata = metadata.set_index('sample_id')
for col in metadata.columns:
    if np.issubdtype(metadata[col].dtype, np.datetime64):
        metadata[col] = metadata[col].astype(str)
adata_combined.obs = adata_combined.obs.join(metadata, on='sample_id')

adata_combined.write(args.output)
