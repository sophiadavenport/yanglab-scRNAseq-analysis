import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse
from pathlib import Path
import h5py

parser = argparse.ArgumentParser()
parser.add_argument("--data_dir", required=True, help="Path to directory with 10x sample folders")
parser.add_argument("--meta", required=True, help="Path to metadata file (should be in excel format)")
parser.add_argument("--output", required=True, help="Path to output .h5ad file")
args = parser.parse_args()

data_dir = Path(args.data_dir)
sample_files = list(data_dir.glob("*_cellbender.h5"))

adatas = []

for f in sample_files:
    sample_id = f.name.replace("_cellbender_filtered.h5", "")
    adata = sc.read_10x_h5(f)
    adata.var_names_make_unique()
    adata.obs['sample_id'] = sample_id
    adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]

    adatas.append(adata)

adata_combined = adatas[0].concatenate(*adatas[1:], batch_key="sample_id", batch_categories=[a.obs['sample_id'][0] for a in adatas])

metadata=pd.read_excel(args.meta)
metadata = metadata.set_index('sample_id')
for col in metadata.columns:
    if np.issubdtype(metadata[col].dtype, np.datetime64):
        metadata[col] = metadata[col].astype(str)
adata_combined.obs = adata_combined.obs.join(metadata, on='sample_id')

adata_combined.write(args.output)
