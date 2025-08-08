#must use TACCO_env
import pandas as pd
import scanpy as sc
import numpy as np
import anndata as an
import tacco as tc # type: ignore
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--reference", required=True, help="Full path to mouse brain reference h5ad file.")
parser.add_argument("--formatted", required=True, help="Full path to h5ad file with QC done.")
parser.add_argument("--output", required=True, help="Full path to where the new adata should be written")
args = parser.parse_args()

adata=sc.read_h5ad(args.formatted)

adata.layers['normalized_counts']=adata.X
adata.X=adata.layers['counts'] #resetting to raw counts
adata.X=adata.X.astype('float32')

reference_adata=sc.read_h5ad(args.reference) #automatically loaded as counts and float32
subclass_counts=reference_adata.obs.subclass.value_counts()
keep_subclass=pd.DataFrame(subclass_counts[subclass_counts > 100]).index.tolist()
reference_adata=reference_adata[reference_adata.obs.subclass.isin(keep_subclass), :]

tc.tl.annotate(adata=adata, reference=reference_adata, annotation_key='subclass', result_key='Celltype')

try:
    adata.obs['Celltype']
except:
    try:
        print('manually adding unified celltype to obs')
        unified_celltype_matrix = adata.obsm['Celltype']
        max_indices = np.argmax(unified_celltype_matrix, axis=1)
        column_names = unified_celltype_matrix.columns
        adata.obs['Celltype'] = column_names[max_indices]
    except:
        print('issue adding obs from obsm')

adata.write(args.output)
