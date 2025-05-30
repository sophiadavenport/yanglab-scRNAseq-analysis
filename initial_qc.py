import scanpy as sc
import pandas as pd
import numpy as np
import anndata as an
import argparse
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import umaps as u
import matplotlib
import matplotlib.colors as mcolors

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Path to input .h5ad file")
parser.add_argument("--output", required=True, help="Path to output .h5ad file")
parser.add_argument("--report", required=True, help="Path to output report .txt")
args = parser.parse_args()

sc.settings.figdir = "results/qc_figures"

#Load initial data
adata = sc.read_h5ad(args.input)
start_shape=adata.shape
max_counts = adata.X.max()

adata.var_names_make_unique()
if np.issubdtype(adata.obs['batch'].dtype, np.integer):
    adata.obs['batch'] = 'Batch' + adata.obs['batch'].astype(str) #ensuring batch is not type int

#ensure mitochondrial genes are properly captured
mt_genes = [gene for gene in adata.var_names if gene.startswith('mt-')]
if len(mt_genes)!=13:
    print('mitochondrial genes not found using mt-')
adata.var["mt"] = adata.var_names.str.startswith("mt-") #lowercase since mouse
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True, size=0.1, save='violin_qc_metrics.png')

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save='scatter_pctcountsmt.png')

sc.pp.filter_cells(adata, min_genes=200) #trying with 200 genes removed
cell_filter_shape=adata.shape
sc.pp.filter_cells(adata, max_genes=6000) #upper genes expressed limit (for round 2 of filtering lowered to 6000)
cell_filter_shape_upper=adata.shape
sc.pp.filter_genes(adata, min_cells=100) #upping the number of cells required for genes to be expressed in (100)
gene_filter_shape=adata.shape
adata = adata[(adata.obs.pct_counts_mt < 5)] #reducing mitochondrial genes
mt_filter_shape=adata.shape

sc.pp.scrublet(adata, batch_key='batch')
predicted_doublet_idx=list(adata[adata.obs.predicted_doublet==True].obs.index)
if len(predicted_doublet_idx) > 0:
    adata = adata[~adata.obs['predicted_doublet']==True] 
doublet_filter_shape=adata.shape

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save='scatter_pctcountsmt_postfiltering.png')

sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True, size=0.1, save='violin_qc_metrics_postfiltering.png')
for batch in list(adata.obs.batch.unique()): #violin plots for each batch individually
    cur_adata=adata[adata.obs.batch==batch, :]
    sc.pl.violin(cur_adata, ["n_genes_by_counts", "total_counts"], jitter=0.4, multi_panel=True, size=0.1, save=f"{batch}_violin_qc_metrics_postfiltering.png")

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

if 'batch' in adata.obs and not pd.api.types.is_categorical_dtype(adata.obs['batch']): #highly variable genes requires category type for batch
    adata.obs['batch'] = adata.obs['batch'].astype('category')
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')

adata.obs['batch'] = adata.obs['batch'].astype('category')

sc.pl.highly_variable_genes(adata, save='scatter_hv_genes.png')

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save='pca_plot.png')
sc.pl.pca(adata, color=["batch", "batch", "id", "id"], dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)], ncols=2, size=2, save='pc_variation_plots.png')

sc.external.pp.scanorama_integrate(adata, key='batch')
sc.external.pp.harmony_integrate(adata, key='batch') #adding to try to ensure that clusters look better

sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.leiden(adata, flavor="igraph", n_iterations=10)
sc.tl.umap(adata)

try:
    u.create_umaps(adata=adata, adata_name='Yang', 
                   colnames=['sample_id', 'genotype', 'age', 'batch', 'date_processed', "leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"], 
                   date='5_5')
    e='no error'
except Exception as err:
    print('error in create umaps')
    e=str(err)

adata.write(args.output)

#Creating QC doc...
with open(args.report, "w") as f:
    f.write(f"Initial Cells by Genes: {start_shape}\n")
    f.write(f"Max Raw Counts: {max_counts}\n")
    f.write(f"Shape after filtering out cells with less than 200 genes expressed: {cell_filter_shape} ({start_shape[0]-cell_filter_shape[0]} cells removed)\n")
    f.write(f"Shape after filtering out cells with more than 6000 genes expressed: {cell_filter_shape_upper} ({cell_filter_shape[0]-cell_filter_shape_upper[0]} cells removed)\n")
    f.write(f"Shape after filtering out genes expressed in less than 100 cells: {gene_filter_shape} ({start_shape[1]-gene_filter_shape[1]} genes removed)\n")
    f.write(f"Shape after filtering out cells with over 5% mitochondrial genes: {mt_filter_shape} ({gene_filter_shape[0]-mt_filter_shape[0]} cells removed due to mt content)\n")
    f.write(f"Number of predicted doublets: {len(predicted_doublet_idx)} (Post doublet filter: {doublet_filter_shape})")
    f.write(f"Final adata shape: {adata.shape}")
    f.write(f"Adata Structure: {adata}")
    f.write(f"Remaining cells: {adata.n_obs}, genes: {adata.n_vars}\n")
    f.write(f"Mean total counts: {adata.obs.total_counts.mean():.2f}\n")
    f.write(f"Mean % mitochondrial counts: {adata.obs.pct_counts_mt.mean():.2f}\n")
    f.write(f"Max % mitochondrial counts: {adata.obs.pct_counts_mt.max():.2f}\n")
    f.write(f"Error in Creating Umaps: {e}\n")