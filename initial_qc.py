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

def releiden(adata):
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch")
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50) #n_neighbors default is 15
    sc.tl.leiden(adata, flavor="igraph", n_iterations=10, key_added='leiden')
    sc.tl.umap(adata)
    return adata

sc.settings.figdir = "results/qc_figures"

#Load initial data
adata = sc.read_h5ad(args.input)
start_shape=adata.shape
max_counts = adata.X.max()

adata.var_names_make_unique()
if adata.obs['batch'].dtype.name != 'category':
    adata.obs['batch'] = adata.obs['batch'].apply(lambda x: f'Batch{x}' if adata.obs["batch"].dtype.name != 'category' else x).astype('category')

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

print('trying scrublet...')
sc.pp.scrublet(adata, batch_key='batch')
predicted_doublet_idx=list(adata[adata.obs.predicted_doublet==True].obs.index)
if len(predicted_doublet_idx) > 0:
    print('number of predicted doublets...:', len(predicted_doublet_idx))
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

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')

sc.pl.highly_variable_genes(adata, save='scatter_hv_genes.png')

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save='pca_plot.png')
sc.pl.pca(adata, color=["batch", "batch", "id", "id"], dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)], ncols=2, size=2, save='pc_variation_plots.png')

sc.external.pp.harmony_integrate(adata, key='batch') #adding to try to ensure that clusters look better

sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.leiden(adata, flavor="igraph", n_iterations=10)
sc.tl.umap(adata)

try:
    u.create_umaps(adata=adata, adata_name='Yang', 
                   colnames=['sample_id', 'genotype', 'age', 'batch', 'date_processed', "leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"], 
                   date='7_21')
    e='no error'
except Exception as err:
    print('error in create umaps')
    e=str(err)

#additional leiden cluster filters...
cluster_counts = adata.obs['leiden'].value_counts(normalize=True)
small_clusters_005 = cluster_counts[cluster_counts < 0.005].index.tolist() #clusters with less than 0.05% of total cells
small_clusters_01 = cluster_counts[cluster_counts < 0.01].index.tolist() #clusters with less than 0.01% of total cells

sc.pl.umap(adata, color='leiden', groups=small_clusters_005, title='0.5% Removal', save='leiden_smallclusters005.png')
sc.pl.umap(adata, color='leiden', groups=small_clusters_01, title='1% Removal', save='leiden_smallclusters01.png')
sc.pl.umap(adata, color='doublet_score', title="Scrublet: Doublet Scores",save='doubletscores.png')

pre_filter2=adata.shape
pre_removal_cells=adata.shape[0]
cells_insmallclusters=adata[adata.obs.leiden.isin(small_clusters_005)].shape[0]
cells_highdoubletpredict=adata[adata.obs.doublet_score<0.25].shape[0]
adata=adata[(adata.obs.doublet_score<0.25) & ~(adata.obs.leiden.isin(small_clusters_005))]
print('cells removed: ', pre_removal_cells-adata.shape[0])

adata=releiden(adata)
sc.tl.rank_genes_groups(adata, groupby='leiden',rankby_abs=False)
sc.pl.umap(adata, color='leiden', title="Leiden Clusters",save='leiden_filtered.png')
sc.pl.umap(adata, color='leiden', title="Leiden Clusters",legend_loc='on data', legend_fontsize=8,save='leiden_filtered_ondata.png')

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
    f.write(f"Adata shape after standard initial QC: {pre_filter2}")
    f.write(f"Number of cells removed due to leiden cluster with less than 0.05% of total cells: {cells_insmallclusters}")
    f.write(f"Number of cells removed due to doublet prediction score over 0.25: {cells_highdoubletpredict}")
    f.write(f"Final adata shape: {pre_filter2}")
    f.write(f"Adata Structure: {adata}")
    f.write(f"Remaining cells: {adata.n_obs}, genes: {adata.n_vars}\n")
    f.write(f"Mean total counts: {adata.obs.total_counts.mean():.2f}\n")
    f.write(f"Mean % mitochondrial counts: {adata.obs.pct_counts_mt.mean():.2f}\n")
    f.write(f"Max % mitochondrial counts: {adata.obs.pct_counts_mt.max():.2f}\n")
    f.write(f"Error in Creating Umaps: {e}\n")