import pandas as pd
import scanpy as sc
import argparse
import re
import seaborn as sns
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser()
parser.add_argument("--adata", required=True, help="Path to input .h5ad file")
parser.add_argument("--save_folder", required=True, help="Path to folder for saving the plots")
parser.add_argument("--celltype_counts", required=True, help="Name of celltype counts bar plot. (Only tracked plot for this rule)")
parser.add_argument("--excel", required=True, help="Path to excel file.")
args = parser.parse_args()

sc.settings.figdir = args.save_folder

adata=sc.read_h5ad(args.adata)

####UMAP of Celltypes:
sc.pl.umap(adata, color='Celltype', save='celltype.png')

####UMAP of Broad Celltype Classes:
def assign_broadtype(subclass):
    if isinstance(subclass, str) and (subclass.lower().endswith('glut')):
        return "Neuron - Glut"
    elif isinstance(subclass, str) and (subclass.lower().endswith('gaba')):
        return "Neuron - Gaba"
    elif isinstance(subclass, str) and (subclass.lower().endswith('dopa')):
        return "Neuron - Dopa"
    elif isinstance(subclass, str) and (subclass.lower().endswith('imn')):
        return "Neuron - IMN"
    elif isinstance(subclass, str) and ('astro' in subclass.lower()):
        return "Astrocyte"
    else:
        cleaned = re.sub(r'\d+', '', subclass)
        cleaned = re.sub(r'NN$', '', cleaned)
        return cleaned.strip()

adata.obs['celltype_class'] = adata.obs.apply(
    lambda row: assign_broadtype(row['Celltype']), axis=1
)

cell_types= adata.obs['celltype_class'].unique()
neurons=[ct for ct in cell_types if ct.startswith('Neuron')]
nonneurons=[ct for ct in cell_types if ct not in neurons]
n_palette = sns.color_palette("Greys", n_colors=len(neurons))
other_palette = sns.color_palette("husl", n_colors=len(nonneurons))
color_map = dict(zip(neurons, n_palette))
color_map.update(dict(zip(nonneurons, other_palette)))
sc.pl.umap(adata, color='celltype_class', palette=color_map, title='Broad Cell Type', save='celltype_broad.png')

####UMAP of For Each Individual Broad Celltype Class:
def plot_umap_highlight(adata, obs_column, subtype_column,highlight_color='blue', background_color='lightgray'):
    unique_labels = adata.obs[obs_column].unique()

    for label in unique_labels:
        temp_col = f'{obs_column}_highlight_temp'
        adata.obs[temp_col] = adata.obs[obs_column].apply(
            lambda x: str(label) if x == label else 'Other'
        )
        color_map = {
            'Other': background_color,
            str(label): highlight_color
        }
        cur_save_path=f"{obs_column.replace(' ', '_')}_{label.replace(' ', '_')}.png"
        sc.pl.umap(adata,color=temp_col, palette=color_map, title=f'{label.replace('_', ' ').title()}', save=cur_save_path)

        temp_col_2 = f'{obs_column}_subtype_temp'
        adata.obs[temp_col_2] = adata.obs.apply(
            lambda row: row[subtype_column] if row[obs_column] == label else 'Other', axis=1
        )
        subtypes = adata.obs[temp_col_2].unique().tolist()
        subtype_colors = sns.color_palette(n_colors=len(subtypes))
        palette = dict(zip(subtypes, subtype_colors))
        palette['Other'] = background_color

        save_path_2=f"{obs_column.replace(' ', '_')}_{label.replace(' ', '_')}_subtypes.png"
        sc.pl.umap(adata, color=temp_col_2, palette=palette, title=f'{label.replace('_', ' ').title()}: Subtypes', save=save_path_2)
        
        del adata.obs[temp_col] #removing temp column from adata
        del adata.obs[temp_col_2]


plot_umap_highlight(adata, 'celltype_class', 'Celltype')

####Writing Counts to Excel:
with pd.ExcelWriter(args.excel) as writer:
    celltype_counts = adata.obs['Celltype'].value_counts().reset_index()
    celltype_counts.columns = ['Celltype', 'count']
    celltype_counts.to_excel(writer, sheet_name='Celltype_Counts', index=False)

    celltype_class_counts = adata.obs['celltype_class'].value_counts().reset_index()
    celltype_class_counts.columns = ['celltype_class', 'count']
    celltype_class_counts.to_excel(writer, sheet_name='Celltype_Class_Counts', index=False)

####Bar Chart: Celltype Counts
def plot_celltype_counts(adata, groupby_col, count_col, save_path):
    counts = adata.obs.groupby(groupby_col, observed=True)[count_col].value_counts()
    counts_unstacked = counts.unstack(fill_value=0)

    counts_unstacked.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='tab20')
    plt.ylabel('Cell Count')
    plt.title(f'Count of {count_col.replace('_', ' ').title()} per {groupby_col.replace('_', ' ').title()}')
    plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(save_path)

plot_celltype_counts(adata, 'batch', 'celltype_class', args.celltype_counts)
plot_celltype_counts(adata, 'sample_id', 'celltype_class', os.path.join(args.save_folder, 'celltype_counts_sampleid.png'))
plot_celltype_counts(adata, 'genotype', 'celltype_class', os.path.join(args.save_folder, 'celltype_counts_genotype.png'))

#####Total Counts UMAPs by Celltype:
for celltype in list(adata.obs.celltype_class.unique()):
    cur_adata=adata[adata.obs.celltype_class==celltype, :]
    if cur_adata.shape[0]>10: #only plotting for more than 10 cells
        sc.pp.neighbors(cur_adata, use_rep='X_pca_harmony')
        sc.tl.leiden(cur_adata, flavor="igraph", n_iterations=10)
        sc.tl.umap(cur_adata)

        sc.pl.umap(cur_adata, color='Celltype', title=f"{celltype}: Subtypes", save=f"{celltype}_subtypes.png")
        sc.pl.umap(cur_adata, color='total_counts', title=f"{celltype}: Total Counts", save=f"{celltype}_total_counts.png", vmax=60000)

    else:
        print(celltype, "has less than 10 cells, not plotting indiv umaps")