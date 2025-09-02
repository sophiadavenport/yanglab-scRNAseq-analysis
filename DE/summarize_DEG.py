import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--input_files', nargs='+', required=True)
parser.add_argument('--report', required=True)
parser.add_argument('--figures_dir', required=True)
args=parser.parse_args()

os.makedirs(args.figures_dir, exist_ok=True)

def plot_logFC_heatmap(matrix, tag, output_path, cluster_rows=True, cluster_cols=True):
    num_genes=matrix.shape[0]
    if num_genes <= 5:
        fig_height=8
    elif num_genes <= 10:
        fig_height=6
    else:
        fig_height=min(0.25 * num_genes, 30)
    font_size=max(3, min(6, 10 - 0.02 * num_genes))
    clean_index=[str(g).replace("_", " ") for g in matrix.index]
    clean_columns=[str(c).replace("_", " ") for c in matrix.columns]
    matrix.index=clean_index
    matrix.columns=clean_columns

    g=sns.clustermap(matrix, cmap="bwr", center=0, row_cluster=cluster_rows, col_cluster=cluster_cols, cbar_kws={'shrink': 0.75}, linecolor='lightgray')
    
    g.ax_row_dendrogram.set_visible(True)
    g.ax_col_dendrogram.set_visible(True)
    g.ax_heatmap.set_xlabel('Cell Type', fontsize=10)
    g.ax_heatmap.set_ylabel('Gene', fontsize=10)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=65, fontsize=8)
    g.ax_heatmap.set_yticks(np.arange(matrix.shape[0]))
    g.ax_heatmap.set_yticklabels(matrix.index, fontsize=font_size)
    g.figure.set_size_inches(10, fig_height)
    g.figure.suptitle(f"{tag} Significant Genes", fontsize=14, y=1.05, ha='center')
    cbar_ax=g.cax
    cbar_ax.set_title('logFC', fontsize=10, pad=10, loc='center')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

comparison_celltype_logfc={}
comparison_sig_genes={}
summary_records=[]

for filepath in args.input_files:
    df=pd.read_csv(filepath)
    filename=os.path.basename(filepath)
    celltype=filename.replace("_edger_results.csv", "")

    if df.empty:
        summary_records.append({"celltype": celltype,"comparison": "NA","upregulated": 0,"downregulated": 0})
        continue

    if "comparison" in df.columns:
        groups=df.groupby("comparison")
    else:
        groups=[("NA", df)]

    for comparison, subdf in groups:
        sig_df=subdf[(subdf["FDR"] < 0.05) & (subdf["logFC"].abs() >= 1)]
        up=sig_df[sig_df["logFC"] >= 1].shape[0]
        down=sig_df[sig_df["logFC"] <= -1].shape[0]

        summary_records.append({"celltype": celltype,"comparison": comparison,"upregulated": up,"downregulated": down})

        if comparison not in comparison_celltype_logfc:
            comparison_celltype_logfc[comparison]={}
            comparison_sig_genes[comparison]=set()

        sig_genes=sig_df["gene"].tolist()
        comparison_sig_genes[comparison].update(sig_genes)
        comparison_celltype_logfc[comparison].setdefault(celltype, {})
        comparison_celltype_logfc[comparison][celltype].update(
            sig_df.set_index("gene")["logFC"].to_dict()
        )

#dataframe for CSV
summary_df=pd.DataFrame(summary_records)
summary_csv=os.path.join(args.figures_dir, "DEG_summary.csv")
summary_df.to_csv(summary_csv, index=False)

report_lines=["Summary Report", "=" * 60 + "\n"]

for comparison, comp_df in summary_df.groupby("comparison"):
    report_lines.append(f"\ncomparison: {comparison}")
    report_lines.append("-" * 40)
    for _, row in comp_df.iterrows():
        report_lines.append(f"Cell Type: {row['celltype']}")
        report_lines.append(
            f"  Total significant genes: {row['upregulated'] + row['downregulated']}"
        )
        report_lines.append(f"  Upregulated: {row['upregulated']}")
        report_lines.append(f"  Downregulated: {row['downregulated']}\n")

with open(args.report, "w") as f:
    f.write("\n".join(report_lines))

#one heatmap per comparison:
for comparison, celltype_logfc in comparison_celltype_logfc.items():
    sorted_genes=sorted(comparison_sig_genes[comparison])
    sorted_celltypes=sorted(celltype_logfc.keys())
    if not sorted_genes or not sorted_celltypes:
        continue

    heatmap_matrix=pd.DataFrame(index=sorted_genes, columns=sorted_celltypes)
    for celltype in sorted_celltypes:
        for gene in sorted_genes:
            logfc=celltype_logfc[celltype].get(gene, 0)
            heatmap_matrix.loc[gene, celltype]=logfc

    heatmap_matrix=heatmap_matrix.astype(float)

    heatmap_path=os.path.join(
        args.figures_dir, f"logFC_heatmap_{comparison}.png"
    )
    tag=f"comparison: {comparison}"
    plot_logFC_heatmap(heatmap_matrix, tag, heatmap_path)