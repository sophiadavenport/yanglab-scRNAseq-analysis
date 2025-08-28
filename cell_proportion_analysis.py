import scanpy as sc
from scanpro import scanpro
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import itertools
import os
import re

parser=argparse.ArgumentParser()
parser.add_argument("--save_directory", required=True)
parser.add_argument("--adata_path", required=True)
parser.add_argument("--celltype_column", required=True)
parser.add_argument("--celltype_name", required=True)
parser.add_argument("--condition_column", required=True)
parser.add_argument("--samples_column", required=True)
args=parser.parse_args()

save_directory=args.save_directory
os.makedirs(save_directory, exist_ok=True)

adata=sc.read_h5ad(args.adata_path)

def run_scanpro(adata, clusters, condition, cluster_name, samples):
    out=scanpro(adata,clusters_col=clusters,conds_col=condition,samples_col=samples, robust=True, n_reps=10)
    
    #Scanpro Scatterplot results:
    scatterplotpath=os.path.join(save_directory,"proportions_scatterplot.png")
    out.plot(save=scatterplotpath)
    boxplotpath=os.path.join(save_directory,"proportions_boxplot.png")
    out.plot(kind='boxplot', save=boxplotpath)

    #Print counts:
    counts=adata.obs.groupby(condition)[clusters].value_counts()
    print(f"{cluster_name} Counts:\n", counts)

    #Plot celltype counts:
    counts_unstack=counts.unstack(fill_value=0)
    counts_unstack.plot(kind='bar',stacked=True, figsize=(10, 6),colormap='tab20')
    plt.ylabel('Cell Count')
    plt.xlabel(condition.replace("_", " ").title())
    plt.title(f'Count of {cluster_name} per {condition.replace("_", " ").title()}')
    plt.legend(title=cluster_name, bbox_to_anchor=(1.05, 1), loc='best', fontsize=6)
    plt.tight_layout()
    barplotpath=os.path.join(save_directory,"celltypecounts_barplot.png")
    plt.savefig(barplotpath)
    plt.close()

    #Plot percentages:
    percentages=counts_unstack.div(counts_unstack.sum(axis=1), axis=0) * 100
    percentages.plot(kind='bar', stacked=True, figsize=(10, 6),colormap='tab20')
    plt.ylabel('Percentage of Cells')
    plt.xlabel(condition.replace("_", " ").title())
    plt.title(f'Proportion of Cells per {condition.replace("_", " ").title()}')
    plt.legend(title=cluster_name, bbox_to_anchor=(1.05, 1), loc='best', fontsize=6)
    plt.tight_layout()
    percentbarplotpath=os.path.join(save_directory,"proportioncelltypes_barplot.png")
    plt.savefig(percentbarplotpath)
    plt.close()

    return out

results=run_scanpro(adata=adata, clusters=args.celltype_column, condition=args.condition_column, cluster_name=args.celltype_name, samples=args.samples_column)


conditions=adata.obs[args.condition_column].unique()

pairwise_results={}
for cond_pair in itertools.combinations(conditions, 2):
    out=scanpro(adata=adata, clusters_col=args.celltype_column, samples_col=args.samples_column, conds_col=args.condition_column,conditions=list(cond_pair), robust=True, n_reps=10)

    pairwise_results[f"{cond_pair[0]}_vs_{cond_pair[1]}"]=out.results

excel_path=f"save_directory/{args.celltype_name}_{args.condition_column}_comparison_results.xlsx"
with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
    for comparison_name, df in pairwise_results.items():
        sheet_name=re.sub(r'[:\\/?*\[\]]', '', comparison_name)[:31]
        df.to_excel(writer, sheet_name=sheet_name, index=True)