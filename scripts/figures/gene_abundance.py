#%%
import scanpy as sc
import pandas as pd
from pathlib import Path
from vectorseq.utils import check_gene_abundance, create_dir
from vectorseq.marker_constants import BrainGenes

data_dir = Path("/spare_volume/vectorseq-data")
figure_save_dir = create_dir(data_dir / "gene_abundance")

#%% [markdown]
# ## Gene Abundance Table for Experiment: 3250, Brain Region: v1
#%%
experiment_id = "3250"
brain_region = "v1"
run_dir = data_dir / experiment_id / brain_region
all_cells_output_dir = create_dir(run_dir / "all_cells")
adata = sc.read_h5ad(all_cells_output_dir / "filter" / "adata.h5ad")

filtered_tg_list = [
    gene for gene in BrainGenes.TG_MARKERS if gene.upper() in adata.obs.columns
]
endogenous_genes_list = [
    "Snap25",
    "Rbfox3",
    "Slc17a6",
    "Camk2a",
    "Gad1",
    "Gad2",
    "Mog",
    "Flt1",
]
gene_list = filtered_tg_list + endogenous_genes_list

count_fractions_df = pd.DataFrame()
for gene in gene_list:
    temp = check_gene_abundance(adata, gene_of_interest=gene)
    if not temp.empty:
        count_fractions_df = count_fractions_df.append(
            pd.DataFrame.from_dict(
                {
                    "gene": gene,
                    "number_of_expressing_cells": temp.shape[0],
                    "number_of_reads": temp.goi_counts.sum(),
                    "abundance_in_expressing_cells": f"{round(temp.percent_count_goi.mean(),2)} +/- {round(temp.percent_count_goi.std(),2)}",
                },
                orient="index",
            ).T
        )
        print(f"{gene} detected.")
    else:
        print(f"{gene} not detected.")

count_fractions_df.set_index(keys="gene", drop=True, inplace=True)
count_fractions_df.to_csv(
    figure_save_dir / f"{experiment_id}_{brain_region}_all_cells_gene_abundance.csv"
)

# %%
#%% [markdown]
# ## Gene Abundance Table for Experiment: 3382, Brain Region: snr
#%%
experiment_id = "3382"
brain_region = "snr"
run_dir = data_dir / experiment_id / brain_region
all_cells_output_dir = create_dir(run_dir / "all_cells")

adata = sc.read_h5ad(all_cells_output_dir / "filter" / "adata.h5ad")

filtered_tg_list = [
    gene for gene in BrainGenes.TG_MARKERS if gene.upper() in adata.obs.columns
]
endogenous_genes_list = [
    "Snap25",
    "Rbfox3",
    "Slc17a6",
    "Camk2a",
    "Gad1",
    "Gad2",
    "Mog",
    "Flt1",
]
gene_list = filtered_tg_list + endogenous_genes_list

count_fractions_df = pd.DataFrame()
for gene in gene_list:
    temp = check_gene_abundance(adata, gene_of_interest=gene)
    if not temp.empty:
        count_fractions_df = count_fractions_df.append(
            pd.DataFrame.from_dict(
                {
                    "gene": gene,
                    "number_of_expressing_cells": temp.shape[0],
                    "number_of_reads": temp.goi_counts.sum(),
                    "abundance_in_expressing_cells": f"{round(temp.percent_count_goi.mean(),2)} +/- {round(temp.percent_count_goi.std(),2)}",
                },
                orient="index",
            ).T
        )
        print(f"{gene} detected.")
    else:
        print(f"{gene} not detected.")

count_fractions_df.set_index(keys="gene", drop=True, inplace=True)
count_fractions_df.to_csv(
    figure_save_dir / f"{experiment_id}_{brain_region}_all_cells_gene_abundance.csv"
)
# %%
#%% [markdown]
# ## Gene Abundance Table for Experiment: 3454, Brain Region: sc
#%%
data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3454"
brain_region = "sc"
run_dir = data_dir / experiment_id / brain_region
all_cells_output_dir = create_dir(run_dir / "all_cells")

adata = sc.read_h5ad(all_cells_output_dir / "filter" / "adata.h5ad")

filtered_tg_list = [
    gene for gene in BrainGenes.TG_MARKERS if gene.upper() in adata.obs.columns
]
endogenous_genes_list = [
    "Snap25",
    "Rbfox3",
    "Slc17a6",
    "Camk2a",
    "Gad1",
    "Gad2",
    "Mog",
    "Flt1",
]
gene_list = filtered_tg_list + endogenous_genes_list

count_fractions_df = pd.DataFrame()
for gene in gene_list:
    temp = check_gene_abundance(adata, gene_of_interest=gene)
    if not temp.empty:
        count_fractions_df = count_fractions_df.append(
            pd.DataFrame.from_dict(
                {
                    "gene": gene,
                    "number_of_expressing_cells": temp.shape[0],
                    "number_of_reads": temp.goi_counts.sum(),
                    "abundance_in_expressing_cells": f"{round(temp.percent_count_goi.mean(),2)} +/- {round(temp.percent_count_goi.std(),2)}",
                },
                orient="index",
            ).T
        )
        print(f"{gene} detected.")
    else:
        print(f"{gene} not detected.")

count_fractions_df.set_index(keys="gene", drop=True, inplace=True)
count_fractions_df.to_csv(
    figure_save_dir / f"{experiment_id}_{brain_region}_all_cells_gene_abundance.csv"
)
#%%
