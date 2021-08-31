# %%
from vectorseq.utils import create_dir, check_gene_abundance
import scanpy as sc
import matplotlib
import re
from pathlib import Path
import pandas as pd

# %%
data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3454"
brain_region = "sc"
run_dir = data_dir / experiment_id / brain_region / "all_cells"
figure_save_dir = create_dir(data_dir / "quality_control")

adata_filter = sc.read_h5ad(run_dir / "filter/adata.h5ad")
adata = sc.read_h5ad(run_dir / "cluster/adata.h5ad")

# %%
n_neighbors = 55
leiden_resolution = 0.55
nonneuron_cluster_ids = [21, 24]
cluster_naming_scheme = f"neighbors_{n_neighbors}_leiden_{leiden_resolution}"
# %% find nonneuron clusters
nn_mask = adata.obs[cluster_naming_scheme].isin(
    [str(item) for item in nonneuron_cluster_ids]
)
adata_nn = adata_filter[nn_mask, :]
# %% find neuron clusters
adata_neu = adata_filter[~nn_mask, :]

# %% create placeholder dataframe
df = pd.DataFrame(
    [],
    columns=[
        "experiment_id",
        "brain_region",
        "gene",
        "cell_type",
        "percent_count_goi_mean",
        "percent_count_goi_std",
        "sum_of_all_goi_counts",
        "num_cells_goi_positive",
        "total_num_cells",
        "avg_count_per_positive_cell",
    ],
)

GENES_LIST = ["AAVrg-CAG-tdTomato", "AAVrg-CAG-GFP", "HSV-Cre"]
for gene in GENES_LIST:
    df_nn = check_gene_abundance(adata=adata_nn, gene_of_interest=gene)
    df_neu = check_gene_abundance(adata=adata_neu, gene_of_interest=gene)
    df = df.append(
        pd.DataFrame.from_dict(
            {
                "experiment_id": experiment_id,
                "brain_region": brain_region,
                "gene": gene,
                "cell_type": "nonneuron",
                "percent_count_goi_mean": df_nn.percent_count_goi.mean(),
                "percent_count_goi_std": df_nn.percent_count_goi.std(),
                "sum_of_all_goi_counts": df_nn.goi_counts.sum(),
                "num_cells_goi_positive": df_nn.shape[0],
                "total_num_cells": adata_nn.obs.shape[0],
                "avg_count_per_positive_cell": df_nn.goi_counts.sum() / df_nn.shape[0],
            },
            orient="index",
        ).T
    ).append(
        pd.DataFrame.from_dict(
            {
                "experiment_id": experiment_id,
                "brain_region": brain_region,
                "gene": gene,
                "cell_type": "neuron",
                "percent_count_goi_mean": df_neu.percent_count_goi.mean(),
                "percent_count_goi_std": df_neu.percent_count_goi.std(),
                "sum_of_all_goi_counts": df_neu.goi_counts.sum(),
                "num_cells_goi_positive": df_neu.shape[0],
                "total_num_cells": adata_neu.obs.shape[0],
                "avg_count_per_positive_cell": df_neu.goi_counts.sum()
                / df_neu.shape[0],
            },
            orient="index",
        ).T
    )
df.reset_index(drop=True, inplace=True)
df.to_csv(
    figure_save_dir / f"{experiment_id}_{brain_region}_all_cells_ambient_RNA_check.csv"
)


# %%
