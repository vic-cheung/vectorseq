#%%
import numpy as np
from pathlib import Path
from vectorseq.utils import create_dir
from vectorseq.pipeline.pipelines import (
    ExpressionPlots,
    ExpressionPlotsPipeline,
    CreateUMAP,
)
from vectorseq.marker_constants import CortexGenes
import scanpy as sc
import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3250"
brain_region = "v1"
run_dir = data_dir / experiment_id / brain_region
all_cells_output_dir = run_dir / "all_cells"
figure_save_dir = create_dir(data_dir / "figures")
figure_save_subdir = create_dir(f"{experiment_id}_{brain_region}_all_cells")

# Choose clustering scheme to use
n_neighbors = 15
leiden_resolution = 1.2

# %% [markdown]
# ## Expression Plots of Top 50 ranked genes for each all_cells cluster
#%%
ranked_genes_expression = ExpressionPlots(
    source_path=all_cells_output_dir / "cluster" / "adata.h5ad",
    output_dir=all_cells_output_dir,
    output_subdir=f"expression_plots_top50",
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    rank_genes=True,
    n_genes=50,
)
results = ranked_genes_expression()
adata = results["adata"]

# %%
gene_expression_pipe = ExpressionPlotsPipeline(
    output_dir=all_cells_output_dir,
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    brain_region=brain_region,
)
adata = gene_expression_pipe()
# %% create cell-type mapping based off of n_neighbors=15, leiden_resolution=1.2,
cortex_dict = {
    "microglia": ["7", "10", "17", "21"],
    "mural cells": ["11"],
    "endothelia": ["18", "20", "5", "12", "3", "1", "2"],
    "astrocytes": ["16"],
    "excitatory": ["13", "15"],
    "inhibitory": ["22"],
    "oligodendrocytes": ["19", "8", "9", "4", "0", "6", "14"],
}

invert_cortex_dict = {
    value: key for key, values in cortex_dict.items() for value in values
}
clustering_scheme_name = f"neighbors_{n_neighbors}_leiden_{leiden_resolution}"
rename_clusters = adata.obs[clustering_scheme_name].apply(
    lambda x: invert_cortex_dict[x]
)

adata.obs["labeled_clusters"] = rename_clusters

# %%
sc._settings.ScanpyConfig.figdir = figure_save_dir / figure_save_subdir

# %%
sc.pl.dotplot(
    adata,
    var_names=list(CortexGenes().V1_MARKERS_WITH_TG),
    groupby="labeled_clusters",
    use_raw=False,
    dendrogram=True,
    standard_scale="var",
    swap_axes=True,
    cmap="viridis",
    figsize=None,
    colorbar_title="Relative expression",
    var_group_rotation=90,
    log=False,
    save="labeled_clusters.pdf",
)
# %%
sc.pl.dotplot(
    adata,
    var_names=["DreO", "Ef1a-FLPo", "Ef1a-mCherry-IRES-Cre"]
    + list(CortexGenes().V1_MARKERS),
    groupby="labeled_clusters",
    use_raw=False,
    dendrogram=True,
    standard_scale="var",
    swap_axes=True,
    cmap="viridis",
    figsize=None,
    colorbar_title="Relative expression",
    var_group_rotation=90,
    log=False,
    save="labeled_clusters_final_figure.pdf",
)

# %%
# custom umap plot
random_seed = 0
sc.tl.umap(
    adata,
    min_dist=0.5,
    spread=0.8,
    n_components=2,
    neighbors_key=f"neighbors_{n_neighbors}",
    random_state=random_seed,
)
sc.pl.umap(
    adata,
    neighbors_key=f"neighbors_{n_neighbors}",
    color="labeled_clusters",
    save=f"_labeled_clusters.svg",
    show=False,
)
# %%
