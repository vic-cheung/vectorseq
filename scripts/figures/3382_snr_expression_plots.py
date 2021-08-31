#%%
from pathlib import Path
from vectorseq.pipeline.stages import ExpressionPlots
from vectorseq.pipeline import ExpressionPlotsPipeline
from vectorseq.utils import create_dir
from vectorseq.marker_constants import VentralMidbrainGenes

data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3382"
brain_region = "snr"
run_dir = data_dir / experiment_id / brain_region
all_cells_output_dir = create_dir(run_dir / "all_cells")
inhibitory_output_dir = run_dir / "inhibitory"
figure_save_dir = create_dir(data_dir / "figures")
figure_save_subdir = create_dir(f"{experiment_id}_{brain_region}_inhibitory_subset")

# Choose clustering scheme to use
n_neighbors = 45
leiden_resolution = 0.6

# %% [markdown]
# ## Expression Plots of Top 50 ranked genes for each inhibitory cluster

# %%
ranked_genes_expression = ExpressionPlots(
    source_path=inhibitory_output_dir / "cluster" / "adata.h5ad",
    output_dir=figure_save_dir,
    output_subdir=Path(figure_save_subdir) / "expression_plots_top50",
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    rank_genes=True,
    n_genes=50,
)
results = ranked_genes_expression()
adata = results["adata"]

# %%
gene_expression_pipe = ExpressionPlotsPipeline(
    source_data_path=inhibitory_output_dir / "cluster" / "adata.h5ad",
    output_dir=figure_save_dir / figure_save_subdir,
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    brain_region=brain_region,
)
adata = gene_expression_pipe()

# %%
final_figure = ExpressionPlots(
    source_path=inhibitory_output_dir / "cluster" / "adata.h5ad",
    output_dir=figure_save_dir,
    output_subdir=Path(figure_save_subdir) / "expression_plots",
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    var_names={
        "Genes of Interest": list(
            VentralMidbrainGenes().VENTRAL_MIDBRAIN_MARKERS_WITH_TG_SUBSET
        )
    },
    save_name=f"neighbors_{n_neighbors}_leiden_{leiden_resolution}_final_figure",
    n_genes=50,
)
results = final_figure()
adata = results["adata"]

# %%
