#%%
from pathlib import Path
from vectorseq.pipeline.stages import ExpressionPlots
from vectorseq.pipeline import ExpressionPlotsPipeline
from vectorseq.utils import create_dir
from vectorseq.marker_constants import SuperiorColliculusGenes

data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3454"
brain_region = "sc"
run_dir = data_dir / experiment_id / brain_region
excitatory_output_dir = run_dir / "excitatory"
figure_save_dir = create_dir(data_dir / "figures")
figure_save_subdir = create_dir(f"{experiment_id}_{brain_region}_excitatory_subset")

# Choose clustering scheme to use
n_neighbors = 55
leiden_resolution = 0.55

# %% [markdown]
# ## Expression Plots of Top 50 ranked genes for each excitatory cluster
#%%
ranked_genes_expression = ExpressionPlots(
    source_path=excitatory_output_dir / "cluster" / "adata.h5ad",
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
    source_data_path=excitatory_output_dir / "cluster" / "adata.h5ad",
    output_dir=figure_save_dir / figure_save_subdir,
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    brain_region=brain_region,
)
adata = gene_expression_pipe()

# %%
final_figure = ExpressionPlots(
    source_path=excitatory_output_dir / "cluster" / "adata.h5ad",
    output_dir=figure_save_dir,
    output_subdir=Path(figure_save_subdir) / "expression_plots",
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    var_names={
        "Genes of Interest": list(
            SuperiorColliculusGenes().SCRG_EXCITATORY_MARKERS_WITH_TG
        )
    },
    save_name=f"neighbors_{n_neighbors}_leiden_{leiden_resolution}_final_figure",
    n_genes=50,
)
results = final_figure()
adata = results["adata"]

# %%
