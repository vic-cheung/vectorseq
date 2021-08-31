#%%
from pathlib import Path
from vectorseq.pipeline.stages import CreateUMAP
from vectorseq.utils import create_dir

data_dir = Path("/spare_volume/vectorseq-data")
figure_save_dir = create_dir(data_dir / "figures")

#%% [markdown]
# ## Overlay Genes of Interest on UMAP plots for Experiment: 3454, Brain Region: sc
#%%
experiment_id = "3454"
brain_region = "sc"
run_dir = data_dir / experiment_id / brain_region
all_cells_output_dir = create_dir(run_dir / "all_cells")

# Choose clustering scheme to use
n_neighbors = 55
leiden_resolution = 0.55

gene_list = [
    "Flt1",
    "Mog",
    "C1qa",
    "Aqp4",
    "Slc17a6",
    "Gad2",
    "AAVrg-CAG-GFP",
    "AAVrg-CAG-tdTomato",
    "HSV-Cre",
]

for gene in gene_list:
    create_umap = CreateUMAP(
        source_path=all_cells_output_dir / "cluster" / "adata.h5ad",
        output_dir=figure_save_dir / f"{experiment_id}_{brain_region}_all_cells",
        n_neighbors=n_neighbors,
        leiden_resolution=leiden_resolution,
        color=gene,
        save_adata=False,
    )
    results = create_umap()

#%% [markdown]
# ## Overlay Genes of Interest on UMAP plots for Experiment: 3382, Brain Region: snr
#%%
experiment_id = "3382"
brain_region = "snr"
run_dir = data_dir / experiment_id / brain_region
all_cells_output_dir = create_dir(run_dir / "all_cells")

# Choose clustering scheme to use
n_neighbors = 60
leiden_resolution = 0.6

gene_list = [
    "Flt1",
    "Mog",
    "C1qa",
    "Aqp4",
    "Slc17a6",
    "Gad1",
    "Ef1a-FLPo",
    "Ef1a-mCherry-IRES-Cre",
    "tdTomato",
]

for gene in gene_list:
    create_umap = CreateUMAP(
        source_path=all_cells_output_dir / "cluster" / "adata.h5ad",
        output_dir=figure_save_dir / f"{experiment_id}_{brain_region}_all_cells",
        n_neighbors=n_neighbors,
        leiden_resolution=leiden_resolution,
        color=gene,
        save_adata=False,
    )
    results = create_umap()
