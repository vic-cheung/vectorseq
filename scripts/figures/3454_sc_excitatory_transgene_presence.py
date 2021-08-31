#%%
import scanpy as sc
from pathlib import Path
from vectorseq.utils import create_dir
from vectorseq.pipeline.stages import Cluster, CreateUMAP

data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3454"
brain_region = "sc"
figure_save_dir = create_dir(data_dir / "figures")
figure_save_subdir = f"{experiment_id}_{brain_region}_excitatory_subset"
# %%
adata = sc.read_h5ad(
    figure_save_dir / figure_save_subdir / "expression_plots" / "adata.h5ad"
)

# Set Random Seed
random_seed = 0
# Choose clustering scheme to use
n_neighbors = 55
leiden_resolution = 0.55
clustering_scheme_name = f"neighbors_{n_neighbors}_leiden_{leiden_resolution}"
# Choose cluster to visualize
cluster_ids = [str(item) for item in [0, 13, 18, 21, 4, 11, 10, 6]]

for cluster_id in cluster_ids:
    group_adata = adata[(adata.obs[clustering_scheme_name] == cluster_id), :]
    n_pcs = 50 if group_adata.n_obs > 50 else group_adata.n_obs - 1
    cluster = Cluster(
        adata=group_adata,
        output_dir=figure_save_dir / figure_save_subdir,
        output_subdir=f"cluster_{cluster_id}",
        n_neighbors=n_neighbors,
        leiden_resolution=leiden_resolution,
    )
    results = cluster()
    group_adata = results["adata"]

    create_umap = CreateUMAP(
        adata=group_adata,
        output_dir=figure_save_dir / figure_save_subdir / f"cluster_{cluster_id}",
        output_subdir=f"umap_{cluster_id}",
        n_neighbors=n_neighbors,
        leiden_resolution=leiden_resolution,
        save_adata=False,
    )
    results = create_umap()
    group_adata = results["adata"]

    # Save UMAP with Transgene Present
    create_umap = CreateUMAP(
        adata=group_adata,
        output_dir=figure_save_dir / figure_save_subdir / f"cluster_{cluster_id}",
        output_subdir=f"umap_{cluster_id}",
        n_neighbors=n_neighbors,
        leiden_resolution=leiden_resolution,
        color="transgene_present",
        save_adata=False,
    )
    results = create_umap()
    group_adata = results["adata"]

#%%
