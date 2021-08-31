#%%
import scanpy as sc
from pathlib import Path
from vectorseq.utils import create_dir
from vectorseq.pipeline.stages import Cluster, CreateUMAP
from tqdm.auto import tqdm

data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3382"
brain_region = "snr"
figure_save_dir = create_dir(data_dir / "figures")
figure_save_subdir = f"{experiment_id}_{brain_region}_inhibitory_subset"
# %%
adata = sc.read_h5ad(
    figure_save_dir / figure_save_subdir / "expression_plots" / "adata.h5ad"
)

# Set Random Seed
random_seed = 0
# Choose clustering scheme to use
n_neighbors = 45
leiden_resolution = 0.6
clustering_scheme_name = f"neighbors_{n_neighbors}_leiden_{leiden_resolution}"
# Choose cluster to visualize
cluster_id = "4"

group_adata = adata[(adata.obs[clustering_scheme_name] == cluster_id), :]
n_pcs = 50 if group_adata.n_obs > 50 else group_adata.n_obs - 1
cluster = Cluster(
    adata=group_adata,
    output_dir=figure_save_dir / figure_save_subdir,
    output_subdir=f"cluster_{cluster_id}",
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    n_pcs=n_pcs,
)
results = cluster()
group_adata = results["adata"]

sc.tl.rank_genes_groups(
    group_adata,
    groupby=clustering_scheme_name,
    use_raw=False,
    n_genes=50,
    method="wilcoxon",
    corr_method="benjamini-hochberg",
    tie_correct=True,
    pts=True,
    key_added=f"ranked_{clustering_scheme_name}",
)

#%%
subcluster_ids = tqdm(list(group_adata.obs[clustering_scheme_name].unique()))
for subcluster_id in subcluster_ids:
    subcluster_ids.set_description(f"Subcluster: {subcluster_id}")
    genes_df = sc.get.rank_genes_groups_df(
        group_adata, key=f"ranked_{clustering_scheme_name}", group=subcluster_id
    )
    k = 10
    top_k_genes = tqdm(genes_df.names.head(k))
    for gene in top_k_genes:
        top_k_genes.set_description(f"Gene: {gene}")
        # Save UMAP with cells expressing specific gene colored
        create_umap = CreateUMAP(
            adata=group_adata,
            output_dir=figure_save_dir / figure_save_subdir / f"cluster_{cluster_id}",
            output_subdir=f"subcluster_{subcluster_id}",
            n_neighbors=n_neighbors,
            leiden_resolution=leiden_resolution,
            color=gene,
            save_adata=False,
        )
        create_umap()

#%%
