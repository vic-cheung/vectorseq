#%%
import numpy as np
from pathlib import Path
from vectorseq.pipeline import SubsetPipeline, ClusterPipeline, ExpressionPlotsPipeline

data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3454"
brain_region = "sc"
run_dir = data_dir / experiment_id / brain_region
all_cells_output_dir = run_dir / "all_cells"
excitatory_output_dir = run_dir / "excitatory"

#%% [markdown]
# ## Subset Excitatory Clusters
#
# Excitatory clusters are identified from expression plots based on known
# markers.  These excitatory clusters are isolated and re-clustered again
# in order to identify specific sub-clusters of excitatory neurons with
# different gene expression profiles.
#%%
# Choose clustering scheme to use
n_neighbors = 55
leiden_resolution = 0.55

# Define Cluster identity based on gene expression profile & known markers
inhibitory_cluster_ids = [17, 18, 6, 3, 19, 5, 16, 1, 14]
excitatory_cluster_ids = [7, 20, 9, 26, 23, 4, 13, 15, 10, 11, 25, 0, 22, 8, 2, 12]
nonneuron_cluster_ids = [21, 24]

# Choose cluster subset to retain
pipeline = SubsetPipeline(
    source_data_path=all_cells_output_dir / "cluster" / "adata.h5ad",
    output_dir=excitatory_output_dir,
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    include_values=[str(x) for x in excitatory_cluster_ids],
)
adata = pipeline()

#%% [markdown]
# ## Cluster only Excitatory Subset

#%%
# Generate Clusters & Metrics w/o silhouette scores. Recommend 8-16 core CPU, 64GB RAM
pipeline = ClusterPipeline(
    source_data_path=excitatory_output_dir / "subset" / "adata.h5ad",
    output_dir=excitatory_output_dir,
    n_neighbors=np.arange(start=5, stop=105, step=5),
    leiden_resolution=np.round(np.arange(start=0.05, stop=1.25, step=0.05), decimals=2),
    metrics=[
        "num-clusters",
        "calinski-harabasz",
        "xie-beni",
        "within-sse",
        "within-variance",
    ],
)
adata = pipeline()

# Compute Davies-Bouldin Index.  Recommend 8-16 core CPU, 64GB RAM
pipeline = ClusterPipeline(
    source_data_path=excitatory_output_dir / "subset" / "adata.h5ad",
    output_dir=excitatory_output_dir,
    n_neighbors=np.arange(start=5, stop=105, step=5),
    leiden_resolution=np.round(np.arange(start=0.05, stop=1.25, step=0.05), decimals=2),
    metrics=["davies-bouldin"],
)
adata = pipeline()
#%%
# Compute Average Silhouette Scores.  Recommend 48-96 core CPU, 64GB RAM
# pipeline = ClusterPipeline(
#     source_data_path=excitatory_output_dir / "subset" / "adata.h5ad",
#     output_dir=excitatory_output_dir,
#     n_neighbors=np.arange(start=5, stop=105, step=5),
#     leiden_resolution=np.round(np.arange(start=0.05, stop=1.25, step=0.05), decimals=2),
#     metrics=["silhouette"],
# )
# adata = pipeline()

#%% [markdown]
# ## Look at Gene of Interest Expressions for specific n_neighbors & leiden_resolution
#%%
gene_expression_pipe = ExpressionPlotsPipeline(
    output_dir=excitatory_output_dir,
    n_neighbors=55,
    leiden_resolution=0.55,
    brain_region=brain_region,
)
adata = gene_expression_pipe()
#%%
