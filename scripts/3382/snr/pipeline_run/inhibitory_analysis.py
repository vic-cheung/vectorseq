#%%
import numpy as np
from pathlib import Path
from vectorseq.pipeline import SubsetPipeline, ClusterPipeline, ExpressionPlotsPipeline

data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3382"
brain_region = "snr"
run_dir = data_dir / experiment_id / brain_region
all_cells_output_dir = run_dir / "all_cells"
inhibitory_output_dir = run_dir / "inhibitory"

#%% [markdown]
# ## Subset Inhibitory Clusters
#
# Inhibitory clusters are identified from expression plots based on known
# markers.  These inhibitory clusters are isolated and re-clustered again
# in order to identify specific sub-clusters of inhibitory neurons with
# different gene expression profiles.
#%%
# Choose clustering scheme to use
n_neighbors = 60
leiden_resolution = 0.6

# Define Cluster identity based on gene expression profile & known markers
inhibitory_cluster_ids = [5, 4, 11, 18, 0, 14, 15]
excitatory_cluster_ids = [1, 6, 10, 17, 8, 19, 3, 2]
undetermined_neuron_cluster_ids = [7]
nonneuron_cluster_ids = [9, 12, 13, 16]

# Choose cluster subset to retain
pipeline = SubsetPipeline(
    source_data_path=all_cells_output_dir / "cluster" / "adata.h5ad",
    output_dir=inhibitory_output_dir,
    n_neighbors=n_neighbors,
    leiden_resolution=leiden_resolution,
    include_values=[str(x) for x in inhibitory_cluster_ids],
)
adata = pipeline()

#%% [markdown]
# ## Cluster only Inhibitory Subset

#%%
# Generate Clusters & Metrics w/o silhouette scores. Recommend 8-16 core CPU, 64GB RAM
pipeline = ClusterPipeline(
    source_data_path=inhibitory_output_dir / "subset" / "adata.h5ad",
    output_dir=inhibitory_output_dir,
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
    source_data_path=inhibitory_output_dir / "subset" / "adata.h5ad",
    output_dir=inhibitory_output_dir,
    n_neighbors=np.arange(start=5, stop=105, step=5),
    leiden_resolution=np.round(np.arange(start=0.05, stop=1.25, step=0.05), decimals=2),
    metrics=["davies-bouldin"],
)
adata = pipeline()

#%%
# Compute Average Silhouette Scores.  Recommend 48-96 core CPU, 64GB RAM
# pipeline = ClusterPipeline(
#     source_data_path=inhibitory_output_dir / "subset" / "adata.h5ad",
#     output_dir=inhibitory_output_dir,
#     n_neighbors=np.arange(start=5, stop=105, step=5),
#     leiden_resolution=np.round(np.arange(start=0.05, stop=1.25, step=0.05), decimals=2),
#     metrics=["silhouette"],
# )
# adata = pipeline()

#%% [markdown]
# ## Look at Gene of Interest Expressions for specific n_neighbors & leiden_resolution
#%%
gene_expression_pipe = ExpressionPlotsPipeline(
    output_dir=inhibitory_output_dir,
    n_neighbors=45,
    leiden_resolution=0.6,
    brain_region=brain_region,
)
adata = gene_expression_pipe()


# %%
