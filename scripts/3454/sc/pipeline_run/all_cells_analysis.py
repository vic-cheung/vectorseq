#%%
import numpy as np
from pathlib import Path
from vectorseq.pipeline.pipelines import (
    ReformatPipeline,
    PreprocessPipeline,
    ClusterPipeline,
    ExpressionPlotsPipeline,
)

data_dir = Path("/spare_volume/vectorseq-data")
experiment_id = "3454"
brain_region = "sc"
run_dir = data_dir / experiment_id / brain_region
raw_data_path = run_dir / "raw" / "filtered_feature_bc_matrix.h5"
all_cells_output_dir = run_dir / "all_cells"

#%% [markdown]
# ## Process data & Hyperparameter Sweep n_neighbors & leiden_resolution

#%%
# Convert raw data into AnnData format & normalize.  Recommend 8-16 core CPU, 64GB RAM
pipeline = ReformatPipeline(
    source_data_path=raw_data_path,
    output_dir=all_cells_output_dir,
    brain_region=brain_region,
    num_seq_lanes=4,
)
adata = pipeline()

pipeline = PreprocessPipeline(
    output_dir=all_cells_output_dir,
    min_genes_per_cell=200,
    min_counts_per_cell=500,
    mito_pct=0.05,
)
adata = pipeline()

#%%
# Generate Clusters & Metrics w/o silhouette scores. Recommend 8-16 core CPU, 64GB RAM
pipeline = ClusterPipeline(
    output_dir=all_cells_output_dir,
    n_neighbors=np.arange(start=5, stop=105, step=20),
    leiden_resolution=np.round(np.arange(start=0.05, stop=1.25, step=0.2), decimals=2),
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
    output_dir=all_cells_output_dir,
    n_neighbors=np.arange(start=5, stop=105, step=5),
    leiden_resolution=np.round(np.arange(start=0.05, stop=1.25, step=0.05), decimals=2),
    metrics=["davies-bouldin"],
)
adata = pipeline()
#%%
# Compute Average Silhouette Scores.  Recommend 48-96 core CPU, 64GB RAM
# pipeline = ClusterPipeline(
#     output_dir=all_cells_output_dir,
#     n_neighbors=np.arange(start=5, stop=105, step=5),
#     leiden_resolution=np.round(np.arange(start=0.05, stop=1.25, step=0.05), decimals=2),
#     metrics=["silhouette"],
# )
# adata = pipeline()

#%% [markdown]
# ## Look at Gene Expressions for specific n_neighbors & leiden_resolution
#%%
gene_expression_pipe = ExpressionPlotsPipeline(
    output_dir=all_cells_output_dir,
    n_neighbors=15,
    leiden_resolution=0.45,
    brain_region="general",
)
adata = gene_expression_pipe()
#%%
gene_expression_pipe = ExpressionPlotsPipeline(
    output_dir=all_cells_output_dir,
    n_neighbors=15,
    leiden_resolution=0.65,
    brain_region="general",
)
adata = gene_expression_pipe()
#%%
gene_expression_pipe = ExpressionPlotsPipeline(
    output_dir=all_cells_output_dir,
    n_neighbors=55,
    leiden_resolution=0.55,
    brain_region="general",
)
adata = gene_expression_pipe()
#%%
gene_expression_pipe = ExpressionPlotsPipeline(
    output_dir=all_cells_output_dir,
    n_neighbors=60,
    leiden_resolution=0.60,
    brain_region="general",
)
adata = gene_expression_pipe()
#%%
gene_expression_pipe = ExpressionPlotsPipeline(
    output_dir=all_cells_output_dir,
    n_neighbors=70,
    leiden_resolution=0.35,
    brain_region="general",
)
adata = gene_expression_pipe()
#%%
gene_expression_pipe = ExpressionPlotsPipeline(
    output_dir=all_cells_output_dir,
    n_neighbors=95,
    leiden_resolution=0.7,
    brain_region="general",
)
adata = gene_expression_pipe()

#%%
