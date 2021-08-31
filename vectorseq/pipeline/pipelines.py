from pathlib import Path
import numpy as np
from vectorseq.pipeline.stages.reformat import Reformat
from vectorseq.pipeline.stages.distribution_plots import DistributionPlots
from vectorseq.pipeline.stages.filter import Filter
from vectorseq.pipeline.stages.normalize import Normalize
from vectorseq.pipeline.stages.cluster import Cluster
from vectorseq.pipeline.stages.cluster_metrics import ClusterMetrics
from vectorseq.pipeline.stages.expression_plots import ExpressionPlots
from vectorseq.pipeline.stages.create_umap import CreateUMAP
from vectorseq.pipeline.stages.subset import Subset
from vectorseq.marker_constants import (
    BrainGenes,
    VentralMidbrainGenes,
    SuperiorColliculusGenes,
    CortexGenes,
)
from vectorseq.utils.file_utils import create_dir
from typing import List, Union


class ReformatPipeline:
    """
    Pipeline for reformatting raw data into AnnData & generating distribution plots
    for assessing data quality and filtering cut-offs.

    If pipeline stage has already been computed, will load from cache & resume.
    """

    def __init__(
        self,
        source_data_path: str = None,
        output_dir: str = None,
        brain_region: str = "",
        num_seq_lanes: int = 4,
    ):
        self.source_data_path = Path(source_data_path)
        self.output_dir = create_dir(output_dir)
        self.stages = [
            Reformat(
                source_path=source_data_path,
                output_dir=output_dir,
                brain_region=brain_region,
                num_seq_lanes=num_seq_lanes,
            ),
            DistributionPlots(
                source_path=output_dir / "reformat" / "adata.h5ad",
                output_dir=output_dir,
            ),
        ]

    def __call__(self):
        # Sequentially execute each stage
        for stage in self.stages:
            # Execute Pipeline Stage
            results = stage()
        # adata with all annotations
        adata = results["adata"]
        return adata


class PreprocessPipeline:
    """
    Pipeline for filtering & normalizing data prior to further analysis.

    If pipeline stage has already been computed, will load from cache & resume.
    """

    def __init__(
        self,
        source_data_path: str = None,
        output_dir: str = None,
        min_genes_per_cell: int = 200,
        min_counts_per_cell: int = 500,
        mito_pct: float = 0.05,
        tg_markers: Union[List[str]] = None,
        gene_markers: Union[List[str]] = None,
    ):
        self.source_data_path = (
            Path(source_data_path)
            if source_data_path
            else Path(output_dir) / "distribution_plots" / "adata.h5ad"
        )
        self.output_dir = create_dir(output_dir)
        self.stages = [
            Filter(
                source_path=self.source_data_path,
                output_dir=output_dir,
                min_genes_per_cell=min_genes_per_cell,
                min_counts_per_cell=min_counts_per_cell,
                mito_pct=mito_pct,
            ),
            Normalize(
                source_path=output_dir / "filter" / "adata.h5ad",
                output_dir=output_dir,
                tg_markers=tg_markers,
                gene_markers=gene_markers,
            ),
        ]

    def __call__(self):
        # Sequentially execute each stage
        for stage in self.stages:
            # Execute Pipeline Stage
            results = stage()
        # adata with all annotations
        adata = results["adata"]
        return adata


class ClusterPipeline:
    """
    Pipeline for clustering cells, generating clustering metrics and visualizations.

    If pipeline stage has already been computed, will load from cache & resume.
    """

    def __init__(
        self,
        source_data_path: str = None,
        output_dir: str = None,
        n_neighbors: List[int] = np.arange(55, 65, 5),
        leiden_resolution: List[float] = np.round(
            np.arange(0.55, 0.65, 0.05), decimals=2
        ),
        metrics: List[str] = [
            "num-clusters",
            "calinski-harabasz",
            "davies-bouldin",
            "xie-beni",
            "silhouette",
            "within-sse",
            "within-variance",
        ],
    ):
        self.source_data_path = (
            Path(source_data_path)
            if source_data_path
            else Path(output_dir) / "normalize" / "adata.h5ad"
        )
        self.output_dir = create_dir(output_dir)
        self.n_neighbors = n_neighbors
        self.leiden_resolution = leiden_resolution
        self.stages = [
            Cluster(
                source_path=self.source_data_path,
                output_dir=output_dir,
                n_neighbors=n_neighbors,
                leiden_resolution=leiden_resolution,
            ),
            ClusterMetrics(
                source_path=output_dir / "cluster" / "adata.h5ad",
                output_dir=output_dir,
                metrics=metrics,
            ),
            CreateUMAP(
                source_path=output_dir / "cluster" / "adata.h5ad",
                output_dir=output_dir,
            ),
        ]

    def __call__(self):
        # Sequentially execute each stage
        for stage in self.stages:
            # Execute Pipeline Stage
            results = stage()
        # adata with all annotations
        adata = results["adata"]
        return adata


class ExpressionPlotsPipeline:
    """
    Pipeline for creating gene expression plots.
    Only accepts single values for `n_neighbors` and `leiden_resolution`
    If pipeline stage has already been computed, will load from cache.
    """

    def __init__(
        self,
        source_data_path: str = None,
        output_dir: str = None,
        n_neighbors: int = 15,
        leiden_resolution: float = 0.6,
        brain_region: str = "",
    ):
        self.source_data_path = (
            Path(source_data_path)
            if source_data_path
            else Path(output_dir) / "cluster" / "adata.h5ad"
        )
        self.output_dir = output_dir
        self.n_neighbors = n_neighbors
        self.leiden_resolution = leiden_resolution
        self.brain_region = brain_region
        if brain_region == "general":
            genes = BrainGenes()
            kwargs = {
                "n_neighbors": n_neighbors,
                "leiden_resolution": leiden_resolution,
                "rank_genes": False,
                "var_names": {
                    "Neuron": genes.NEURON_MARKERS,
                    "Non-Neuron": genes.NON_NEURON_MARKERS,
                    "Excitatory": genes.EXCITATORY_MARKERS,
                    "Inhibitory": genes.INHIBITORY_MARKERS,
                },
            }
        elif brain_region == "sc":
            genes = SuperiorColliculusGenes()
            kwargs = {
                "n_neighbors": n_neighbors,
                "leiden_resolution": leiden_resolution,
                "rank_genes": False,
                "var_names": {
                    "Transgenes": genes.SCRG_TG,
                    "Genes of Interest": genes.SCRG_EXCITATORY_MARKERS,
                },
            }
        elif brain_region == "snr":
            genes = VentralMidbrainGenes()
            kwargs = {
                "n_neighbors": n_neighbors,
                "leiden_resolution": leiden_resolution,
                "rank_genes": False,
                "var_names": {
                    "Transgenes": genes.VENTRAL_MIDBRAIN_TG,
                    "Genes of Interest": genes.VENTRAL_MIDBRAIN_MARKERS,
                },
            }
        elif brain_region == "v1":
            genes = CortexGenes()
            kwargs = {
                "n_neighbors": n_neighbors,
                "leiden_resolution": leiden_resolution,
                "rank_genes": False,
                "var_names": {
                    "Transgenes": genes.V1_TG,
                    "Genes of Interest": genes.V1_MARKERS,
                },
            }
        else:
            raise ValueError("Unknown value for argument `brain_region`.")

        self.stage = ExpressionPlots(
            source_path=self.source_data_path, output_dir=output_dir, **kwargs
        )

    def __call__(self):
        results = self.stage()
        return results["adata"]


class SubsetPipeline:
    """
    Pipeline for extracting subset of cells from data based on values in adata.obs
    Only accepts single values for `n_neighbors` and `leiden_resolution`
    """

    def __init__(
        self,
        source_data_path: str = None,
        output_dir: str = None,
        n_neighbors: int = 15,
        leiden_resolution: float = 0.6,
        include_values: List[str] = [],
    ):
        self.source_data_path = (
            Path(source_data_path)
            if source_data_path
            else Path(output_dir) / "cluster" / "adata.h5ad"
        )
        self.output_dir = output_dir
        self.n_neighbors = n_neighbors
        self.leiden_resolution = leiden_resolution
        self.include_values = include_values
        self.stage = Subset(
            source_path=self.source_data_path,
            output_dir=output_dir,
            n_neighbors=n_neighbors,
            leiden_resolution=leiden_resolution,
            include_values=include_values,
        )

    def __call__(self):
        results = self.stage()
        return results["adata"]
