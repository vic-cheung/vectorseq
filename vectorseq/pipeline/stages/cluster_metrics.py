import scanpy as sc
import numpy as np
import pandas as pd
from itertools import product
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score,
)
from vectorseq.utils.metrics import (
    xie_beni_index,
    within_cluster_sum_of_squares,
    within_cluster_variance,
)
from vectorseq.utils.plotting import heatmap_plot, surface_plot_3d
from vectorseq.utils import (
    parallel_process,
    create_dir,
    hyphen_to_underscore,
    extract_neighbors_leiden,
)
from tqdm.auto import tqdm
from anndata._core.anndata import AnnData
from typing import Union, List, Iterable


class ClusterMetrics:
    """
    Compute metrics for cluster hyperparameter sweep.

    Args:
        adata: annotated data object.  If `None`, will instead load adata from
            `source_path` arg.
        source_path: should be a .h5ad file with saved AnnData object
        output_dir: directory to create subdir in which resultant files are saved.
        output_subdir: default is "cluster".
        n_neighbors: number of neighbors used to construct the neighborhood graph.
            The input to this argument can be either a single value or a list
            of values.  If a list of values is passed, this stage will sweep
            all unique combinations of n_neighbors x leiden_resolution and generate
            a plot for each.  Default `None` will sweep all combinations in adata.
        leiden_resolution: resolution for Leiden community detection algorithm.
            The input to this argument can be either a single value or a list
            of values.  If a list of values is passed, this stage will sweep
            all unique combinations of n_neighbors x leiden_resolution and generate
            a plot for each.  Default `None` will sweep all combinations in adata.
        metrics: which metrics to compute.  By default, computes all.  Options are:
            num-clusters, calinski-harabasz, davies-bouldin, xie-beni, silhouette,
            within-sse, within-variance
        save_format: file format for saving plots
        random_seed: random seed for scanpy and numpy
        save_adata: whether to save resultant adata at end of stage
    """

    def __init__(
        self,
        adata: AnnData = None,
        source_path: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        output_subdir: str = "cluster_metrics",
        n_neighbors: Union[int, List[int]] = None,
        leiden_resolution: Union[float, List[float]] = None,
        random_seed: int = 0,
        metrics: List[str] = [
            "num-clusters",
            "calinski-harabasz",
            "davies-bouldin",
            "xie-beni",
            "silhouette",
            "within-sse",
            "within-variance",
        ],
        save_format: str = "svg",
        save_adata: bool = True,
        **kwargs,
    ):
        self.adata = adata
        self.source_path = source_path
        self.output_dir = Path(output_dir) / output_subdir
        self.n_neighbors = n_neighbors
        self.leiden_resolution = leiden_resolution
        self.random_seed = random_seed
        self.save_format = save_format
        self.save_adata = save_adata
        self.metrics = metrics

    def __call__(self):
        self._setup()
        if all([path.exists() for path in self.metric_paths.values()]):
            print("CLUSTER METRICS stage outputs already exist.")
            try:
                d = {
                    "adata": self.adata
                    if self.adata
                    else sc.read_h5ad(self.output_dir / "adata.h5ad"),
                }
            except Exception:
                print(
                    "Cannot load adata.  No file at path: ",
                    self.output_dir / "adata.h5ad",
                )
                d = {}
            # Load metrics from saved filepaths
            for metric, path in self.metric_paths.items():
                d.update({metric: pd.read_parquet(path)})
            return d
        else:
            print("Running stage: CLUSTER METRICS")
            return self._run()

    def _setup(self):
        """
        Additional run-time pipeline configuration.

        If either `n_neighbors` or `leiden_resolution` is `None`, then all values
        for `n_neighbors` and `leiden_resolution` that exist as annotations on
        `adata` will be used.
        """
        create_dir(self.output_dir)

        if not self.adata:
            self.adata = sc.read_h5ad(self.source_path)

        # Infer n_neighbors & leiden_resolution from adata.obs
        if (self.n_neighbors is None) or (self.leiden_resolution is None):
            lst = [
                extract_neighbors_leiden(x)
                for x in self.adata.obs.columns
                if ("neighbors" or "leiden") in x
            ]
            self.n_neighbors = sorted(set([t[0] for t in lst]))
            self.leiden_resolution = sorted(set([t[1] for t in lst]))
        # Ensure specific n_neighbors & leiden_resolution wrapped in list
        else:
            if not isinstance(self.n_neighbors, Iterable):
                self.n_neighbors = [self.n_neighbors]
            if not isinstance(self.leiden_resolution, Iterable):
                self.leiden_resolution = [self.leiden_resolution]

        self.metric_paths = {
            metric: self.output_dir / f"{hyphen_to_underscore(metric)}.parquet"
            for metric in self.metrics
        }

        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        matplotlib.rcParams["svg.fonttype"] = "none"

    def _run(self):
        adata = self.adata

        # Metrics for all combinations of n_neighbors & leiden_resolution in adata
        if (self.n_neighbors[0] or self.leiden_resolution[0]) is None:
            d = self.metrics_and_plots(
                adata,
                random_seed=self.random_seed,
                save_dir=self.output_dir,
                save_format=self.save_format,
            )
        # Metrics only for specific values of n_neighbors & leiden_resolution
        else:
            d = self.metrics_and_plots(
                adata,
                n_neighbors=self.n_neighbors,
                leiden_resolution=self.leiden_resolution,
                random_seed=self.random_seed,
                save_dir=self.output_dir,
                save_format=self.save_format,
            )

        if self.save_adata and not (self.output_dir / "adata.h5ad").exists():
            adata.write(self.output_dir / "adata.h5ad")
            print("Saved AnnData used to compute Metrics & generate Plots.")
        return d

    def metrics_and_plots(
        self,
        adata: AnnData,
        n_neighbors: List[int] = None,
        leiden_resolution: List[float] = None,
        random_seed: int = 0,
        save_dir: str = None,
        save_format: str = "svg",
        **kwargs,
    ):
        """
        Jointly sweeps all combinations of n_neighbors & leiden_resolution.
        This method expects list inputs for `n_neighbors` and `leiden_resolution`.
        """
        metric_fns = {
            "num-clusters": self.compute_num_clusters,
            "calinski-harabasz": self.calinski_harabasz_index,
            "davies-bouldin": self.davies_bouldin_index,
            "xie-beni": self.xie_beni_index,
            "silhouette": self.average_silhouette_score,
            "within-sse": self.within_cluster_sse,
            "within-variance": self.within_cluster_variance,
        }

        d = {"adata": adata}
        for metric in self.metric_paths.keys():
            kwargs = {
                "adata": adata,
                "n_neighbors": n_neighbors,
                "leiden_resolution": leiden_resolution,
                "seed": random_seed,
                "save_dir": save_dir,
                "save_name": hyphen_to_underscore(metric),
                "save_format": save_format,
            }
            d.update({metric: metric_fns[metric](**kwargs)})
        return d

    def compute_num_clusters(
        self,
        adata: AnnData,
        n_neighbors: List[int],
        leiden_resolution: List[float],
        save_dir: Union[str, Path],
        save_name: str = "num_clusters",
        save_format: str = "svg",
        **kwargs,
    ):
        hyperparam_permutations = list(product(n_neighbors, leiden_resolution))
        N, R = len(n_neighbors), len(leiden_resolution)
        heatmap_plot_path = Path(save_dir) / f"{save_name}_heatmap.{save_format}"
        surface_plot_path = Path(save_dir) / f"{save_name}_surface3d.{save_format}"
        data_path = Path(save_dir) / f"{save_name}.parquet"

        if data_path.exists():
            print("Loading Cached Num Clusters")
            df = pd.read_parquet(data_path)
        else:
            print("Computing Num Clusters")
            num_clusters_list = []
            for n, resolution in tqdm(
                hyperparam_permutations, desc="Computing Num Clusters"
            ):
                clustering_scheme_name = f"neighbors_{n}_leiden_{resolution}"
                cluster_ids = adata.obs[clustering_scheme_name].cat.categories
                num_clusters_list += [len(cluster_ids)]
            df = pd.DataFrame(
                data=np.array(num_clusters_list).reshape(N, R),
                index=n_neighbors,
                columns=[f"{resolution}" for resolution in leiden_resolution],
            )
            df.to_parquet(data_path)
        # Heatmap
        heatmap_plot(
            df,
            title="Number of Clusters",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            figsize=(10, 8),
        )
        plt.savefig(heatmap_plot_path)
        # 3D Plot
        surface_plot_3d(
            df,
            title="Number of Clusters",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            zlabel="Number of Clusters",
            azimuth=120,
        )
        plt.savefig(surface_plot_path)
        plt.close("all")
        return df

    def within_cluster_sse(
        self,
        adata: AnnData,
        n_neighbors: List[int],
        leiden_resolution: List[float],
        save_dir: Union[str, Path],
        save_name: str = "within_sse",
        save_format: str = "svg",
        **kwargs,
    ):
        hyperparam_permutations = list(product(n_neighbors, leiden_resolution))
        N, R = len(n_neighbors), len(leiden_resolution)
        heatmap_plot_path = Path(save_dir) / f"{save_name}_heatmap.{save_format}"
        surface_plot_path = Path(save_dir) / f"{save_name}_surface3d.{save_format}"
        data_path = Path(save_dir) / f"{save_name}.parquet"

        if data_path.exists():
            print("Loading Cached Within-Cluster SSE")
            df = pd.read_parquet(data_path)
        else:
            print("Computing Within-Cluster SSE")
            within_cluster_sse_list = parallel_process(
                iterable=(
                    {
                        "X": adata.X,
                        "labels": adata.obs[f"neighbors_{n}_leiden_{resolution}"],
                    }
                    for n, resolution in hyperparam_permutations
                ),
                function=within_cluster_sum_of_squares,
                use_kwargs=True,
                desc="Computing Within-Cluster SSE",
            )
            df = pd.DataFrame(
                data=np.array(within_cluster_sse_list).reshape(N, R),
                index=n_neighbors,
                columns=[f"{resolution}" for resolution in leiden_resolution],
            )
            df.to_parquet(data_path)
        # Heatmap
        heatmap_plot(
            df,
            title="Total Within-Cluster Sum of Squared Error",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            figsize=(36, 8),
            annotate_format=".4E",
        )
        plt.savefig(heatmap_plot_path)
        # 3D Plot
        surface_plot_3d(
            df,
            title="Total Within Cluster Sum of Squared Error",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            zlabel="Total Within-Cluster SSE",
            azimuth=300,
        )
        plt.savefig(surface_plot_path)
        plt.close("all")
        return df

    def within_cluster_variance(
        self,
        adata: AnnData,
        n_neighbors: List[int],
        leiden_resolution: List[float],
        save_dir: Union[str, Path],
        save_name: str = "within_variance",
        save_format: str = "svg",
        **kwargs,
    ):
        hyperparam_permutations = list(product(n_neighbors, leiden_resolution))
        N, R = len(n_neighbors), len(leiden_resolution)
        heatmap_plot_path = Path(save_dir) / f"{save_name}_heatmap.{save_format}"
        surface_plot_path = Path(save_dir) / f"{save_name}_surface3d.{save_format}"
        data_path = Path(save_dir) / f"{save_name}.parquet"

        if data_path.exists():
            print("Loading Cached Within-Cluster Variance")
            df = pd.read_parquet(data_path)
        else:
            print("Computing Within-Cluster Variance")
            within_cluster_variance_list = parallel_process(
                iterable=(
                    {
                        "X": adata.X,
                        "labels": adata.obs[f"neighbors_{n}_leiden_{resolution}"],
                    }
                    for n, resolution in hyperparam_permutations
                ),
                function=within_cluster_variance,
                use_kwargs=True,
                desc="Computing Within-Cluster Variance",
            )
            df = pd.DataFrame(
                data=np.array(within_cluster_variance_list).reshape(N, R),
                index=n_neighbors,
                columns=[f"{resolution}" for resolution in leiden_resolution],
            )
            df.to_parquet(data_path)
        # Heatmap
        heatmap_plot(
            df,
            title="Total Within-Cluster Variance",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            figsize=(36, 8),
            annotate_format=".2f",
        )
        plt.savefig(heatmap_plot_path)
        # 3D Plot
        surface_plot_3d(
            df,
            title="Total Within Cluster Variance",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            zlabel="Total Within-Cluster Variance",
            azimuth=120,
        )
        plt.savefig(surface_plot_path)
        plt.close("all")
        return df

    def calinski_harabasz_index(
        self,
        adata: AnnData,
        n_neighbors: List[int],
        leiden_resolution: List[float],
        save_dir: Union[str, Path],
        save_name: str = "calinski-harabasz",
        save_format: str = "svg",
        **kwargs,
    ):
        """
        Compute Calinski-Harabasz Index (Pseudo-F Statistic)
        Ratio of between-cluster dispersion & within-cluster dispersion.
        Higher values indicates clusters are dense & well-separated.
        """
        hyperparam_permutations = list(product(n_neighbors, leiden_resolution))
        N, R = len(n_neighbors), len(leiden_resolution)
        heatmap_plot_path = Path(save_dir) / f"{save_name}_heatmap.{save_format}"
        surface_plot_path = Path(save_dir) / f"{save_name}_surface3d.{save_format}"
        data_path = Path(save_dir) / f"{save_name}.parquet"

        if data_path.exists():
            print("Loading Cached Calinski-Harabasz Index")
            df = pd.read_parquet(data_path)
        else:
            print("Computing Calinski-Harabasz Index")
            ch_index_list = parallel_process(
                iterable=(
                    {
                        "X": adata.X,
                        "labels": adata.obs[f"neighbors_{n}_leiden_{resolution}"],
                    }
                    for n, resolution in hyperparam_permutations
                ),
                function=calinski_harabasz_score,
                use_kwargs=True,
                desc="Computing Calinski-Harabasz Index",
            )
            df = pd.DataFrame(
                data=np.array(ch_index_list).reshape(N, R),
                index=n_neighbors,
                columns=[f"{resolution}" for resolution in leiden_resolution],
            )
            df.to_parquet(data_path)
        # Heatmap
        heatmap_plot(
            df.astype(int),
            title="Calinski-Harabasz Index (Pseudo-F statistic)",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            figsize=(16, 8),
            annotate_format="d",
        )
        plt.savefig(heatmap_plot_path)
        # 3D Plot
        surface_plot_3d(
            df,
            title="Calinski-Harabasz Index (Pseudo-F Statistic)",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            zlabel="Calinski-Harabasz Index",
            azimuth=300,
        )
        plt.savefig(surface_plot_path)
        plt.close("all")
        return df

    def davies_bouldin_index(
        self,
        adata: AnnData,
        n_neighbors: List[int],
        leiden_resolution: List[float],
        save_dir: Union[str, Path],
        save_name: str = "davies-bouldin",
        save_format: str = "svg",
        **kwargs,
    ):
        """
        Compute Davies-Bouldin Index
        Average "similarity" between clusters, taking into account size of clusters and
        distance between clusters.  Lower values indicate better separation between
        clusters.  Zero is lowest possible score.
        """
        hyperparam_permutations = list(product(n_neighbors, leiden_resolution))
        N, R = len(n_neighbors), len(leiden_resolution)
        heatmap_plot_path = Path(save_dir) / f"{save_name}_heatmap.{save_format}"
        surface_plot_path = Path(save_dir) / f"{save_name}_surface3d.{save_format}"
        data_path = Path(save_dir) / f"{save_name}.parquet"

        if data_path.exists():
            print("Loading Cached Davies-Bouldin Index")
            df = pd.read_parquet(data_path)
        else:
            print("Computing Davies-Bouldin Index")
            db_index_list = parallel_process(
                iterable=(
                    {
                        "X": adata.X,
                        "labels": adata.obs[f"neighbors_{n}_leiden_{resolution}"],
                    }
                    for n, resolution in hyperparam_permutations
                ),
                function=davies_bouldin_score,
                use_kwargs=True,
                desc="Computing Davies-Bouldin Index",
            )
            df = pd.DataFrame(
                data=np.array(db_index_list).reshape(N, R),
                index=n_neighbors,
                columns=[f"{resolution}" for resolution in leiden_resolution],
            )
            df.to_parquet(data_path)
        # Heatmap
        heatmap_plot(
            df,
            title="Davies-Bouldin Index",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            figsize=(16, 8),
            annotate_format=".2f",
        )
        plt.savefig(heatmap_plot_path)
        # 3D Plot
        surface_plot_3d(
            df,
            title="Davies-Bouldin Index",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            zlabel="Davies-Bouldin",
            azimuth=120,
        )
        plt.savefig(surface_plot_path)
        plt.close("all")
        return df

    def xie_beni_index(
        self,
        adata: AnnData,
        n_neighbors: List[int],
        leiden_resolution: List[float],
        save_dir: Union[str, Path],
        save_name: str = "xie-beni",
        save_format: str = "svg",
        **kwargs,
    ):
        """
        Compute Xie-Beni Index
        Measures compactness & separation.  It is the sum of cluster variances
        divided by minimum distance between all cluster centers, normalized by
        number of data points. Smaller values indicate more compact clusters.
        """
        hyperparam_permutations = list(product(n_neighbors, leiden_resolution))
        N, R = len(n_neighbors), len(leiden_resolution)
        heatmap_plot_path = Path(save_dir) / f"{save_name}_heatmap.{save_format}"
        surface_plot_path = Path(save_dir) / f"{save_name}_surface3d.{save_format}"
        data_path = Path(save_dir) / f"{save_name}.parquet"

        if data_path.exists():
            print("Loading Cached Xie-Beni Index")
            df = pd.read_parquet(data_path)
        else:
            print("Computing Xie-Beni Index")
            xb_index_list = parallel_process(
                iterable=(
                    {
                        "X": adata.X,
                        "labels": adata.obs[f"neighbors_{n}_leiden_{resolution}"],
                    }
                    for n, resolution in hyperparam_permutations
                ),
                function=xie_beni_index,
                use_kwargs=True,
                desc="Computing Xie-Beni Index",
            )
            df = pd.DataFrame(
                data=np.array(xb_index_list).reshape(N, R),
                index=n_neighbors,
                columns=[f"{resolution}" for resolution in leiden_resolution],
            )
            df.to_parquet(data_path)
        # Heatmap
        heatmap_plot(
            df,
            title="Xie-Beni Index",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            figsize=(16, 8),
            annotate_format=".2f",
        )
        plt.savefig(heatmap_plot_path)
        # 3D Plot
        surface_plot_3d(
            df,
            title="Xie-Beni Index",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            zlabel="Xie-Beni",
            azimuth=60,
        )
        plt.savefig(surface_plot_path)
        plt.close("all")
        return df

    def average_silhouette_score(
        self,
        adata: AnnData,
        n_neighbors: List[int],
        leiden_resolution: List[float],
        save_dir: Union[str, Path],
        save_name: str = "xie-beni",
        save_format: str = "svg",
        seed: int = 0,
        **kwargs,
    ):
        """
        Compute Average Silhouette Score
        For each data point, computes a score that takes into account closeness of
        points within its cluster, and also the distance to the next nearest cluster.
        Average silhouette for all data points is the final score.
        Higher average silhouette score indicates better defined clusters.
        Ranges from -1 to +1.
        note: computationally intense; sklearn implementation already multi-core.
        """
        hyperparam_permutations = list(product(n_neighbors, leiden_resolution))
        N, R = len(n_neighbors), len(leiden_resolution)
        heatmap_plot_path = Path(save_dir) / f"{save_name}_heatmap.{save_format}"
        surface_plot_path = Path(save_dir) / f"{save_name}_surface3d.{save_format}"
        data_path = Path(save_dir) / f"{save_name}.parquet"

        if data_path.exists():
            print("Loading Cached Average Silhouette Score")
            df = pd.read_parquet(data_path)
        else:
            print("Computing Average Silhouette Score")
            avg_silhouette_score_list = []
            for n, resolution in tqdm(
                hyperparam_permutations, desc="Computing Average Silhouette Score"
            ):
                clustering_scheme_name = f"neighbors_{n}_leiden_{resolution}"
                avg_silhouette_score_list += [
                    silhouette_score(
                        adata.X, adata.obs[clustering_scheme_name], random_state=seed
                    )
                ]
            df = pd.DataFrame(
                data=np.array(avg_silhouette_score_list).reshape(N, R),
                index=n_neighbors,
                columns=[f"{resolution}" for resolution in leiden_resolution],
            )
            df.to_parquet(data_path)
        # Heatmap
        heatmap_plot(
            df,
            title="Mean Silhouette Coefficient",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            figsize=(24, 8),
            annotate_format=".4f",
        )
        plt.savefig(heatmap_plot_path)
        # 3D Plot
        surface_plot_3d(
            df,
            title="Mean Silhouette Coefficient",
            xlabel="Leiden Resolution",
            ylabel="Number of Neighbors",
            zlabel="Mean Silhouette Coefficient",
            azimuth=300,
            elevation=45,
        )
        plt.savefig(surface_plot_path)
        plt.close("all")
        return df
