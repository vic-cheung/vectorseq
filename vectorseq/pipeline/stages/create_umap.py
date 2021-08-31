import scanpy as sc
from pathlib import Path
from itertools import product
from tqdm.auto import tqdm
import matplotlib
from vectorseq.utils import create_dir, extract_neighbors_leiden
from anndata._core.anndata import AnnData
from typing import Union, List, Iterable
from multiprocessing import cpu_count


class CreateUMAP:
    """
    Creates UMAP plots.  Prerequisite is Cluster pipeline stage.

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
        leiden_resolution: resolution for Leiden community detection algorithm
            The input to this argument can be either a single value or a list
            of values.  If a list of values is passed, this stage will sweep
            all unique combinations of n_neighbors x leiden_resolution and generate
            a plot for each.  Default `None` will sweep all combinations in adata.
        random_seed: random seed for scanpy and numpy
        save_format: file format for saving plots
        color: string name of column in adata.obs or adata.var for coloration of umap plot
        save_adata: whether to save resultant adata at end of stage

    Returns:
        dict with AnnData object
    """

    def __init__(
        self,
        adata: AnnData = None,
        source_path: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        output_subdir: str = "umap",
        n_neighbors: Union[int, List[int]] = None,
        leiden_resolution: Union[float, List[float]] = None,
        random_seed: int = 0,
        save_format: str = "svg",
        color: str = None,
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
        self.color = color
        self.save_adata = save_adata

    def __call__(self):
        self._setup()
        if all([path.exists() for path in self.umap_paths.values()]):
            print("UMAP stage outputs already exist.")
            try:
                return {
                    "adata": self.adata
                    if self.adata
                    else sc.read_h5ad(self.output_dir / "adata.h5ad"),
                }
            except Exception:
                print(
                    "Cannot load adata.  No file at path: ",
                    self.output_dir / "adata.h5ad",
                )
                return {}
        else:
            print("Running stage: UMAP")
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
        if (self.n_neighbors or self.leiden_resolution) is None:
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

        self.clustering_scheme_names = [
            f"neighbors_{n}_leiden_{resolution}"
            for n, resolution in product(self.n_neighbors, self.leiden_resolution)
        ]
        # Define save path names
        self.umap_paths = {
            scheme: (
                self.output_dir / f"umap_{scheme}.{self.save_format}"
                if (self.color is None)
                else self.output_dir / f"umap_{scheme}_{self.color}.{self.save_format}"
            )
            for scheme in self.clustering_scheme_names
        }

        sc._settings.ScanpyConfig.n_jobs = cpu_count()
        sc.settings.figdir = create_dir(self.output_dir)
        sc.settings.autoshow = False
        sc.settings.autosave = True
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        matplotlib.rcParams["svg.fonttype"] = "none"

    def _run(self):
        adata = self.adata

        # Plot UMAP for all combinations of n_neighbors & leiden_resolution in adata
        if (self.n_neighbors[0] or self.leiden_resolution[0]) is None:
            self.plot_umap(
                adata, random_seed=self.random_seed, save_format=self.save_format
            )
        # Plot UMAP only for specific values of n_neighbors & leiden_resolution
        else:
            self.plot_umap(
                adata,
                n_neighbors=self.n_neighbors,
                leiden_resolution=self.leiden_resolution,
                random_seed=self.random_seed,
                save_format=self.save_format,
            )
        if self.save_adata and not (self.output_dir / "adata.h5ad").exists():
            adata.write(self.output_dir / "adata.h5ad")
            print("Saved AnnData used to create UMAP.")
        return {"adata": self.adata}

    def plot_umap(
        self,
        adata: AnnData,
        n_neighbors: List[int] = None,
        leiden_resolution: List[float] = None,
        random_seed: int = 0,
        save_format: str = "svg",
        **kwargs,
    ):
        """
        Jointly sweeps all combinations of n_neighbors & leiden_resolution.
        This method expects list inputs for `n_neighbors` and `leiden_resolution`.
        """
        pbar = tqdm(list(product(n_neighbors, leiden_resolution)))
        for n, resolution in pbar:
            pbar.set_description(f"n_neighbors={n} | resolution={resolution}")
            # Compute UMAP projection
            sc.tl.umap(
                adata,
                min_dist=0.5,
                spread=0.8,
                n_components=2,
                neighbors_key=f"neighbors_{n}",
                random_state=random_seed,
            )
            # Plot with Leiden Clusters
            clustering_scheme_name = f"neighbors_{n}_leiden_{resolution}"
            if self.color:
                sc.pl.umap(
                    adata,
                    neighbors_key=f"neighbors_{n}",
                    color=self.color,
                    save=f"_{clustering_scheme_name}_{self.color}.{save_format}",
                    show=False,
                )
            else:
                sc.pl.umap(
                    adata,
                    neighbors_key=f"neighbors_{n}",
                    color=clustering_scheme_name,
                    save=f"_{clustering_scheme_name}.{save_format}",
                    show=False,
                )
        return adata
