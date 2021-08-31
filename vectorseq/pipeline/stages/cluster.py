import scanpy as sc
from pathlib import Path
from itertools import product
from tqdm.auto import tqdm
from vectorseq.utils import create_dir
from anndata._core.anndata import AnnData
from typing import Union, List, Iterable
from multiprocessing import cpu_count


class Cluster:
    """
    Creates neighborhood graph and computes Leiden community clusters.

    Args:
        adata: annotated data object.  If `None`, will instead load adata from
            `source_path` arg.
        source_path: should be a .h5ad file with saved AnnData object
        output_dir: directory to create subdir in which resultant files are saved.
        output_subdir: default is "cluster".
        random_seed: random seed for scanpy and numpy
        n_pcs: number of principal components for PCA
        use_highly_variable: only use highly variable genes for PCA.  Default in
            scanpy is to use them if they have been computed.
        n_neighbors: number of neighbors used to construct the neighborhood graph,
            directly controlling balance between local vs. global structure in the
            data.  Low values emphasize on learning local data structure whereas
            higher values emphasize on learning global data structure.  It is
            analogous to setting the level of fidelity in which we want to
            approximate a high-dimensional manifold/surface.  Both Leiden algorithm
            and UMAP visualization is dependent on this neighborhood graph.
            The input to this argument can be either a single value or a list
            of values.  If a list of values is passed, this stage will sweep
            all unique combinations of n_neighbors x leiden_resolution.
        leiden_resolution: resolution for Leiden community detection algorithm
            that affects the size and number of clusters.  Higher resolutions lead to
            more clusters, lower resolution leads to fewer clusters.
            The input to this argument can be either a single value or a list
            of values.  If a list of values is passed, this stage will sweep
            all unique combinations of n_neighbors x leiden_resolution.
            save_adata: whether to save resultant adata at end of stage

    Returns:
        dict with AnnData object
    """

    def __init__(
        self,
        adata: AnnData = None,
        source_path: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        output_subdir: str = "cluster",
        n_pcs: int = 50,
        use_highly_variable: bool = True,
        n_neighbors: Union[int, List[int]] = 15,
        leiden_resolution: Union[float, List[float]] = 0.6,
        random_seed: int = 0,
        save_adata: bool = True,
        **kwargs,
    ):
        self.adata = adata
        self.source_path = source_path
        self.output_dir = Path(output_dir) / output_subdir
        self.n_pcs = n_pcs
        self.use_highly_variable = use_highly_variable
        self.n_neighbors = n_neighbors
        self.leiden_resolution = leiden_resolution
        self.random_seed = random_seed
        self.save_adata = save_adata

    def __call__(self):
        self._setup()
        if (self.output_dir / "adata.h5ad").exists():
            print("CLUSTER stage outputs already exist.  Loading existing files.")
            return {
                "adata": sc.read_h5ad(self.output_dir / "adata.h5ad"),
            }
        else:
            print("Running stage: CLUSTER")
            return self._run()

    def _setup(self):
        "Additional run-time pipeline configuration."
        create_dir(self.output_dir)

        if not self.adata:
            self.adata = sc.read_h5ad(self.source_path)

        # Ensure specific n_neighbors & leiden_resolution wrapped in list
        if not isinstance(self.n_neighbors, Iterable):
            self.n_neighbors = [self.n_neighbors]
        if not isinstance(self.leiden_resolution, Iterable):
            self.leiden_resolution = [self.leiden_resolution]

        sc._settings.ScanpyConfig.n_jobs = cpu_count()

    def _run(self):
        adata = self.adata

        print("Calculating PCA")
        sc.pp.pca(
            adata,
            n_comps=self.n_pcs,
            svd_solver="arpack",
            use_highly_variable=self.use_highly_variable,
            random_state=self.random_seed,
        )

        print("Computing Neighborhood Graph & Leiden Clusters")
        self.sweep_neighbors(
            adata,
            n_neighbors=self.n_neighbors,
            leiden_resolution=self.leiden_resolution,
            n_pcs=self.n_pcs,
            random_seed=self.random_seed,
        )

        if self.save_adata and not (self.output_dir / "adata.h5ad").exists():
            adata.write(self.output_dir / "adata.h5ad")
            print("Saved AnnData with Neighborhood Graph & Leiden Clusters.")
        return {"adata": self.adata}

    def sweep_neighbors(
        self,
        adata: AnnData,
        n_neighbors: List[int],
        leiden_resolution: List[float],
        n_pcs: int,
        random_seed: int,
    ):
        """
        Jointly sweeps all combinations of n_neighbors & leiden_resolution.
        This method expects list inputs for `n_neighbors` and `leiden_resolution`.
        """
        pbar = tqdm(list(product(n_neighbors, leiden_resolution)))
        for n, resolution in pbar:
            pbar.set_description(f"n_neighbors={n} | resolution={resolution}")
            sc.pp.neighbors(
                adata,
                n_neighbors=n,
                n_pcs=n_pcs,
                method="umap",
                metric="euclidean",
                knn=True,
                key_added=f"neighbors_{n}",
                random_state=random_seed,
            )
            sc.tl.leiden(
                adata,
                resolution=resolution,
                directed=True,
                use_weights=True,
                n_iterations=-1,
                neighbors_key=f"neighbors_{n}",
                key_added=f"neighbors_{n}_leiden_{resolution}",
                random_state=random_seed,
            )
        return adata
