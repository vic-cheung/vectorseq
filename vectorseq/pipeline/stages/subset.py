import scanpy as sc
from pathlib import Path
from vectorseq.utils import create_dir
from anndata._core.anndata import AnnData
from typing import Union, List


class Subset:
    """
    Split annotated data into subset based on values in adata.obs.

    Args:
        adata: annotated data object.  If `None`, will instead load adata from
            `source_path` arg.
        source_path: should be a .h5ad file with saved AnnData object
        output_dir: directory to create subdir in which resultant files are saved.
        output_subdir: default is "cluster".
        n_neighbors: name of column in adata.obs to create subset
        leiden_resolution:
        include_values: values in column of adata.obs to include
        save_adata: whether to save resultant adata at end of stage
    """

    def __init__(
        self,
        adata: AnnData = None,
        source_path: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        output_subdir: str = "subset",
        n_neighbors: int = 15,
        leiden_resolution: float = 0.6,
        include_values: Union[str, List[str]] = None,
        save_adata: bool = True,
        **kwargs,
    ):
        self.adata = adata
        self.source_path = source_path
        self.output_dir = Path(output_dir) / output_subdir
        if not n_neighbors:
            raise ValueError("Must specify n_neighbors")
        else:
            self.n_neighbors = n_neighbors
        if not leiden_resolution:
            raise ValueError("Must specify leiden_resolution")
        else:
            self.leiden_resolution = leiden_resolution
        if not include_values:
            raise ValueError(
                "Must specify string or list of string value for `include_values`."
            )
        else:
            self.include_values = (
                include_values if isinstance(include_values, list) else [include_values]
            )
        self.save_adata = save_adata

    def __call__(self):
        self._setup()
        if (self.output_dir / "adata.h5ad").exists():
            print("SUBSET stage outputs already exist.")
            return {"adata": sc.read_h5ad(self.output_dir / "adata.h5ad")}
        else:
            print("Running stage: SUBSET")
            return self._run()

    def _setup(self):
        "Additional run-time pipeline configuration."
        create_dir(self.output_dir)

        if not self.adata:
            self.adata = sc.read_h5ad(self.source_path)

    def _run(self):
        adata = self.adata

        clustering_scheme_name = (
            f"neighbors_{self.n_neighbors}_leiden_{self.leiden_resolution}"
        )
        # Subset Cells
        subset_adata = adata[
            adata.obs[clustering_scheme_name].isin(self.include_values)
        ]
        # Reset all cluster-related annotations; preprocess annotations remain
        subset_adata.obs = subset_adata.obs.drop(
            columns=[
                col
                for col in subset_adata.obs.columns
                if ("leiden" or "neighbors") in col
            ]
        )
        subset_adata.uns = {}
        subset_adata.obsm, subset_adata.obsp, subset_adata.varm = [], [], []
        print(
            f"Adata: {adata.shape} --> Subset of Adata: {subset_adata.shape}.  Cluster annotations reset."
        )
        if self.save_adata and not (self.output_dir / "adata.h5ad").exists():
            subset_adata.write(filename=self.output_dir / "adata.h5ad")
            print("Subset AnnData saved.")
        return {"adata": subset_adata}
