import scanpy as sc
from vectorseq.utils import create_dir
from pathlib import Path
from typing import Union
from anndata._core.anndata import AnnData


class Filter:
    """
    Filter cells in annotated data based on cut-offs selected from inspecting
    distribution plots, filter counts, genes, and mitochondrial DNA.

    Args:
        adata: annotated data object.  If `None`, will instead load adata from
            `source_path` arg.
        source_path: should be a .h5ad file with saved AnnData object
        output_dir: directory to create subdir in which resultant files are saved.
        output_subdir: default is "filter".
        min_genes_per_cell: cut-off for filtering out genes
        min_counts_per_cell: cut-off for filtering out cells
        mito_pct: cut-off for filtering out mitochondrial DNA
        save_adata: whether to save resultant adata at end of stage

    Returns:
        dict with AnnData object
    """

    def __init__(
        self,
        adata: AnnData = None,
        source_path: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        output_subdir: str = "filter",
        min_genes_per_cell: int = 200,
        min_counts_per_cell: int = 500,
        mito_pct: float = 0.05,
        save_adata: bool = True,
        **kwargs,
    ):
        self.adata = adata
        self.source_path = source_path
        self.output_dir = Path(output_dir) / output_subdir
        self.min_genes_per_cell = min_genes_per_cell
        self.min_counts_per_cell = min_counts_per_cell
        self.mito_pct = mito_pct
        self.save_adata = save_adata

    def __call__(self):
        self._setup()
        if (self.output_dir / "adata.h5ad").exists():
            print("FILTER stage outputs already exist.  Loading existing files.")
            return {
                "adata": sc.read_h5ad(self.output_dir / "adata.h5ad"),
            }
        else:
            print("Running stage: FILTER")
            return self._run()

    def _setup(self):
        "Additional run-time pipeline configuration."
        create_dir(self.output_dir)

        if not self.adata:
            self.adata = sc.read_h5ad(self.source_path)

    def _run(self):
        adata = self.adata
        # Filter out cells by gene & cell count
        sc.pp.filter_cells(adata, min_counts=self.min_counts_per_cell)
        sc.pp.filter_cells(adata, min_genes=self.min_genes_per_cell)
        # Filter out cells by % mitochondrial DNA
        adata = adata[adata.obs["percent_mito"] < self.mito_pct, :]

        if self.save_adata and not (self.output_dir / "adata.h5ad").exists():
            adata.write(filename=self.output_dir / "adata.h5ad")
            print("Filtered AnnData saved.")
        return {"adata": adata}
