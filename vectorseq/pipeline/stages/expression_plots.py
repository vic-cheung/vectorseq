import scanpy as sc
from pathlib import Path
from itertools import chain
import matplotlib
from vectorseq.utils import create_dir
from anndata._core.anndata import AnnData
from typing import Union, List, Dict, Tuple
from multiprocessing import cpu_count
import warnings


class ExpressionPlots:
    """
    Creates Heatmap, Dot Plot, Matrix Plot & Dendrogram Plots that show gene expression
    of interest for each cluster.

    Args:
        adata: annotated data object.  If `None`, will instead load adata from
            `source_path` arg.
        source_path: should be a .h5ad file with saved AnnData object
        output_dir: directory to create subdir in which resultant files are saved.
        output_subdir: default is "cluster".
        n_neighbors: number of neighbors used to construct the neighborhood graph.
            Used to specify cluster labels.
        leiden_resolution: resolution for Leiden community detection algorithm.
            Used to specify cluster labels.
        rank_genes: if True, will rank genes for each cluster.  if False, will
            show genes in `var_names` argument.
        n_genes: if `rank_genes` is True, specifies k-top ranked genes for each cluster
            to visualize in plot.  No effect if `rank_genes` is False.
        var_names: dict with string key and values that are list of string gene names.
            Key is used to label the group of values on the plot.
        save_name: string to append to saved files.  By default, output filename
            for plots is `{plot_type}_neighbors_{n}_leiden_{resolution}_(ranked_top_{n_genes})`.
            Specifying a save_name will change output filename to `{plot_type}_{save_name}`.
        save_format: file format for saving plots
        cmap: color map for plot
        show_warnings: whether to show warnings
        save_adata: whether to save resultant adata at end of stage

    Returns:
        dict with AnnData object
    """

    def __init__(
        self,
        adata: AnnData = None,
        source_path: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        output_subdir: str = "expression_plots",
        n_neighbors: int = None,
        leiden_resolution: float = None,
        rank_genes: bool = False,
        n_genes: int = None,
        var_names: Dict[str, List] = None,
        save_name: str = None,
        save_format: str = "pdf",
        cmap: str = "viridis",
        show_warnings: bool = False,
        save_adata: bool = True,
        **kwargs,
    ):
        self.adata = adata
        self.source_path = source_path
        self.output_dir = Path(output_dir) / output_subdir
        self.n_neighbors = n_neighbors
        self.leiden_resolution = leiden_resolution
        self.rank_genes = rank_genes
        self.n_genes = n_genes
        self.var_names = var_names
        self.cmap = cmap
        self.show_warnings = show_warnings
        self.save_name = save_name
        self.save_format = save_format
        self.save_adata = save_adata

    def __call__(self):
        self._setup()
        if (
            self.heatmap_path.exists()
            and self.matrixplot_path.exists()
            and self.dotplot_path.exists()
        ):
            print("EXPRESSION PLOTS stage outputs already exist.")
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
            print("Running stage: EXPRESSION PLOTS")
            return self._run()

    def _setup(self):
        "Additional run-time pipeline configuration."
        create_dir(self.output_dir)

        if not self.adata:
            self.adata = sc.read_h5ad(self.source_path)

        # Only Keep Var Names that exist in adata
        if self.var_names:
            if self.show_warnings:
                diff = set(self.adata.var_names) - set(self.var_names)
                if len(diff) > 0:
                    warnings.warn(
                        f"AnnData `adata` does not have genes {list(diff)} in adata.var_names"
                    )

            # Get Transgenes Annotated in adata.obs
            # (recall in normalize stage transgenes moved from adata.var_names to adata.obs)
            self.obs_names_interest = [
                x
                for x in self.adata.obs.columns
                if x in chain(*self.var_names.values())
            ]

            # Only Keep Var_names that exist in adata
            self.var_names_filtered = {
                k: [
                    x
                    for x in v
                    if x in (list(self.adata.var_names) + self.obs_names_interest)
                ]
                for k, v in self.var_names.items()
            }

        self.clustering_scheme_name = (
            f"neighbors_{self.n_neighbors}_leiden_{self.leiden_resolution}"
        )
        # Generate default file save name
        if not self.save_name:
            self.save_name = (
                self.clustering_scheme_name
                if not self.rank_genes
                else f"{self.clustering_scheme_name}_ranked_top_{self.n_genes}"
            )
        self.suffix = f"{self.save_name}.{self.save_format}"
        self.dendrogram_path = self.output_dir / f"dendrogram_{self.suffix}"
        self.heatmap_path = self.output_dir / f"heatmap_{self.suffix}"
        self.matrixplot_path = self.output_dir / f"matrixplot_{self.suffix}"
        self.dotplot_path = self.output_dir / f"dotplot_{self.suffix}"

        sc._settings.ScanpyConfig.n_jobs = cpu_count()
        sc.settings.figdir = create_dir(self.output_dir)
        sc.settings.autoshow = False
        sc.settings.autosave = True
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        matplotlib.rcParams["svg.fonttype"] = "none"

    def _run(self):
        adata = self.adata

        self.dendrogram(
            adata, groupby=self.clustering_scheme_name, save_format=self.save_format
        )

        if self.rank_genes:
            self.rank_genes_and_plot(
                adata,
                groupby=self.clustering_scheme_name,
                n_genes=self.n_genes,
                cmap=self.cmap,
                save_suffix=self.suffix,
            )
        else:
            self.plot(
                adata,
                groupby=self.clustering_scheme_name,
                var_names=self.var_names_filtered,
                cmap=self.cmap,
                save_suffix=self.suffix,
            )

        if self.save_adata and not (self.output_dir / "adata.h5ad").exists():
            adata.write(self.output_dir / "adata.h5ad")
            print("Saved AnnData used to generate Plots.")
        return {"adata": self.adata}

    def dendrogram(self, adata, groupby, save_format="svg"):
        print("Computing & Creating Dendrogram Plot.")
        sc.tl.dendrogram(adata, groupby=groupby, use_raw=False)
        sc.pl.dendrogram(adata, groupby=groupby, save=f"_{self.suffix}")

    def rank_genes_and_plot(
        self,
        adata: AnnData,
        groupby: str,
        n_genes: int = 50,
        cmap: str = "viridis",
        figsize: Tuple[float, float] = None,
        save_suffix: str = None,
    ):
        print("Ranking Genes Clusters.")
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            use_raw=False,
            n_genes=n_genes,
            method="wilcoxon",
            corr_method="benjamini-hochberg",
            tie_correct=True,
            pts=True,
            key_added=f"ranked_{groupby}",
        )
        print("Creating Heatmap.")
        sc.pl.rank_genes_groups_heatmap(
            adata,
            groupby=groupby,
            key=f"ranked_{groupby}",
            use_raw=False,
            n_genes=n_genes,
            dendrogram=True,
            standard_scale="var",
            show_gene_labels=True,
            swap_axes=True,
            cmap=cmap,
            figsize=figsize,
            save=f"_{save_suffix}" if save_suffix else save_suffix,
        )
        print("Creating Matrix Plot.")
        sc.pl.rank_genes_groups_matrixplot(
            adata,
            groupby=groupby,
            key=f"ranked_{groupby}",
            use_raw=False,
            n_genes=n_genes,
            dendrogram=True,
            standard_scale="var",
            swap_axes=True,
            cmap=cmap,
            figsize=figsize,
            save=f"{save_suffix}" if save_suffix else save_suffix,
        )
        print("Creating Dot Plot.")
        sc.pl.rank_genes_groups_dotplot(
            adata,
            groupby=groupby,
            key=f"ranked_{groupby}",
            use_raw=False,
            n_genes=n_genes,
            dendrogram=True,
            standard_scale="var",
            swap_axes=True,
            cmap=cmap,
            figsize=figsize,
            save=f"{save_suffix}" if save_suffix else save_suffix,
        )

    def plot(
        self,
        adata: AnnData,
        groupby: str,
        var_names: Dict[str, str],
        cmap: str = "viridis",
        figsize: Tuple[float, float] = None,
        save_suffix: str = None,
    ):
        print("Creating Heatmap.")
        sc.pl.heatmap(
            adata,
            var_names=var_names,
            groupby=groupby,
            use_raw=False,
            dendrogram=True,
            standard_scale="var",
            show_gene_labels=True,
            swap_axes=True,
            cmap=cmap,
            figsize=figsize,
            var_group_rotation=90,
            save=f"_{save_suffix}" if save_suffix else save_suffix,
        )
        print("Creating Matrix Plot.")
        sc.pl.matrixplot(
            adata,
            var_names=var_names,
            groupby=groupby,
            use_raw=False,
            dendrogram=True,
            standard_scale="var",
            swap_axes=True,
            cmap=cmap,
            figsize=figsize,
            colorbar_title="Relative expression",
            var_group_rotation=90,
            log=False,
            save=f"{save_suffix}" if save_suffix else save_suffix,
        )
        print("Creating Dot Plot.")
        sc.pl.dotplot(
            adata,
            var_names=var_names,
            groupby=groupby,
            use_raw=False,
            dendrogram=True,
            standard_scale="var",
            swap_axes=True,
            cmap=cmap,
            figsize=figsize,
            colorbar_title="Relative expression",
            var_group_rotation=90,
            log=False,
            save=f"{save_suffix}" if save_suffix else save_suffix,
        )
