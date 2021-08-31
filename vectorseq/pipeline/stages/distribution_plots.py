import scanpy as sc
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re
from vectorseq.utils import create_dir
from pathlib import Path
from typing import Union
from anndata._core.anndata import AnnData


class DistributionPlots:
    """
    Generates distribution plots for number of genes, number of counts, and
    percent mitochondria.  Plots serve as an assessment of sequencing quality
    and can be used to choose filtering cut-offs (e.g. dead cells, doublets, multiplets).

    Plots are only saved as files and not returned.

    Args:
        adata: annotated data object.  If `None`, will instead load adata from
            `source_path` arg.
        source_path: should be a .h5ad file with saved AnnData object
        output_dir: directory to create subdir in which resultant files are saved.
        output_subdir: default is "distplots".
        save_format: file format for saving plots
        save_adata: whether to save resultant adata at end of stage

    Returns:
        dict with AnnData object
    """

    def __init__(
        self,
        adata: AnnData = None,
        source_path: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        output_subdir: str = "distribution_plots",
        save_format: str = "svg",
        save_adata: bool = True,
        **kwargs,
    ):
        self.adata = adata
        self.source_path = source_path
        self.output_dir = Path(output_dir) / output_subdir
        self.save_format = save_format
        self.save_adata = save_adata

    def __call__(self):
        self._setup()
        if (
            self.dist_gene_count_path.exists()
            and self.scatter_gene_count_path.exists()
            and self.scatter_pct_mito_path.exists()
            and self.violin_pct_mito_ribo_path.exists()
        ):
            print("DISTRIBUTION PLOTS stage outputs already exist.")
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
            print("Running stage: DISTRIBUTION PLOTS")
            return self._run()

    def _setup(self):
        "Additional run-time pipeline configuration."
        create_dir(self.output_dir)

        if not self.adata:
            self.adata = sc.read_h5ad(self.source_path)

        self.dist_gene_count_path = self.output_dir / "dist_gene_count.pdf"
        self.scatter_gene_count_path = (
            self.output_dir / f"scatter_gene_count.{self.save_format}"
        )
        self.scatter_pct_mito_path = (
            self.output_dir / f"scatter_pct_mito.{self.save_format}"
        )
        self.violin_pct_mito_ribo_path = (
            self.output_dir / f"violin_pct_mito_ribo.{self.save_format}"
        )

        sc.settings.figdir = create_dir(self.output_dir)
        sc.settings.autoshow = False
        sc.settings.autosave = True
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        matplotlib.rcParams["svg.fonttype"] = "none"

    def _run(self):
        adata = self.adata

        sc.pp.filter_genes(adata, min_cells=1, inplace=True)
        sc.pp.calculate_qc_metrics(adata, inplace=True)

        # Remove malat1, ribosomal genes, mitochondrial genes
        malat1 = adata.var.index.str.fullmatch("malat1", flags=re.IGNORECASE)
        remove = adata.var["mito_mask"] | adata.var["ribo_mask"] | malat1
        adata = adata[:, ~remove]

        # Gene Count Distribution Plots
        with PdfPages(self.dist_gene_count_path) as pdf:
            sns.distplot(adata.obs["total_counts"], kde=True)
            pdf.savefig()
            plt.close()

            sns.distplot(
                adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
                kde=True,
                bins=1000,
            )
            pdf.savefig()
            plt.close()

            sns.distplot(adata.obs["n_genes_by_counts"], kde=True, bins=60)
            pdf.savefig()
            plt.close()

            sns.distplot(
                adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 1000],
                kde=True,
                bins=60,
            )
            pdf.savefig()
            plt.close()
            print("Saved Gene Count Distribution Plots")

        # Gene Count Scatter Plot
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="n_genes_by_counts",
            save=f"_gene_count.{self.save_format}",
        )
        print("Saved Gene Count Scatter Plot")

        # Percent Mito Scatter Plot
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="percent_mito",
            save=f"_pct_mito.{self.save_format}",
        )
        print("Saved Percent Mito Scatter Plot")

        # Percent Mito Ribo Violin Plot
        sc.pl.violin(
            adata,
            ["percent_mito", "percent_ribo"],
            groupby="transgene_present",
            jitter=0.4,
            multi_panel=True,
            save=f"_pct_mito_ribo.{self.save_format}",
        )
        print("Saved Percent Mito Ribo Violin Plot")

        if self.save_adata and not (self.output_dir / "adata.h5ad").exists():
            adata.write(filename=self.output_dir / "adata.h5ad")
            print("AnnData saved.")
        return {"adata": self.adata}
