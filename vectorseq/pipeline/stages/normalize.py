import scanpy as sc
import pandas as pd
import numpy as np
from vectorseq.utils import create_dir
from vectorseq.utils.adata_utils import tfidf
from vectorseq.marker_constants import (
    BrainGenes,
    CortexGenes,
    SuperiorColliculusGenes,
    VentralMidbrainGenes,
)
from scipy.sparse import csr_matrix
from pathlib import Path
from typing import Union, List
from anndata._core.anndata import AnnData
from multiprocessing import cpu_count


class Normalize:
    """
    Normalizes gene expressions.  Count normalize, Log normalize & TFIDF.

    Args:
        adata: annotated data object.  If `None`, will instead load adata from
            `source_path` arg.
        source_path: should be a .h5ad file with saved AnnData object
        output_dir: directory to create subdir in which resultant files are saved.
        output_subdir: default is "normalize".
        tg_markers: list of transgene markers
        gene_markers: list or tuple of gene markers to whitelist and make sure
            is included in normalization even if gene is not highly variable
        save_adata: whether to save resultant adata at end of stage

    Returns:
        dict with AnnData object
    """

    def __init__(
        self,
        adata: AnnData = None,
        source_path: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        output_subdir: str = "normalize",
        tg_markers: Union[List[str]] = None,
        gene_markers: Union[List[str]] = None,
        save_adata: bool = True,
        **kwargs
    ):
        self.adata = adata
        self.tg_markers = tg_markers
        self.gene_markers = gene_markers
        self.source_path = source_path
        self.output_dir = Path(output_dir) / output_subdir
        self.save_adata = save_adata

    def __call__(self):
        self._setup()
        if (self.output_dir / "adata.h5ad").exists():
            print("NORMALIZE stage outputs already exist.  Loading existing files.")
            return {
                "adata": sc.read_h5ad(self.output_dir / "adata.h5ad"),
            }
        else:
            print("Running stage: NORMALIZE")
            return self._run()

    def _setup(self):
        "Additional run-time pipeline configuration."
        create_dir(self.output_dir)

        if not self.adata:
            self.adata = sc.read_h5ad(self.source_path)

        if not self.tg_markers:
            self.tg_markers = BrainGenes().TG_MARKERS

        if not self.gene_markers:
            self.gene_markers = list(
                set(
                    [
                        *BrainGenes().all_genes(),
                        *CortexGenes().all_genes(),
                        *SuperiorColliculusGenes().all_genes(),
                        *VentralMidbrainGenes().all_genes(),
                    ]
                )
            )

        sc._settings.ScanpyConfig.n_jobs = cpu_count()

    def _run(self):
        adata = self.adata

        # Normalize total counts so they are comparable among cells, then logmarithmize to find HVG
        sc.pp.normalize_total(adata, target_sum=1e4)
        print("Normalized counts per cell to 1e4")
        sc.pp.log1p(adata)
        print("Log normalized counts")

        # TFIDF Transform, then move Transgene Expressions into adata.obs
        data_frame, df_bool, idf, df_tfidf = tfidf(adata)
        df_tfidf = df_tfidf.sort_index(axis=1)
        # adata_temp used to calculate approximate transgene expressions
        # before moving into adata.obs.  We compute TFIDF, then re-normalize.
        adata_temp = sc.AnnData(df_tfidf)
        adata_temp.var = adata.var.sort_index(axis=0).copy()
        adata_temp.X = csr_matrix(adata_temp.X)
        sc.pp.normalize_total(adata_temp, target_sum=1e4, exclude_highly_expressed=True)
        sc.pp.log1p(adata_temp)
        print("Computed TFIDF on adata with transgenes.")

        # Move Transformed Transgene expression into adata.obs as annotation
        # rather than gene expression.  All adata.X will now be endogenous
        # gene expression.
        filtered_tg_list = list(
            adata.var_names[adata.var_names.isin(list(BrainGenes().TG_MARKERS))]
        )
        adata_temp_tg_indices = (
            pd.DataFrame(
                np.arange(0, len(adata_temp.var)),
                index=adata_temp.var.index,
            )
            .loc[filtered_tg_list]
            .values.squeeze()
        )
        adata.obs = pd.concat(
            [
                adata.obs,
                pd.DataFrame(
                    data=csr_matrix(adata_temp[:, adata_temp_tg_indices].X).todense(),
                    columns=filtered_tg_list,
                    index=adata_temp.obs.index,
                ),
            ],
            axis=1,
        )
        # Remove transgenes from the matrix
        adata_tg_indices = (
            pd.DataFrame(np.arange(0, len(adata.var)), index=adata.var.index)
            .loc[filtered_tg_list]
            .values.squeeze()
        )
        adata = adata[:, ~np.isin(np.arange(len(adata.var_names)), adata_tg_indices)]
        print("Moved transgene expression into adata.obs")

        # Get highly variable genes (endogenous genes only)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        print("Highly variable genes prior to TFIDF determined...")

        # TFIDF on adata.X (endogenous genes only, transgenes moved into adata.obs)
        data_frame, df_bool, idf, df_tfidf = tfidf(adata)
        df_tfidf = df_tfidf.sort_index(axis=1)
        adata_tfidf = sc.AnnData(df_tfidf)
        adata_tfidf.var = adata.var.sort_index(axis=0).copy()
        adata_tfidf.obs = adata.obs.copy()
        adata_tfidf.raw = adata
        adata_tfidf.X = csr_matrix(adata_tfidf.X)
        sc.pp.normalize_total(
            adata_tfidf, target_sum=1e4, exclude_highly_expressed=True
        )
        sc.pp.log1p(adata_tfidf)
        print(
            "Recompute TFIDF & Re-Normalize now that transgenes are moved "
            "into adata.obs and adata.X includes only endogenous genes."
        )

        # Get highly variable genes after TFIDF transformation
        sc.pp.highly_variable_genes(adata_tfidf, n_top_genes=2000)
        print("Highly variable genes after TFIDF determined...")

        # Take Union of raw highly variable genes before & after TFIDF
        adata_variable_genes = adata.var.loc[:, "highly_variable"]
        adata_tfidf_variable_genes = adata_tfidf.var.loc[:, "highly_variable"]
        variable_genes_summed = adata_tfidf_variable_genes | adata_variable_genes
        adata_tfidf.var.loc[:, "highly_variable"] = variable_genes_summed
        # print(adata_tfidf.var.loc[filtered_tg_list, "highly_variable"])
        print("Take Union of HVGs from before and after TFIDF transform...")

        # Keep these reference genes in case HVG eliminates them (Force "Highly Variable" Annotation)
        gene_whitelist = list(
            adata_tfidf.var_names[adata_tfidf.var_names.isin(self.gene_markers)]
        )
        for i, gene in enumerate(gene_whitelist):
            if not adata_tfidf.var.loc[gene_whitelist[i], "highly_variable"]:
                adata_tfidf.var.loc[gene_whitelist[i], "highly_variable"] = True
        print(
            "Markers of interest that are not statistically HVG are manually labeled "
            "highly_variable so that expression can be plotted downstream..."
        )

        # Save TFIDF-transformed adata, only highly-variable genes
        adata_tfidf_high_var = adata_tfidf[:, adata_tfidf.var.loc[:, "highly_variable"]]

        # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
        sc.pp.scale(adata_tfidf_high_var, zero_center=True, max_value=10)
        print(
            "Zero-scaled and clipped gene expressions exceeding 10 standard deviations."
        )

        # Save scaled adata
        if self.save_adata and not (self.output_dir / "adata.h5ad").exists():
            adata_tfidf_high_var.write(filename=self.output_dir / "adata.h5ad")
            print("Save fully normalized and transformed AnnData object.")
        return {"adata": adata_tfidf_high_var}
