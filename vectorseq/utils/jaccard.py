import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import gaussian_kde
from itertools import combinations, compress
from pathlib import Path
from typing import Union
from anndata._core.anndata import AnnData
from vectorseq.utils.math_utils import create_jitter


def jaccard(kde0_x: np.ndarray, kde1_x: np.ndarray, x: np.ndarray):
    "find the jaccard index of two KDEs. np.ndarray  must be of same dimensions"
    set_x = np.minimum(kde0_x, kde1_x)
    union_x = np.maximum(kde0_x, kde1_x)
    area_set_x = np.trapz(set_x, x)
    area_union_x = np.trapz(union_x, x)
    jaccard_index = area_set_x / area_union_x
    area_outer_x = area_union_x - area_set_x
    return area_set_x, set_x, area_union_x, union_x, area_outer_x, jaccard_index


def pct_overlap(adata, gp0, gp1):
    "get the percent of genes that overlap in the top genes that define each group"
    rank_hvg0 = sc.get.rank_genes_groups_df(adata, group=str(gp0))
    rank_hvg1 = sc.get.rank_genes_groups_df(adata, group=str(gp1))
    set_genes = list(set(rank_hvg0.names) & set(rank_hvg1.names))
    union_genes = list(set(rank_hvg0.names.append(rank_hvg1.names).values))
    jaccard_genes_overlap = len(set_genes) / len(union_genes)
    return jaccard_genes_overlap, set_genes, union_genes


def jaccarding_kde2(
    gp0: int,
    gp1: int,
    adata_cluster_folder: Union[
        Path, str
    ] = "/home/victoria/Feinberg-VC-3454-analysis/sc-rg/data/processed/adata_clusters",
    resolution: str = "leiden_0.6",
):
    gp0_path = Path(adata_cluster_folder) / f"adata_cluster_{gp0}.h5ad"
    gp1_path = Path(adata_cluster_folder) / f"adata_cluster_{gp1}.h5ad"

    adata_gp0 = sc.read_h5ad(gp0_path)
    adata_gp1 = sc.read_h5ad(gp1_path)
    genes = adata_gp0.var.index

    jaccard_kde_table = []
    for gene_idx, gene in zip(np.arange(len(genes)), genes):
        # Get gene expression for all cells in each cluster
        x0 = adata_gp0.X[:, gene_idx]
        x1 = adata_gp1.X[:, gene_idx]

        # create jitter. Avoid singular matrix inv covariance error on gene expressions with only one repeated value
        x0 = x0 + create_jitter(x0)
        x1 = x1 + create_jitter(x1)

        kde0 = gaussian_kde(x0, bw_method=0.25)
        kde1 = gaussian_kde(x1, bw_method=0.25)

        xmin = min(x0.min(), x1.min())
        xmax = max(x0.max(), x1.max())

        # add a 20% margin, as the kde is wider than the data
        dx = 0.20 * (xmax - xmin)
        xmin -= dx
        xmax += dx

        x = np.linspace(xmin.item(), xmax.item(), 1000)
        kde0_x = kde0(x)
        kde1_x = kde1(x)

        (
            area_set_x,
            set_x,
            area_union_x,
            union_x,
            area_outer_x,
            jaccard_index,
        ) = jaccard(kde0_x, kde1_x, x)

        row = {
            "gene": gene,
            "gp0": gp0,
            "gp1": gp1,
            "jaccard_index": jaccard_index,
        }
        jaccard_kde_table += [row]

    df = pd.DataFrame(jaccard_kde_table)
    return df


#%%
def jaccarding_kde(
    adata: AnnData,
    gp0: int,
    gp1: int,
    resolution: str = "leiden_0.6",
    union_genes_only: bool = True,
    metrics_only: bool = False,
    pdfpages_filepath: Path = None,
    save_fig: bool = False,
):
    if union_genes_only:
        jaccard_genes_overlap, set_genes, union_genes = pct_overlap(adata, gp0, gp1)
        genes = union_genes
    else:
        genes = adata.var.index
        jaccard_genes_overlap = 1.0
        # set_genes = np.nan

    var = adata.var.reset_index().loc[:, "idx"]
    # all_var_genes = sc.get.rank_genes_groups_df(adata, group=None)
    # var_unique_genes = all_var_genes.names.unique()
    var_idx = list(var.index[var.isin(genes)])

    obs = adata.obs.loc[:, resolution].reset_index(drop=True)
    x0_idx = list(compress(range(len(obs)), obs == str(gp0)))
    x1_idx = list(compress(range(len(obs)), obs == str(gp1)))

    jaccard_kde_table = []
    for idx, gene in zip(np.arange(len(genes)), genes):
        x0 = adata[x0_idx, :].X[:, var_idx[idx]]
        x1 = adata[x1_idx, :].X[:, var_idx[idx]]

        # create jitter. Avoid singular matrix inv covariance error on gene expressions with only one repeated value
        x0 = x0 + create_jitter(x0)
        x1 = x1 + create_jitter(x1)

        kde0 = gaussian_kde(x0, bw_method=0.25)
        kde1 = gaussian_kde(x1, bw_method=0.25)

        xmin = min(x0.min(), x1.min())
        xmax = max(x0.max(), x1.max())

        # add a 20% margin, as the kde is wider than the data
        dx = 0.20 * (xmax - xmin)
        xmin -= dx
        xmax += dx

        x = np.linspace(xmin.item(), xmax.item(), 1000)
        kde0_x = kde0(x)
        kde1_x = kde1(x)

        (
            area_set_x,
            set_x,
            area_union_x,
            union_x,
            area_outer_x,
            jaccard_index,
        ) = jaccard(kde0_x, kde1_x, x)

        if metrics_only:
            row = {
                "gene": gene,
                "gp0": gp0,
                "gp1": gp1,
                "jaccard_index": jaccard_index,
                "jaccard_genes_overlap": jaccard_genes_overlap,
            }
        else:
            row = {
                "gene": gene,
                "gp0": gp0,
                "gp1": gp1,
                "gp0_kde": kde0_x,
                "gp1_kde": kde1_x,
                "linspace_x": x,
                "area_set": area_set_x,
                "area_union": area_union_x,
                "area_outer": area_outer_x,
                "jaccard_index": jaccard_index,
                "jaccard_genes_overlap": jaccard_genes_overlap,
            }
        jaccard_kde_table += [row]

        if save_fig:
            save_kde_fig(
                x,
                kde0_x,
                kde1_x,
                genes,
                idx,
                area_set_x,
                set_x,
                area_union_x,
                union_x,
                jaccard_index,
                pdfpages_filepath,
            )

    df = pd.DataFrame(jaccard_kde_table)
    return df


def save_kde_fig(
    x,
    kde0_x,
    kde1_x,
    gp0,
    gp1,
    genes,
    idx,
    area_set_x,
    set_x,
    area_union_x,
    union_x,
    jaccard_index,
    pdfpages_filepath=None,
):
    with PdfPages(pdfpages_filepath) as pdf:
        plt.plot(x, kde0_x, color="b", label=f"cluster_{gp0}")
        plt.fill_between(x, kde0_x, 0, color="b", alpha=0.2)
        plt.plot(x, kde1_x, color="orange", label=f"cluster_{gp1}")
        plt.fill_between(x, kde1_x, 0, color="orange", alpha=0.2)
        plt.plot(x, set_x, color="r")
        plt.fill_between(
            x,
            set_x,
            0,
            facecolor="none",
            edgecolor="r",
            hatch="xx",
            label="set",
        )
        plt.fill_between(
            x,
            union_x,
            0,
            facecolor="none",
            edgecolor="g",
            hatch="O",
            label="union",
        )

        handles, labels = plt.gca().get_legend_handles_labels()
        labels[2] += f": {area_set_x * 100:.1f}"
        labels[3] += f": {area_union_x * 100:.1f}"
        plt.legend(handles, labels, title="kde")
        plt.title(
            f"Cluster {gp0} vs Cluster {gp1} \n gene: {genes[idx]}\nJaccard Index: {round(jaccard_index,4)}"
        )
        plt.tight_layout()
        # plt.show()
        pdf.savefig()
        plt.close()


# %%
def pairwise_jaccard_scoring(adata, resolution: str = "leiden_0.6"):
    test_list = range(0, len(adata.obs[resolution].unique()))
    res = list(combinations(test_list, 2))
    temp = pd.DataFrame(
        index=range(len(res)),
        columns=[
            "pairwise_clusters",
            "jaccard_index",
            "jaccard_index_mean",
            "jaccard_index_median",
            "jaccard_genes_overlap",
        ],
    )
    for i in range(len(res)):
        df = jaccarding_kde(adata, gp0=res[i][0], gp1=res[i][1])
        temp.pairwise_clusters[i] = res[i]
        temp.jaccard_index[i] = df.jaccard_index
        temp.jaccard_index_mean[i] = df.jaccard_index.mean()
        temp.jaccard_index_median[i] = df.jaccard_index.median()
        temp.jaccard_genes_overlap[i] = df.jaccard_genes_overlap.mean()


#%%
# def jaccard_similarity(group_1: pd.DataFrame, group_2: pd.DataFrame):
#     group_1 = pd.DataFrame(np.random.rand(5, 1), columns=["group_1"])
#     group_2 = pd.DataFrame(np.random.rand(5, 1), columns=["group_2"])

#     total = pd.concat([group_1, group_2], axis=1)

#     delta = pd.DataFrame(
#         total.group_1 - total.group_2, columns=["delta"], index=total.index
#     )
#     total = pd.concat([total, delta], axis=1)
#     delta_abs = pd.DataFrame(
#         [abs(value) for value in total.delta], columns=["delta_abs"], index=total.index
#     )
#     intersect = pd.concat(
#         [total.group_1.loc[total.delta < 0], total.group_2.loc[~(total.delta < 0)]],
#         axis=0,
#     ).sort_index()

#     union = intersect + delta_abs.squeeze()
#     total = pd.concat(
#         [
#             total,
#             delta_abs,
#             pd.DataFrame(intersect, columns=["intersect"]),
#             pd.DataFrame(union, columns=["union"]),
#         ],
#         axis=1,
#     )
#     jaccard_index = pd.DataFrame(intersect / union, columns=["jaccard_index"])
#     total = pd.concat([total, jaccard_index], axis=1)
#     return total

# %%
# def get_kde(
#     adataX_col,
#     x_grid=np.linspace(-10.0, 10.0, 1000),
#     bandwidth=0.2,
#     kernel="epanechnikov",
#     **kwargs,
# ):
#     kde = KernelDensity(bandwidth=bandwidth, kernel=kernel).fit(adataX_col)
#     # score_samples() returns the log-likelihood of the samples
#     # log_pdf = kde.score_samples(adataX_col)
#     # return np.exp(log_pdf)
#     return kde


# %%
# obs = adata.obs.loc[:, "leiden_0.6"].reset_index(drop=True)
# var = adata.var.reset_index().loc[:, "idx"]
# gp0 = 7
# gp1 = 21
# var_idx = list(compress(range(len(var)), var == "Cldn2"))
# x0_idx = list(compress(range(len(obs)), obs == str(gp0)))
# x1_idx = list(compress(range(len(obs)), obs == str(gp1)))

# # var_idx = [4]
# x0 = adata[x0_idx, :].X[:, var_idx[0]].copy()
# x1 = adata[x1_idx, :].X[:, var_idx[0]].copy()

# x0 = x0 + create_jitter(x0)
# x1 = x1 + create_jitter(x1)

# kde0 = gaussian_kde(x0, bw_method=0.25)
# kde1 = gaussian_kde(x1, bw_method=0.25)
# xmin = min(x0.min(), x1.min())
# xmax = max(x0.max(), x1.max())
# dx = 0.20 * (xmax - xmin)  # add a 20% margin to accomodate kde spread
# xmin -= dx
# xmax += dx

# x = np.linspace(xmin, xmax, 1000)
# kde0_x = kde0(x)
# kde1_x = kde1(x)

# set_x = np.minimum(kde0_x, kde1_x)
# union_x = np.maximum(kde0_x, kde1_x)
# area_set_x = np.trapz(set_x, x)
# area_union_x = np.trapz(union_x, x)
# area_outer_x = area_union_x - area_set_x

# ji = area_set_x / area_union_x

# plt.plot(x, kde0_x, color="b", label=f"cluster_{gp0}")
# plt.fill_between(x, kde0_x, 0, color="b", alpha=0.2)
# plt.plot(x, kde1_x, color="orange", label=f"cluster_{gp1}")
# plt.fill_between(x, kde1_x, 0, color="orange", alpha=0.2)
# plt.plot(x, set_x, color="r")
# plt.fill_between(x, set_x, 0, facecolor="none", edgecolor="r", hatch="xx", label="set")
# plt.fill_between(
#     x, union_x, 0, facecolor="none", edgecolor="g", hatch="O", label="union"
# )

# handles, labels = plt.gca().get_legend_handles_labels()
# labels[2] += f": {area_set_x * 100:.1f}"
# labels[3] += f": {area_union_x * 100:.1f}"
# plt.legend(handles, labels, title="kde")
# plt.title(f"cluster_{gp0} vs cluster_{gp1} \n JACCARD INDEX {round(ji,2)}")
# plt.tight_layout()
# plt.show()

# # print(f"jaccard index = {ji}")

# # %%
# adata = sc.read_h5ad(
#     Path(
#         "/home/victoria/Feinberg-VC-3454-analysis/sc-rg/data/processed/adata_sc-rg_neighbors_rankgenes_excitatory_leiden_0.6.h5ad"
#     )
# )
# # %%
# gp0 = 7
# gp1 = 21

# # adata_raw = sc.read_h5ad(
# #     Path(
# #         "/home/victoria/Feinberg-VC-3454-analysis/sc-rg/data/interim/adata_sc-rg_filtered.h5ad"
# #     )
# # )
# # var_idx = list(
# #     compress(range(len(adata_raw.var.index)), adata_raw.var.index == "1600022D10Rik")
# # )
# resolution = "leiden_0.6"

# pdfpages_filepath = (
#     Path("/home/victoria/Feinberg-VC-3454-analysis/sc-rg/figures/2021-0703")
#     / "jaccard_kde_cluster7_cluster21.pdf"
# )
