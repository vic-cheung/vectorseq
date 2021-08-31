from pathlib import Path
import re
import pandas as pd
import numpy as np
import scanpy as sc
from enum import Flag
from functools import reduce
from itertools import compress
from tqdm.auto import tqdm
from anndata._core.anndata import AnnData
from vectorseq.utils.utils import stringify_list, hyphen_to_underscore


def obs_rename(adata: AnnData, sample_number: int, brain_region: str):
    """
    sets index to brainregion_seqlane_cellnumber
    include brain location, sequencing lane, cell number
    e.g. sc_01_00001 to sc_03_30560
    """
    filtered_rows = adata.obs.filter(regex="-" + str(sample_number), axis=0).index
    adata.obs = adata.obs.rename(
        index={
            x: y
            for x, y in zip(
                filtered_rows,
                [
                    f"{brain_region}"
                    + "_"
                    + "{0:02}".format(sample_number)
                    + "_"
                    + "{0:05}".format(cell_number)
                    for cell_number in range(1, len(filtered_rows) + 1)
                ],
            )
        }
    )
    return adata


def uniquify(adata_var: pd.DataFrame):
    """
    For ensembl ids that map to the same gene, make gene names unique by appending
    their ensembl id
    e.g. Ptp4a1 --> Ptp4a1_ENSMUSG00000026064.

    Unique gene index is created, original gene_name preserved in next column

    Args:
        adata.var

    Returns:
        adata.var with unique index for gene names called "unique_names" and
            gene_names are saved in the next column over.
    """
    print("Uniquifying names...")
    adata_var = adata_var.reset_index().rename(columns={"index": "gene_name"})
    duplicate_mask = adata_var.gene_name.duplicated(keep=False)
    num_add = (
        adata_var.loc[duplicate_mask].groupby("gene_name").cumcount().add(1).astype(str)
    )
    renamed_duplicates = [
        f"{gene_name}_{gene_ids}"
        for gene_name, gene_ids in zip(
            adata_var.loc[duplicate_mask].gene_name,
            adata_var.loc[duplicate_mask].gene_ids,
        )
    ]
    new_names = (
        pd.DataFrame(zip(adata_var.loc[duplicate_mask].index, renamed_duplicates))
        .set_index(keys=0)
        .rename(columns={1: "unique_names"})
        .squeeze()
    )
    old_names = pd.DataFrame(adata_var.gene_name)
    old_names.gene_name[duplicate_mask] = new_names
    adata_var["unique_names"] = old_names
    adata_var = adata_var.set_index(keys="unique_names")
    print("Uniquified names.")
    return adata_var


def merge_rc_data(adata: AnnData, tg_list: list):
    """
    merge data from the reverse complement together with the sense strand
    do not use this function in reformat.py if you wish to keep reverse compelement
    count information separate.

    During custom reference generation, I added the -RC suffix to all reverse
    complement strands.

    Args:
        adata: AnnData object
        tg_list: list of transgenes without the -RC suffix indicating reverse
                complement strand.

    Returns:
        adata_cat: AnnData object with Complement and Reverse Complement counts summed
            together per cell.

    """
    tg_rc_list = [f"{tg}-RC" for tg in tg_list]
    adata_vars = adata[:, adata.var.index.isin(tg_list)].var
    adata_rc = adata[:, adata.var.index.isin(tg_rc_list)]
    adata_rc.var_names = [item.removesuffix("-RC") for item in adata_rc.var_names]
    adata_c = adata[:, adata.var.index.isin(tg_list)]
    adata = adata[:, ~adata.var.index.isin(tg_list + tg_rc_list)]
    new_adata = sc.AnnData(adata_c.to_df() + adata_rc.to_df())
    new_adata.var = adata_vars
    adata_cat = (
        adata.copy()
        .T.concatenate(new_adata.copy().T, index_unique=None, batch_key="is_transgene")
        .T
    )
    return adata_cat


def gene_list_to_flag(adata, gene_list):
    """
    Generates Flag enum for all combinations of gene_list that exist in adata,
    and counts all cells for each combination of genes.

    Args:
        adata: AnnData object
        gene_list: tuple or list of genes to tally up counts

    Returns:
        Tuple with first element being a pandas series that contains cell counts
        for all gene combinations in gene_list that are expressed by cells.
        Second element is a pandas dataframe with index of list of cells that
        expressed at least one gene in gene_list and column values of Flag.
        Third element is the enum Flag object.
    """
    # Find all genes in gene_list that exist in adata
    GENE_LIST = list(adata.var_names[adata.var_names.isin(gene_list)])
    # Generate Flags for genes in gene_list
    flags = Flag("Genes", hyphen_to_underscore(GENE_LIST))

    # Dataframe of cells x genes in list that exist in adata
    gene_presence_df = adata[:, adata.var_names.isin(GENE_LIST)].to_df().astype(bool)
    # Melt dataframe to Long Format & filter out (cell, gene) combinations not present
    melted_df = gene_presence_df.melt(ignore_index=False, value_name="present")
    gene_present_df = melted_df[melted_df.present].variable
    # Convert hyphens to underscore (Enum/Flags don't support hyphens in enum values)
    gene_present_df = gene_present_df.apply(hyphen_to_underscore)
    # Convert gene string to Flag
    gene_present_df = gene_present_df.apply(
        lambda gene_str: flags._member_map_[gene_str]
    )
    # Aggregate all gene Flags for each cell
    cell_gene_lists = gene_present_df.reset_index().groupby("index").agg(list)

    def reduce_gene_flag_list(lst: list):
        return reduce(lambda x, y: x | y, lst)

    # Reduce gene Flags list for each cell down to a single Flag value
    cell_gene_flags = cell_gene_lists.variable.apply(
        lambda gene_flag_list: reduce_gene_flag_list(gene_flag_list)
    )

    # For each gene Flag combination, tally up cell count
    gene_flag_cell_counts = cell_gene_flags.value_counts()
    return gene_presence_df, gene_flag_cell_counts, cell_gene_flags, flags


def gene_mask(adata: AnnData, genes_of_interest: str, col_name: str):
    """
    Creates a transgenes boolean mask and then append transgenes data to adata

    Arg:
        adata: AnnData
        genes of interest: str, can be regex expression
        col_name: specify name of the column to be generated
    """
    # make transgenes mask and then append transgenes data to adata
    gene_bool = adata.var.index.str.contains(
        genes_of_interest,
        flags=re.IGNORECASE,
    )
    adata.var[f"{col_name}"] = pd.Series(gene_bool, index=adata.var.index)
    print(f"""items in {col_name} analyzed: {sum(gene_bool)}""")
    return adata, gene_bool


def reformat_adata(
    adata: AnnData, brain_region: str, num_seq_lanes: int, transgenes_list: str
):
    """
    script that takes in user specified inputs in the data_reformat script
    transforms dataframe input to usable AnnData output with group cell count labels,
    df_obs

    it also makes genes in the index since multiple ensembl IDs can map onto the same gene
    """
    for i in range(1, num_seq_lanes + 1):
        adata = obs_rename(adata, i, brain_region)

    obs_seq_lanes_keys = [
        int(seq_lane[1]) for seq_lane in adata.obs.index.str.split("_")
    ]
    obs_seq_lanes_df = pd.DataFrame(
        obs_seq_lanes_keys, index=adata.obs.index, columns=["seq_lane_number"]
    )

    print("Num seq_lanes parsed...")

    # create bit labels for each transgene and its possible combinations.
    gene_presence_df, _, cell_gene_flags, _ = gene_list_to_flag(adata, transgenes_list)
    adata.obs[[col.upper() for col in gene_presence_df.columns]] = gene_presence_df
    adata.obs["which_transgenes"] = cell_gene_flags
    adata.obs["transgene_present"] = (
        adata.obs["which_transgenes"].notnull().astype("str")
    )
    group_cell_count_labels = adata.obs["which_transgenes"].value_counts(dropna=False)
    adata.obs["seq_lane"] = obs_seq_lanes_df
    print("Group cell count labels generated")

    if adata.var.index.has_duplicates:
        print(f"Duplicate gene names in index (T/F): {adata.var.index.has_duplicates}")
        adata.var = uniquify(adata.var)
    else:
        print(f"Duplicate gene names in index (T/F): {adata.var.index.has_duplicates}")

    adata, __ = gene_mask(
        adata, stringify_list(transgenes_list), col_name="transgene_mask"
    )
    adata, ribo_mask = gene_mask(adata, "^rp[sl][0-9]", col_name="ribo_mask")
    adata, mito_mask = gene_mask(adata, "^mt*-", col_name="mito_mask")

    adata.obs["percent_ribo"] = np.sum(adata[:, ribo_mask].X, axis=1) / np.sum(
        adata.X, axis=1
    )
    adata.obs["percent_mito"] = np.sum(adata[:, mito_mask].X, axis=1) / np.sum(
        adata.X, axis=1
    )
    adata.obs = adata.obs.drop(
        columns=adata.obs.columns[adata.obs.columns.str.contains("temp")]
    )
    return (group_cell_count_labels, adata)


def tfidf(adata: AnnData):
    """
    My implementation of the tf-idf idea and using probabilities as a cleaner way to
    cluster cells. This may help better classify cells regardless of the fact that a
    fraction of a different cell could contaminate it.
    """
    data_frame = adata.to_df()
    df_bool = data_frame.astype(bool)
    idf = df_bool.sum(axis=0).apply(lambda x: np.log(len(df_bool.index) / (x + 1))) + 1
    df_tfidf = data_frame.multiply(idf, axis=1)
    return (data_frame, df_bool, idf, df_tfidf)


def adata_split(
    adata,
    resolution: str = "leiden_0.6",
    save_folder: Path = Path(
        "/home/victoria/Feinberg-VC-3454-analysis/sc-rg/data/processed/adata_clusters"
    ),
):
    "split adata clusters ids into separate files"
    obs = adata.obs.loc[:, resolution].reset_index(drop=True)
    cluster_ids = list(obs.cat.categories)

    for cluster_id in tqdm(cluster_ids, desc="Splitting adata clusters"):
        cluster_row_indices = list(compress(range(len(obs)), obs == cluster_id))
        adata_cluster = adata[cluster_row_indices, :]
        adata_cluster.write(filename=save_folder / f"adata_cluster_{cluster_id}.h5ad")
