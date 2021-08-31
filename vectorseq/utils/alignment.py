import warnings
import pandas as pd
from pathlib import Path
from anndata._core.anndata import AnnData
from scipy.sparse import csr_matrix
import pysam

warnings.simplefilter(action="ignore", category=FutureWarning)


def get_master_transgene_sequence(
    path=Path("/home/victoria/Feinberg-VC-3250-analysis/transgenes_ref.fa"),
):
    """
    grabs transgene name and sequence and puts it into a dict
    """
    with open(file=path, mode="r") as file:
        name = str()
        names = []
        sequences = []
        for line in file:
            # Line is a header
            if line[0] == ">":
                name, length = line.split("\t")
                name = name[1:]
                names.append(name)
            # Line contains part of sequence
            else:
                seq_chunk = (
                    name,
                    line.rstrip(),
                )
                sequences.append(seq_chunk)

    return (
        pd.DataFrame(sequences, columns=["names", "sequences"])
        .groupby("names")["sequences"]
        .apply(lambda x: "".join(x))
        .to_dict()
    )


# %%
def bam_dict(bam_object: pysam.libcalignmentfile.AlignmentFile, genome_name_list: list):
    """
    takes in BAM object and puts it into a dict with genome_name_list
    """
    return {
        gene_name: [aligned_segment for aligned_segment in bam_object.fetch(gene_name)]
        for gene_name in genome_name_list
    }


#%%
def unpack_aligned_segments(
    aligned_segment: pysam.libcalignedsegment.AlignedSegment, gene: str
):
    query_seq = aligned_segment.query_sequence
    ref_pos = aligned_segment.get_reference_positions()
    ref_start = aligned_segment.reference_start
    ref_end = aligned_segment.reference_end
    try:
        barcode = dict(aligned_segment.tags)["CB"]
    except KeyError:
        barcode = dict(aligned_segment.tags)["CR"]
    return {
        "barcode": barcode,
        "gene": gene,
        "query_seq": query_seq,
        "reference_start": ref_start,
        "reference_end": ref_end,
        "reference_positions": ref_pos,
    }


# %%
def rpkm(n_counts: int, gene_length: int, total_n_counts: int):
    return round(n_counts / ((gene_length / 1e3) * (total_n_counts / 1e6)), 4)


#%%
def rpkm_df(
    adata: AnnData, gene_of_interest: str, gene_length_bp: int, is_tg: bool = True
):
    adata_mask = adata.var.index == gene_of_interest
    if is_tg:
        adata_select = adata[
            list(map(lambda ele: ele == "True", adata.obs[gene_of_interest.upper()])), :
        ]
        dense_mat = csr_matrix(adata_select.X[:, adata_mask]).todense().astype(int)
    else:
        adata_select = adata[:, adata_mask]
        dense_mat = csr_matrix(adata_select.X).todense().astype(int)
        dense_mat = dense_mat[dense_mat > 0]

    df = pd.concat(
        [
            adata_select.obs.barcode,
            adata_select.obs.n_counts.astype(int),
        ],
        axis=1,
    )
    df["gene_of_interest"] = gene_of_interest
    goi_counts = pd.DataFrame(
        dense_mat, columns=["gene_counts"], index=adata_select.obs.index
    )

    fraction_goi = pd.DataFrame(
        goi_counts.goi_counts / df.n_counts,
        columns=["fraction_count_goi"],
        index=adata_select.obs.index,
    )

    rpkm = pd.DataFrame(
        dense_mat
        / ((gene_length_bp / 1e3) * (adata.var.total_counts.sum().astype(int) / 1e6)),
        columns=["rpkm"],
        index=adata_select.obs.index,
    )
    count_per_rpkm = pd.DataFrame(
        goi_counts.goi_counts / rpkm.rpkm,
        columns=["count_per_rpkm"],
        index=adata_select.obs.index,
    )
    df = pd.concat([df, goi_counts, fraction_goi, rpkm, count_per_rpkm], axis=1)
    return df


def check_gene_abundance(adata: AnnData, gene_of_interest: str):
    """
    Returns counts for gene of interest per cell, % of gene of interest counts relative
    to total counts per cell, and total counts per cell.

    This function requires gene of interest to be present in adata.var.
    """
    # Select gene of interest in adata
    var_mask = adata.var.index == gene_of_interest
    adata_select = adata[:, var_mask]
    dense_mat = csr_matrix(adata_select.X).todense().A.squeeze().astype(int)
    if not adata_select.var.empty:
        # Keep only positive counts
        adata_select = adata_select[dense_mat > 0]
        dense_mat = dense_mat[dense_mat > 0].T

    # If empty dataframe, then gene doesn't exist in dataset or no cells express it
    if adata_select.var.empty:
        return pd.DataFrame(
            [],
            columns=[
                "barcode",
                "total_counts_per_cell",
                "goi_counts",
                "fraction_count_goi",
                "percent_count_goi",
            ],
        )
    else:
        cell_barcode_ids = adata_select.obs.barcode
        total_counts_per_cell = adata_select.obs.n_counts.astype(int).rename(
            "total_counts_per_cell"
        )
        # Get counts per cell only for gene of interest
        goi_counts = pd.Series(
            dense_mat, name="goi_counts", index=adata_select.obs.index
        )
        fraction_goi = pd.Series(
            goi_counts / total_counts_per_cell,
            name="fraction_count_goi",
            index=adata_select.obs.index,
        )
        percent_goi = (fraction_goi * 100).rename("percent_count_goi")
        return pd.concat(
            [
                cell_barcode_ids,
                total_counts_per_cell,
                goi_counts,
                fraction_goi,
                percent_goi,
            ],
            axis=1,
        )
