import scanpy as sc
import pandas as pd
from vectorseq.utils import (
    pickling,
    unpickling,
    create_dir,
    merge_rc_data,
    reformat_adata,
)
from vectorseq.marker_constants import BrainGenes
from pathlib import Path
from typing import Union, List


class Reformat:
    """
    Converts raw sequencing data into annotated data `adata`.
    Sense and antisense counts are merged.

    Args:
        source_path: should be a .h5 raw data file output from gene sequencer
        output_dir: directory to create subdir in which resultant files are saved.
        output_subdir: default is "reformat".  Output files are:
            `tg_group_cell_counts.csv`
            `ensembl_gene_dict.pkl`
            `adata.h5ad`
        brain_region: brain region targeted by experiment
        num_seq_lanes: each experiment's cell suspension gets split into
            multiple sequencing lanes on the sequencer.  These become different
            sequencing batches that need to be combined to represent the
            cell/gene info for each experiment.  This value is 1-indexed
            and is usually between 1-8.
        tg_markers: list of transgene markers

    Returns:
        dict with values of all pipeline outputs
    """

    def __init__(
        self,
        source_path: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        output_subdir: str = "reformat",
        brain_region: str = "",
        num_seq_lanes: int = 1,
        tg_markers: Union[List[str]] = None,
        save_adata: bool = True,
        **kwargs
    ):
        self.source_path = source_path
        self.output_dir = Path(output_dir) / output_subdir
        self.brain_region = brain_region
        self.num_seq_lanes = num_seq_lanes
        self.tg_markers = BrainGenes().TG_MARKERS if tg_markers is None else tg_markers
        self.save_adata = save_adata

    def __call__(self):
        self._setup()
        if (
            self.tg_group_cell_counts_path.exists()
            and self.ensembl_gene_dict_path.exists()
            and self.adata_path.exists()
        ):
            print("REFORMAT stage outputs already exist.  Loading existing files.")
            return self.load()
        else:
            print("Running stage: REFORMAT")
            return self._run()

    def _setup(self):
        "Additional run-time pipeline configuration."
        create_dir(self.output_dir)

        self.tg_group_cell_counts_path = self.output_dir / "tg_group_cell_counts.csv"
        self.ensembl_gene_dict_path = self.output_dir / "ensembl_gene_dict.pkl"
        self.adata_path = self.output_dir / "adata.h5ad"

    def load(self):
        return {
            "adata": sc.read_h5ad(self.adata_path),
            "tg_group_cell_counts": pd.read_csv(self.tg_group_cell_counts_path),
            "ensembl_gene_dict": unpickling(self.ensembl_gene_dict_path),
        }

    def _run(self):
        # begin unpacking h5 data
        adata = sc.read_10x_h5(self.source_path)
        ensembl_gene_dict = dict((v, k) for k, v in adata.var.gene_ids.iteritems())
        filtered_tg_list = list(
            adata.var_names[adata.var_names.isin(list(self.tg_markers))]
        )
        adata = merge_rc_data(adata, filtered_tg_list)
        adata.obs["barcode"] = adata.obs.index

        # in my alignments I named this CAV-GFP but I want to rename it to Cav2-GFP
        if "CAV-GFP" in filtered_tg_list:
            as_list = adata.var.index.tolist()
            idx = as_list.index("CAV-GFP")
            as_list[idx] = "Cav2-GFP"
            adata.var.index = as_list
            print("Renamed CAV-GFP to Cav2-GFP.")

        group_cell_count_labels, adata = reformat_adata(
            adata, self.brain_region, self.num_seq_lanes, filtered_tg_list
        )
        print("AnnData reformatted.")
        pickling(ensembl_gene_dict, self.ensembl_gene_dict_path)
        print("Ensembl-gene dictionary saved.")
        group_cell_count_labels.to_csv(self.tg_group_cell_counts_path)
        print("Transgene cell-count table saved.")
        adata.write(filename=self.adata_path)
        print("Reformatted AnnData saved.")
        return {
            "adata": adata,
            "tg_group_cell_counts": group_cell_count_labels,
            "ensembl_gene_dict": ensembl_gene_dict,
        }
