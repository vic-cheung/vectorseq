# VECTORseq analysis

Analysis pipelines and scripts used for the paper:

> *Virally Encoded Connectivity Transgenic Overlay RNA sequencing (VECTORseq) defines projection neurons involved in sensorimotor integration.*  Victoria Cheung; Philip Chung; Max Bjorni; Varvara Shvareva; Yesenia Lopez; Evan H Feinberg

## Structure

* Data from our experiments is not packaged in this repo and must be downloaded separately.  [Download data here](https://ucsf.app.box.com/v/vectorseq-data).  Data is organized by sequencing experiment number, then by brain region.
  * `raw_sequencing_output`: raw data from gene sequencer are located in `<experiment_id>/<brain_region>/fastq`.
  * `transgene_sequences`: sequences used to align transgenes.  `RC` files include all transgene sequences as well as their reverse complements.
  * `vectorseq-data`: sparse matricies with cell count tables after raw_sequencing_output is aligned.  This is the starting point for our data analysis pipelines.
* `vectorseq` folder contains library code.  It can be installed using pip or conda (instructions below).
* `scripts` folder contains data analysis scripts organized by experiment number, which parallel structure of data folders.  Location of data folder should be specified in these scripts.  Outputs from these scripts will be generated in the data folders.

  * __3250:__ Sequencing run for primary visual cortex (V1).
    * `all_cells_analysis.py`: data cleaning, normalization, clustering, and overlaying of broad gene marker categories (e.g. excitatory neurons, inhibitory neurons, non-neuronal cells) in defined clusters
  * __3382:__ Sequencing run for ventral midbrain (SNr) and adjacent structures
    * `all_cells_analysis.py`: data cleaning, normalization, clustering, and overlaying of broad gene marker categories (e.g. excitatory neurons, inhibitory neurons, non-neuronal cells) in defined clusters
    * `inhibitory_analysis.py`: subset data using inhibitory gene marker categories, then re-clustering on inhibitory subset to identify cell populations in inhibitory neurons. Dependent on cluster IDs generated from `all_cells_analysis.py`.
  * __3454:__ Sequencing run for superior colliculus (SC)
    * `all_cells_analysis.py`: data cleaning, normalization, clustering, and overlaying of broad gene marker categories (e.g. excitatory neurons, inhibitory neurons, non-neuronal cells) in defined clusters
    * `excitatory_analysis.py`: subset data using excitatory gene marker categories, then re-clustering on excitatory subset to identify cell populations in excitatory neurons. Dependent on cluster IDs generated from `all_cells_analysis.py`.
  * __figures:__ scripts to generate tables and figures in the paper.

## Pipelines
Reusable data analysis pipeline stages are present in the vectorseq package and used across multiple scripts.  These are found in `vectorseq.pipeline.stages`.  Pre-configured pipelines are in `vectorseq.pipeline.pipelines`.

| Pipeline Stages | Description |
| --- | --- |
| reformat | Converts cell count table in .h5 file to AnnData format used by ScanPy |
| distribution_plots | Gene & counts statistics, % mitochondrial genes & % ribosomal genes |
| filter | Remove cells and genes based on thresholds obtained from distribution_plots |
| normalize | Count normalize, log normalize, TF-IDF, select top k highly variable genes, move Transgenes out of expression data into annotations |
| cluster | PCA, Generate leiden clusters *(optional grid search n_neighbors, resolution)* |
| cluster_metrics | Internal cluster validation metrics & plots |
| create_umap | Visualize clusters with UMAP plots |
| expression_plots| Gene expression plots |
| subset | Subset AnnData based on specific clusters |

## Installing `vectorseq`

### Configure Python Development Environment

Install conda package manager.

```sh
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
```

Create conda environment and install necessary packages.  Replace `<env_name>`.

```sh
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --name <env_name> jupyter ipykernel nb_conda anndata==0.7.5 scanpy==1.7.2 leidenalg pysam pynndescent pandas==1.2.3 numpy scipy pytz matplotlib tqdm black flake8 scikit-learn pyarrow fastparquet snappy seaborn
```

How to export package dependencies.

```sh
conda list -e > requirements.txt
```

How to recreate environment from `requirements.txt` file.  Replace `<env_name>`.

```sh
conda create --name <env_name> --file requirements.txt
```

### Make `vectorseq` an importable package

Reusable code for this project is in the src package.  Build/install the package locally in editable mode.

Navigate to vectorseq folder.  To install:

```sh
pip install .
```

To install in editable mode:

```sh
pip install -e .
```

This creates a package in your environment called `vectorseq`.

To uninstall:

```sh
pip uninstall vectorseq
```
