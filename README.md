# Dissecting Gene Regulatory Networks Governing Human Cortical Cell Fate

Scripts for analyses used in [**Dissecting Gene Regulatory Networks Governing Human Cortical Cell Fate**](https://www.biorxiv.org/content/10.1101/2025.09.23.678137v1). The raw reads for the macaque data, and processed data for both species can be found at GEO GSE284197. Processed human data is browsable through [UCSC broswer](http://cortical-lineage-perturb-44tf.cells.ucsc.edu/). 


## Usage

* `preprocessing`: Scripts used to (1) align gene expression, sgRNA and STICR libraries using cellranger, (2) remove ambient RNA in gene count matrix using Cellbender, (3) assign sgRNA using [cellbouncer](https://github.com/nkschaefer/cellbouncer), (4) assign STICR barcodes using a modified [NextClone](https://github.com/cnk113/NextClone) and CloneDetective workflow, (5) individual demultiplexing using Vireo, (6) create scanpy object and load metadata and both barcode assignments, (7) compare and reference map to the Wang et al. 2025 atlas and cell type annotation, (8) calculate knockdown efficiency at the sgRNA and gene level using DEseq2.

* `analysis`: Scripts used to perform (1) trajectory and RNA velocity analysis, (2) compositional analysis between conditions, (3) differential gene expression testing, (4) Euclidean and energy distance testing, (5) lineage coupling, clonal clustering and pseudotime analysis.

* `plot`: Scripts used to plot results from above analyses for figure generation.

* `docs`: Files used for Cellranger and DRAGEN alignment and lists of neurodevelopmental and neuropsychiatric disorders genes.

* `target-selection`: Scripts used to select targets for peturbation in this study.






