# Dissecting Gene Regulatory Networks Governing Human Cortical Cell Fate
Scripts for analyses used in [Dissecting Gene Regulatory Networks Governing Human Cortical Cell Fate]. Contains scripts for analysis, raw imaging manual cell counting data, and the comparative annotation toolkit GTF version of the RheMac10 annotation used for the study. The raw reads for the macaque data, and processed data for both species can be found in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284197 (released upon publication). 

## Usage

* `preprocessing`: Scripts used to (1) align gene expression, sgRNA and STICR libraries using cellranger, (2) remove ambient RNA in gene count matrix using Cellbender, (3) assign sgRNA using [cellbouncer](https://github.com/nkschaefer/cellbouncer), (4) assign STICR barcodes using a modified [NextClone](https://github.com/cnk113/NextClone) and CloneDetective workflow, (5) individual demultiplexing using Vireo, (6) create scanpy object and load metadata and both barcode assignments, (7) compare and reference map to the Wang et al. 2025 atlas and cell type annotation, (8) calculate knockdown efficiency at the sgRNA and gene level using DEseq2.

* `analysis`: Scripts used to perform (1) trajectory and RNA velocity analysis using scFates and scVelo, (2) compositional analysis between conditions using dcats and milo, (3) differential gene expression testing using DEseq2 and TRADE, (4) Euclidean and energy distance testing using pertpy, (5) lineage coupling and clonal clustering analysis using coSpar and scLiTr.

* `plot`: Scripts used to plot results from above analyses for figure generation.

* `docs`: Files used for cellranger alignment and lists of neurological disorder genes.

* `target-selection`: Scripts used to select targets for this study using SEAcells and scenicplus.






