# Interpretation of network-based integration from multi-omics longitudinal data

We propose a generic analytic and integration framework for multi-omics longitudinal datasets that consists of multi-omics kinetic clustering and multi-layer network-based analysis. This framework was successfully applied to two case studies with different experimental designs and omics data collected.
The first case studied transcriptomic and proteomic changes during cell cycle in human HeLa cells, while the second focused on maize transcriptomic and metabolomic response to aphid feeding.
Propagation analysis on multi-layer networks identifies regulatory mechanisms and function prediction during the HeLa cell cycle. 
It was also applied on maize organism to study the dynamic responses to aphid feeding.

In this repository, you will find R scripts to conduct both analysis (`HeLa_Cell_Cycling/` and `Maize_Aphid_Feeding/`).

In each folder, there are:

* an `.Rmd` file: main script
* `data/` folder: raw data and used databases
* `figure/` folder: main figures display in the article
* intermediare `.rda` files

The following R packages are required:

* tidyverse
* mixOmics
* timeOmics
* igraph
* minet
* gprofilet2
* RandomWalkRestartMH
* GO.db
* ComplexHeatmap
* org.Hs.eg.db
