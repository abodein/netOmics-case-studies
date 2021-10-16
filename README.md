# Interpretation of network-based integration from multi-omics longitudinal data

Multi-omics integration is key to fully understand complex biological processes in an holistic manner. Furthermore, multi-omics combined with new longitudinal experimental design can unreveal dynamic relationships between omics layers and identify key players or interactions in system development or complex phenotypes.

However, integration methods have to address various experimental designs and do not guarantee interpretable biological results. The new challenge of multi-omics integration is to solve interpretation and unlock the hidden knowledge within the multi-omics data.

In this paper [@bodein2020interpretation], we go beyond integration and propose a generic approach to face the interpretation problem.
From multi-omics longitudinal data, this approach builds and explores hybrid multi-omics networks composed of both inferred and known relationships within and between omics layers. With smart node labelling and propagation analysis, this approach predicts regulation mechanisms and multi-omics functional modules.
We applied the method on 3 case studies with various multi-omics designs and identified new multi-layer interactions involved in key biological functions that could not be revealed with single omics analysis.
Moreover, we highlighted interplay in the kinetics that could help identify novel biological mechanisms.
This method is available as an R package \texttt{netOmics} to readily suit any application. 

In this repository, you will find R scripts to conduct both analysis (`CS1_HeLa_Cell_Cycling/`, `CS2_Maize_Aphid_Feeding/`, `CS3_Diabetes_seasonal`).

In each folder, there are:

* `.Rmd` and `.html` files: main scripts and compiled version to be viewed in a browser
* `data/` folder: raw data and used databases
* `figure/` folder: main figures display in the article
* `result/` folder: intermediate `.rda` files and output files (`.csv`/`.xlsx`)

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
