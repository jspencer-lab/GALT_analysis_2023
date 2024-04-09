# GALT_analysis_2023
Code for manuscript "Double-negative B cells and DNASE1L3 colocalise with microbiota in gut-associated lymphoid tissue"

Pre-print DOI: 10.1101/2023.08.29.555265

# System requirements

## MacOS requirements
These pipelines were developed on macOS Monterey v12.5.1

## Software requirements
### IMC
#### Python: 
- Conda environment from BodenmillerGroup/ImcSegmentationPipeline. See Usage at https://github.com/BodenmillerGroup/ImcSegmentationPipeline.
#### CellProfiler:
- CellProfiler v4.
#### R:
- cytomapper v1.14.0
- imcRtools v1.8.0
- harmony v1.1.0
- randomForest 4.7-1.1
- cowplot v1.1.1
- bluster v1.12.0
- viridis v0.6.4
- CATALYST v1.26.0
- dittoSeq v1.14.0
- mclust v6.0.0
- scran vv1.30.0
- patchwork v1.1.3
- caret v6.0-94
- tidyverse v2.0.0
- Hmisc v5.1-1
- devtools v2.4.5
- dichromat v2.0-0.1
- Rphenograph v0.99.1.9003


### CyTOF
#### Cytobank: 
- Pre-processing was done on Cytobank (https://mrc.cytobank.org) version 10.3
#### R:
- R packages as per IMC

### RNAScope
#### CellProfiler:
- CellProfiler v4.
#### R:
- cytomapper v1.12.0
- stringr v1.5.0
- imcRtools v1.6.4
- dplyr v1.1.3
- Biostrings v2.68.1
- GenomicAlignments v1.36.0
- scoper v1.3.0
- scater v1.28.0
- batchelor v1.16.0
- Rphenograph v0.99.1.9003
- igraph v1.5.1
- dittoSeq v1.12.1
- tidyr v1.3.0
- tibble v3.2.1
- viridis 0.6.4
- pheatmap v1.0.12

### Visium
#### R:
- dplyr v1.1.4
- ggplot2 v3.4.4
- patchwork v1.2.0
- Seurat v4.4.0
- viridis v0.6.4
- BayesSpace v1.10.1
- Harmony v0.1.1
- ggpubr v0.6.0
- RcmdrMisc v2.9-0
- corrplot v0.92 
