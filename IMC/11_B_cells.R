# Overview ----------------------------------------------------------------

# 11_B_cells
# TODO

library(imcRtools)
library(dplyr)
library(cytomapper)
library(stringr)
library(ggplot2)
library(CATALYST)
library(dittoSeq)
library(viridis)
library(patchwork)
library(corrplot)
library(tidyverse)
library(ggrepel)
library(EBImage)
library(scuttle)
library(mclust)
library(scater)
library(scran)
library(dichromat)
library(caret)
library(doParallel)
library(ggridges)
library(harmony)
library(Rphenograph)
library(igraph)
library(pheatmap)
library(batchelor)
library(bluster)
library(BiocParallel)

source("config/globals.R")
source("functions/avg_exp.R")

DIRECTORIES$OUTPUT_BCELLS <- file.path(DIRECTORIES$OUTPUT, '11_B_cells')
dir.create(DIRECTORIES$OUTPUT_BCELLS, showWarnings = F)

# Load ----------------------------------------------------------------------
spe <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'spe_classified_v2.rds'))
images <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'images.rds'))
masks <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'masks.rds'))

# B cells -------------------------------------------------------------------
bcells <- spe[,spe$cell_type == B_CELLS]

# TODO - check if this was run
# assay(bcells, 'exprs') <- asinh(assay(bcells, 'counts')/0.25)

# Memory checking -----------------------------------------------------------
cytomapperShiny(bcells, masks, images, 'cell_id', 'img_id')

q <- bcells
assay(q, 'exprs_new') <- asinh(assay(q, 'counts')/0.1)

dittoRidgePlot(q, 'IgD', assay='exprs_new', group.by = 'img_id')


# UMAP no correction --------------------------------------------------------
set.seed(101)
bcells <- runUMAP(bcells, subset_row = rowData(spe)$B_cell_marker, 
                  exprs_values = "exprs") 

names(COLOURS$MCD) <- unique(bcells$Metadata_acname)
p <- (plotReducedDim(bcells, "UMAP", colour_by = "Metadata_acname", 
                     point_size = 0.5, point_alpha=0.5) + 
        scale_color_manual(values = COLOURS$MCD) + theme_classic()) | 
  (plotReducedDim(bcells, "UMAP", colour_by = "img_id", 
                  point_size = 0.5, point_alpha=0.5) + theme_classic())
ggsave(file.path(DIRECTORIES$OUTPUT_BCELLS, 'umap_uncorrected.png'), p, width=14)

# Correction harmony ----------------------------------------------------------------
# Harmonise and attach as new reduction
bcells <- runPCA(bcells, 
                 subset_row = rowData(bcells)$B_cell_marker, 
                 exprs_values = 'exprs',
                 ncomponents = 30,
                 BSPARAM = BiocSingular::ExactParam())

ggplot(data.frame(dim = 1:13, percentVar = attr(reducedDim(bcells, 'PCA'), "percentVar")),
       aes(x=dim, y=percentVar)) + geom_line() + geom_point()

set.seed(101)
out <- RunHarmony(bcells, group.by.vars = "img_id", dims.use=1:7)
reducedDim(bcells, "harmony") <- reducedDim(out, 'HARMONY')

# UMAP with correction
set.seed(101)
bcells <- runUMAP(bcells, dimred = "harmony", name = "UMAP_harmony", n_neighbors=20)

p <- (plotReducedDim(bcells, "UMAP_harmony", colour_by = "Metadata_acname",
                     point_size = 0.5, point_alpha=0.5) +
        scale_color_manual(values = COLOURS$MCD) + theme_classic()) |
  (plotReducedDim(bcells, "UMAP_harmony", colour_by = "img_id",
                  point_size = 0.5, point_alpha=0.5) + theme_classic())
ggsave(file.path(DIRECTORIES$OUTPUT_BCELLS, 'umap_corrected_harmony.png'), p, width=14, height=10)


# Output ----------------------------------------------------------------------------
# TODO - work on scaling
plot_list <- multi_dittoDimPlot(bcells, var = PANEL$Actual[PANEL$B_cell_marker], reduction.use = "UMAP_harmony",
                                assay = "exprs", size = 0.2, list.out = TRUE)
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis())
p <- cowplot::plot_grid(plotlist = plot_list, ncol = 5)
ggsave(file.path(DIRECTORIES$OUTPUT_BCELLS, 'umap_features.png'), p, width=10, height=10)

# Clustering ----------------------------------------------------------------
set.seed(101)
clusters <- clusterCells(bcells, 
                         use.dimred = "harmony", 
                         BLUSPARAM = SNNGraphParam(k = 10, 
                                                   cluster.fun = "louvain",
                                                   type = "jaccard"))
bcells$clusters <- clusters
levels(bcells$clusters)

clus_cols <- colorRampPalette(RColorBrewer::brewer.pal(8, 'Set1'))(length(levels(bcells$clusters)))
names(clus_cols) <- levels(bcells$clusters)

p <- dittoHeatmap(bcells[,sample(1:dim(bcells)[2], 1000, replace=F)], 
                  genes=PANEL$Actual[1:35], annot.by = 'clusters', assay = 'exprs')
ggsave(file.path(DIRECTORIES$OUTPUT_BCELLS, 'heatmap_cells_clusters.png'), p)

p <- pheatmap(t(average_expression(bcells, 'clusters', genes=PANEL$Actual[PANEL$B_cell_marker],
                                   scale=T)), 
              breaks = seq(-4, 4, length.out = 101),
              color=colorRampPalette(c('dodgerblue3', 'white', 'firebrick2'))(100))
ggsave(file.path(DIRECTORIES$OUTPUT_BCELLS, 'heatmap_avg_clusters.png'), p)

p <- pheatmap(t(average_expression(bcells, 'clusters', 
                                   genes=PANEL$Actual[!(PANEL$Actual %in% c('CD20', 'TUNEL', 'DNA1', 'DNA2'))], 
                                   scale = F)), 
              color=viridis(100))
ggsave(file.path(DIRECTORIES$OUTPUT_BCELLS, 'heatmap_avg_clusters_unscaled.png'), p)

# Umap clusters
p <- plotReducedDim(bcells, 'UMAP_harmony', colour_by = 'clusters', text_by = 'clusters',
                    point_size=0.5, point_alpha=0.5) + scale_color_manual(values=clus_cols) + theme_classic()
ggsave(file.path(DIRECTORIES$OUTPUT_BCELLS, 'umap_clusters.png'), p)

# Overlays
dir.create(file.path(DIRECTORIES$OUTPUT_BCELLS, "all_clusters"), showWarnings = F)
plotCells(masks, bcells, 'cell_id', 'img_id', 
          colour_by = 'clusters', colour = list(clusters = clus_cols), display = 'single',
          image_title = list(position = "topright",
                             margin = c(5,5),
                             font = 1,
                             cex = 1),
          legend=list(colour_by.legend.cex=3),
          save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_BCELLS, 'all_clusters',
                                                paste0('all_clusters_image', '.png'))))



col <- list(clusters = c('red'))

for(cl in unique(bcells$clusters)){
  names(col$clusters) <- as.character(cl)
  dir.create(file.path(DIRECTORIES$OUTPUT_BCELLS, cl), showWarnings = F)
  plotCells(masks, bcells[,bcells$clusters == cl], 'cell_id', 'img_id', 
            colour_by = 'clusters', colour = col, display = 'single',
            image_title = list(position = "topright",
                               margin = c(5,5),
                               font = 1,
                               cex = 1),
            save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_BCELLS, cl,
                                                  paste0('Cluster_', cl, '.png'))))
  
}

saveRDS(bcells, "data/processed/data_analysis/bcells_clustered.rds")

bcells <- readRDS("data/processed/data_analysis/bcells_clustered.rds")

# Classification ----------------------------------------
dir.create(file.path(DIRECTORIES$OUTPUT_BCELLS, 'classification'), showWarnings = F)
bcells$classification <- 'UNKNOWN'

# 4 and 11?

bcells$classification[bcells$clusters %in% c(17)] <- 'IEDN'
bcells$classification[bcells$clusters %in% c(1,10)] <- 'MZB'
bcells$classification[bcells$clusters %in% c(11,14)] <- 'DN2'
bcells$classification[bcells$clusters %in% c(9)] <- 'Memory'
bcells$classification[bcells$clusters %in% c(16)] <- 'Plasmablast'
bcells$classification[bcells$clusters %in% c(2,4,6)] <- 'Naive'
bcells$classification[bcells$clusters %in% c(8,5,7,13)] <- 'GC/LZ'
bcells$classification[bcells$clusters %in% c(15,12)] <- 'GC/DZ'
bcells$classification[bcells$clusters %in% c(3)] <- 'MZP'

col <- list(classification = c(UNKNOWN = 'white',
                               MZP = 'orange',
                               MZB = 'blue',
                               Memory = 'brown',
                               DN2 = 'green',
                               IEDN = 'purple',
                               `GC/LZ` = 'yellow',
                               Naive = 'cyan',
                               Plasmablast = 'pink',
                               `GC/DZ` = 'red'))

plotCells(masks, bcells, 'cell_id', 'img_id', 
          colour_by = 'classification', colour = col, display = 'single',
          image_title = list(position = "topright",
                             margin = c(5,5),
                             font = 1,
                             cex = 1),
          legend=list(colour_by.legend.cex=3),
          save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_BCELLS, 'classification',
                                                paste0('classification_image', '.png'))))

p <- plotReducedDim(bcells, 'UMAP_harmony', colour_by = 'classification', text_by = 'classification',
                    point_size=0.5, point_alpha=0.5) + scale_color_manual(values=col$classification) + theme_classic()
ggsave(file.path(DIRECTORIES$OUTPUT_BCELLS, 'umap_classification.png'), p)

p <- pheatmap(t(average_expression(bcells, 'classification', genes=PANEL$Actual[PANEL$B_cell_marker],
                                   scale=T)), 
              breaks = seq(-4, 4, length.out = 101),
              color=colorRampPalette(c('dodgerblue3', 'white', 'firebrick2'))(100))
ggsave(file.path(DIRECTORIES$OUTPUT_BCELLS, 'heatmap_avg_classification.png'), p)

saveRDS(bcells, "data/processed/data_analysis/bcells_classified.rds")
