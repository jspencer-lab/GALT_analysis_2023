# Overview ----------------------------------------------------------------

# 12_B_cells_clustering
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

source("config/globals.R")
source("functions/avg_exp.R")

DIRECTORIES$OUTPUT_FIGURE <- file.path(DIRECTORIES$OUTPUT, '14_FIGURE_2')
dir.create(DIRECTORIES$OUTPUT_FIGURE, showWarnings = F)

# Load ----------------------------------------------------------------------
bcells <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS,  "bcells_classified.rds"))
images <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'images.rds'))
masks <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'masks.rds'))

# Panel A - IMC output ------------------------------------------------------
plotPixels(images[11], colour_by = c('E_Cadherin', 'CD20',  'CD3', 'CD68'), 
           colour = list(E_Cadherin = c('black', 'orange'),
                         CD20 = c('black', 'cyan'),
                         CD3 = c('black', 'magenta'),
                         CD68 = c('black', 'yellow')),
           bcg = list(E_Cadherin = c(0,7,1),
                      CD20 = c(0,9,1),
                      CD3 = c(0,5,1),
                      CD68 = c(0,5,1)),
           save_plot = list(filename=file.path(DIRECTORIES$OUTPUT_FIGURE, 'IMC_lineage.png')),
           image_title = NULL, legend=NULL, scale_bar=NULL)

plotPixels(images[13], colour_by = c('E_Cadherin', 'CD20',  'CD3', 'CD68', 'Ki_67'), 
           colour = list(E_Cadherin = c('black', 'orange'),
                         CD20 = c('black', 'magenta'),
                         CD3 = c('black', 'cyan'),
                         CD68 = c('black', 'yellow'),
                         Ki_67 = c('black', 'white')),
           bcg = list(E_Cadherin = c(0,7,1),
                      CD20 = c(0,7,1),
                      CD3 = c(0,5,1),
                      CD68 = c(0,5,1),
                      Ki_67 = c(0,5,1)),
           display='single',
           save_plot = list(filename=file.path(DIRECTORIES$OUTPUT_FIGURE, 'IMC_lineage_2.png')),
           image_title = NULL, scale_bar=NULL)

plotCells(masks[11], bcells[,bcells$img_id == 'S0768_LMontorsi-Visium_I_s0_a2'], 'cell_id', 'img_id', 
          outline_by = 'img_id',
          colour = list('img_id' = c(`S0768_LMontorsi-Visium_I_s0_a2` = 'red')),
          save_plot = list(filename=file.path(DIRECTORIES$OUTPUT_FIGURE, 'bcell_outlines.png')),
          image_title = NULL, legend=NULL, scale_bar=NULL)

names(masks[11])


cols <- c('brown', 'red', 'orange', 'blue', 'green', 'darkgreen', 'cyan', 'yellow', 'purple')
names(cols) <- sort(unique(bcells$classification))

df <- as.data.frame(reducedDim(bcells, 'UMAP_harmony'))
df$classification <- bcells$classification
p <- ggplot(df, aes(x=V1, y=V2, color=classification)) + geom_point(size=0.5) + scale_color_manual(values=cols) + 
  theme_minimal() + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position ='none', 
        axis.line = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) 
ggsave(file.path(DIRECTORIES$OUTPUT_FIGURE, 'umap_classification_v2.png'), p, width = 8, height=8)

p <- plotReducedDim(bcells, 'UMAP_harmony', colour_by = 'classification', point_size=0.3, point_alpha=0.3) + 
  scale_color_manual(values = cols) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position ='none', 
        axis.line = element_blank(), axis.title = element_blank()) 
ggsave(file.path(DIRECTORIES$OUTPUT_FIGURE, 'umap_classification.png'), p, width = 8, height=8)

assay(bcells, 'visual') <- asinh(assay(bcells, 'counts')/0.2)
plot_list <- multi_dittoDimPlot(bcells, var = PANEL$Actual[PANEL$B_cell_marker], reduction.use = "UMAP_harmony",
                                assay = "visual", size = 0.1, list.out = TRUE)
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis())
p <- cowplot::plot_grid(plotlist = plot_list, ncol = 5)
ggsave(file.path(DIRECTORIES$OUTPUT_FIGURE, 'umap_features.png'), p, width=10, height=10)

plotCells(masks[11], bcells, 'cell_id', 'img_id', 
          colour_by = 'classification',
          colour = list('classification' = cols),
          save_plot = list(filename=file.path(DIRECTORIES$OUTPUT_FIGURE, 'bcell_classification_overlay.png')),
          image_title = NULL, legend=NULL, scale_bar=NULL)

p <- pheatmap(t(average_expression(bcells, 'classification', 
                                   genes=c('E_Cadherin', 'FCRL4', 'T_bet', 'IgD', 'CD45RB', 'BTLA',
                                           'CD1c', 'Ki_67', 'CD38', 'CD27', 'CD11c', 'CD11b', 'IgM'),
                                   scale=T)), 
              breaks = seq(-4, 4, length.out = 101),
              color=colorRampPalette(c('dodgerblue3', 'white', 'firebrick2'))(100),
              cluster_rows = F, cluster_cols = F)
ggsave(file.path(DIRECTORIES$OUTPUT_FIGURE, 'heatmap_avg_classification.pdf'), p, width=4, height=4)



plotPixels(images[11], object = bcells[,bcells$classification == 'DN2'], mask = masks[11],
           img_id = 'img_id', cell_id = 'cell_id',
           outline_by = 'classification',
           colour_by=c('CD27', 'IgD', 'CD11c'),
           colour = list(CD27=c('black', 'green'),
                         IgD=c('black', 'red'),
                         CD11c=c('black', 'magenta'),
                         classification = c(DN2='white')),
           save_plot = list(filename=file.path(DIRECTORIES$OUTPUT_FIGURE, 'DN2_overlay_CD11c.png')),
           image_title = NULL,
           scale_bar = list(length=100, label="", lwidth=10),
           display='single',
           bcg=list(CD27=c(0,5,1),
                    IgD=c(0,7,0.9),
                    CD11c=c(0,5,1))
           )
