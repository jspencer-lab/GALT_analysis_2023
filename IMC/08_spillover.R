# Overview ----------------------------------------------------------------

# 8_spillover
# TODO

library(imcRtools)
library(dplyr)
library(cytomapper)
library(stringr)
library(ggplot2)
library(CATALYST)
library(dittoSeq)
library(patchwork)
library(corrplot)

source("config/globals.R")

DIRECTORIES$OUTPUT_SPILLOVER <- file.path(DIRECTORIES$OUTPUT, '08_spillover')
dir.create(DIRECTORIES$OUTPUT_SPILLOVER, showWarnings = F)

# Load -------------------------------------------------------------------
cur_sm <- readRDS("reference/spillover_matrix.rds")

ordered_chs <- data.frame(V1 = PANEL$Metal_tag)
ordered_chs$V1 <- paste0(ordered_chs$V1, "Di")
head(ordered_chs)

cur_sm <- rbind(cur_sm, matrix(data=0, nrow=2, ncol=37, dimnames = list(c('Ir191Di', 'Ir193Di'))))
cur_sm['Ir191Di', 'Ir191Di'] <- 1
cur_sm['Ir193Di', 'Ir193Di'] <- 1

adapted_sm <- adaptSpillmat(cur_sm, out_chs = as.character(ordered_chs$V1))

# Compensate data --------------------------------------------------------
spe <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'spe_load.rds'))

rowData(spe)$channel_name <- paste0(rowData(spe)$Metal_tag, "Di")

isotope_list <- CATALYST::isotope_list

spe <- compCytof(spe, cur_sm, 
                 transform = TRUE, cofactor = 1,
                 isotope_list = CATALYST::isotope_list, 
                 overwrite = FALSE)

# Ressign (raw counts become 'uncomp', comp counts are default)
assay(spe, "counts_uncomp") <- assay(spe, "counts") 
assay(spe, "exprs_uncomp") <- assay(spe, "exprs") 
assay(spe, "counts") <- assay(spe, "compcounts") 
assay(spe, "exprs") <- assay(spe, "compexprs") 
assay(spe, "compcounts") <- assay(spe, "compexprs") <- NULL

# Plots ----------------------------------------------------------------------------
# Plot markers against each other
dir.create(file.path(DIRECTORIES$OUTPUT_SPILLOVER, 'Spillover_scatter'), showWarnings = F)

ds_spe <- spe[,sample(1:dim(spe)[2], 20000)]
for(i in 1:36){
  marker_i <- rownames(spe)[i]
  for(j in (i+1):37){
    marker_j <- rownames(spe)[j]
    p <- dittoScatterPlot(ds_spe, marker_i, marker_j, assay.x = 'exprs_uncomp', assay.y = 'exprs_uncomp', opacity = 0.1) |
      dittoScatterPlot(ds_spe, marker_i, marker_j, assay.x = 'exprs', assay.y = 'exprs', opacity = 0.1)
    ggsave(file.path(DIRECTORIES$OUTPUT_SPILLOVER, 'Spillover_scatter', paste0(marker_i, '_', marker_j, '.png')), p)
  }
}

# Plot cells for markers
masks <- readRDS("data/processed/data_analysis/masks.rds")
dir.create(file.path(DIRECTORIES$OUTPUT_SPILLOVER, 'marker_plots'), showWarnings = F)
for(m in rownames(spe)){
  dir.create(file.path(DIRECTORIES$OUTPUT_SPILLOVER, 'marker_plots', m), showWarnings = F)
  
  plotCells(masks, spe, 'cell_id', 'img_id',
            colour_by = m, exprs_values = 'exprs',
            display = "single",
            image_title = list(position = "topright",
                               margin = c(5,5),
                               font = 1,
                               cex = 1),
            legend=list(colour_by.legend.cex=3),
            save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_SPILLOVER, 'marker_plots', 
                                                  m, paste0(m, '.png'))))
}

# Correlation
df <- as.data.frame(t(assay(spe, 'exprs')))
corr <- cor(df, method = 'pearson')
pdf(file.path(DIRECTORIES$OUTPUT_SPILLOVER, 'correlation.pdf'))
corrplot(corr, type = 'lower', diag = F)
dev.off()

# Save 
saveRDS(spe, file.path(DIRECTORIES$DATA_ANALYSIS, 'spe_comp.rds'))
