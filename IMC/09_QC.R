# Overview ----------------------------------------------------------------

# 9_QC
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

source("config/globals.R")

DIRECTORIES$OUTPUT_QC <- file.path(DIRECTORIES$OUTPUT, '09_QC')
dir.create(DIRECTORIES$OUTPUT_QC, showWarnings = F)

# Load ------------------------------------------------------------------------
spe <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'spe_comp.rds'))
images <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'images.rds'))
masks <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'masks.rds'))

# Segmentation ----------------------------------------------------------------

cur_images <- images
cur_images <- normalize(cur_images, separateImages = TRUE)
cur_images <- normalize(cur_images, inputRange = c(0, 0.2))

dir.create(file.path(DIRECTORIES$OUTPUT_QC, 'segmentation'), showWarnings = F)
plotPixels(cur_images,
           mask = masks,
           display = 'single', 
           save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_QC, 'segmentation', 'segmentation.png')),
           img_id = "img_id",
           missing_colour = "white",
           colour_by = c("DNA1"),
           colour = list(DNA1 = c("black", "blue")),
           image_title = NULL,
           legend = list(colour_by.title.cex = 0.7,
                         colour_by.labels.cex = 0.7))

cur_cells <- sample(seq_len(ncol(spe)), 2000)

dittoHeatmap(spe[,cur_cells], genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", cluster_cols = TRUE, scale = "none",
             heatmap.colors = viridis(100), annot.by = "Metadata_acname")


q <- spe
cur_cells <- sample(seq_len(ncol(q)), 2000)
assay(q, 'new') <- asinh(assay(q, 'counts')/0.1)

dittoHeatmap(q[,cur_cells], genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "new", cluster_cols = TRUE, scale = "none",
             heatmap.colors = viridis(100), annot.by = "Metadata_acname")



# Image level QC ---------------------------------------------------------------

q <- scaleImages(images, 2**16)

cur_snr_img <- lapply(q, function(img){
  mat <- apply(img, 3, function(ch){
    # Otsu threshold
    thres <- otsu(ch, range = c(min(ch), max(ch)))
    # Signal-to-noise ratio (correction for divide by zero)
    snr <- (mean(ch[ch > thres])+0.0000001) / (mean(ch[ch <= thres])+0.0000001)
    # Signal intensity
    ps <- mean(ch[ch > thres])
    
    return(c(snr = snr, ps = ps))
  })
  t(mat) %>% as.data.frame() %>% 
    mutate(marker = colnames(mat)) %>% 
    pivot_longer(cols = c(snr, ps))
})

cur_snr_img <- do.call(rbind, cur_snr_img)

cur_snr_img <- cur_snr_img %>% 
  group_by(marker, name) %>%
  summarize(mean = mean(value),
            ci = qnorm(0.975)*sd(value)/sqrt(n())) %>%
  pivot_wider(names_from = name, values_from = c(mean, ci)) %>%
  mutate(mean_ps_log2 = log2(mean_ps),
         mean_snr_log2 = log2(mean_snr))

p <- ggplot(cur_snr_img) +
  geom_point(aes(log2(mean_ps), log2(mean_snr))) +
  geom_label_repel(aes(log2(mean_ps), log2(mean_snr), label = marker)) +
  theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
  xlab("Signal intensity [log2]")
ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'SNR_image.pdf'), p, width=10, height=10)

write.csv(cur_snr_img, file.path(DIRECTORIES$OUTPUT_QC, 'snr_data_image.csv'))

# Area covered by cells
pc_covered <- list()
for(i in names(images)){
  dim_img <- dim(imageData(images[i][[names(images[i])]]))
  size <- dim_img[1] * dim_img[2]
  
  covered <- sum(spe$AreaShape_Area[spe$img_id == i])
  
  pc_covered[[i]] <- covered/size
}

df <- data.frame(img = names(pc_covered), pc_covered = as.numeric(pc_covered))
p <- ggplot(df, aes(x=img, y=pc_covered)) + geom_point() + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(0,1)
ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'area_covered.pdf'), p)

# Mean expression per image
image_mean <- aggregateAcrossCells(spe, 
                                   ids = spe$img_id, 
                                   statistics="mean",
                                   use.assay.type = "counts")
assay(image_mean, "exprs") <- asinh(counts(image_mean))

p <- dittoHeatmap(image_mean, genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", cluster_cols = TRUE, scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("img_id", "tissue", "Metadata_acname"),
             show_colnames = TRUE)
ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'mean_int_per_image.pdf'), p, width=10, height = 10)


# Cell SNR
set.seed(101)
mat <- apply(assay(spe, "exprs"), 1, function(x){
  cur_model <- Mclust(x, G = 2)
  mean1 <- mean(x[cur_model$classification == 1])
  mean2 <- mean(x[cur_model$classification == 2])
  
  signal <- ifelse(mean1 > mean2, mean1, mean2)
  noise <- ifelse(mean1 > mean2, mean2, mean1)
  
  return(c(snr = signal/noise, ps = signal))
})

cur_snr <- t(mat) %>% as.data.frame() %>% 
  mutate(marker = colnames(mat))

p <- cur_snr %>% ggplot() +
  geom_point(aes(log2(ps), log2(snr))) +
  geom_label_repel(aes(log2(ps), log2(snr), label = marker)) +
  theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
  xlab("Signal intensity [log2]")
ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'SNR_cell.pdf'), p, width=10, height=10)

write.csv(cur_snr, file.path(DIRECTORIES$OUTPUT_QC, 'snr_data_cell.csv'))

p <- colData(spe) %>%
  as.data.frame() %>%
  group_by(img_id) %>%
  ggplot() +
  geom_boxplot(aes(img_id, AreaShape_Area)) +
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cell area") + xlab("")
ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'cell_size_per_image.pdf'), p)


# Actual QC -----------------------------------------------------------------------------
dir.create(file.path(DIRECTORIES$OUTPUT_QC, 'pass_qc'), showWarnings = F)

df <- data.frame(img = spe$img_id,
                 size = spe$AreaShape_Area, 
                 int = colSums(assay(spe, 'exprs')[rowData(spe)$use_channel,]),
                 DNA = colMeans(assay(spe, 'exprs')[!rowData(spe)$use_channel,]))

ggplot(df, aes(x=int)) + geom_histogram(bins=100)
ggplot(df, aes(x=DNA)) + geom_histogram(bins=60)
ggplot(df, aes(x=size)) + geom_histogram(bins=60) + scale_x_log10()

ggplot(df, aes(x=int, y=size)) + geom_point()
ggplot(df, aes(x=DNA, y=size, color=img)) + geom_point() + geom_vline(xintercept = 4.5) + geom_hline(yintercept = 15)
ggplot(df, aes(x=DNA, y=int)) + geom_point()

# QC threshold
# df$pass_qc <- df$size >= 10
df$pass_qc <- df$DNA > 4.5 & df$size > 10

spe$pass_qc <- df$pass_qc

# QC plots
p <- ggplot(df, aes(x=DNA, y=size, color=pass_qc)) + geom_point() + scale_color_manual(values=c('red', 'lightgrey')) +
  theme_minimal()
p2 <- ggplot(df %>% group_by(pass_qc) %>% count(), aes(y=n, x='', fill=pass_qc)) + geom_bar(stat='identity') + 
  coord_polar(theta = 'y', direction = -1) + theme_minimal() + scale_fill_manual(values=c('red', 'lightgrey'))
p3 <- ggplot(df %>% group_by(img, pass_qc) %>% count(), aes(y=n, x=img, fill=pass_qc)) + 
  geom_bar(stat='identity', position = 'fill') + scale_fill_manual(values=c('red', 'lightgrey'))  +
  theme_minimal()
ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'pass_qc', 'pass_qc.png'), p | (p2/p3))
  
# plotCells(masks, spe, 'cell_id', 'img_id', colour_by = 'pass_qc',
#           colour = list('pass_qc' = c(`TRUE` = 'lightgrey', `FALSE` = 'red')),
#           display = 'single',
#           image_title = list(position = "topright",
#                              margin = c(5,5),
#                              font = 1,
#                              cex = 1),
#           save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_QC, 'pass_qc', 'pass_qc.png')))

plotPixels(cur_images, spe[,!spe$pass_qc], masks, 'cell_id', 'img_id',
           colour_by = 'DNA1', outline_by = 'pass_qc',
           colour = list('pass_qc' = c(`TRUE` = 'lightgrey', `FALSE` = 'red')),
           display = 'single',
           image_title = list(position = "topright",
                              margin = c(5,5),
                              font = 1,
                              cex = 1),
           save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_QC, 'pass_qc', 'pass_qc.png')))

# Remove cells
spe <- spe[,spe$pass_qc]

# Dimensionality reduction ---------------------------------------------------------

dir.create(file.path(DIRECTORIES$OUTPUT_QC, 'umap'), showWarnings = F)

set.seed(101)
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 

p <- plotReducedDim(spe, 'UMAP', colour_by = 'img_id') + theme_classic()
ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'umap', 'umap_img_id.png'), p)

p <- plotReducedDim(spe, 'UMAP', colour_by = 'Metadata_acname') + theme_classic()+ 
  scale_color_manual(values = COLOURS$MCD)
ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'umap', 'umap_mcd.png'), p)

p <- plotReducedDim(spe, 'UMAP', colour_by = 'tissue') + theme_classic() + 
  scale_color_manual(values = COLOURS$TISSUE)
ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'umap', 'umap_tissue.png'), p)


dir.create(file.path(DIRECTORIES$OUTPUT_QC, 'umap', 'markers'), showWarnings = F)
for(m in PANEL$Actual){
  p <- plotReducedDim(spe, 'UMAP', colour_by = m, by_exprs_values = 'exprs') + theme_classic()
  ggsave(file.path(DIRECTORIES$OUTPUT_QC, 'umap', 'markers', paste0('umap_', m, '.png')), p)
}

# Save -----------------------------------------------------------------------------
saveRDS(spe, file.path(DIRECTORIES$DATA_ANALYSIS, 'spe_QC.rds'))

