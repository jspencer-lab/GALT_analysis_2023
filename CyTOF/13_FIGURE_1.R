# Overview ----------------------------------------------------------------

# Cytof analysis
# TODO

library(flowCore)
library(CATALYST)
library(dplyr)
library(scater)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dittoSeq)
library(viridis)
library(tidyverse)
library(pheatmap)

source("config/globals.R")
source("functions/avg_exp.R")

DIRECTORIES$OUTPUT_FIGURE <- file.path(DIRECTORIES$OUTPUT, '14_FIGURE_2')
dir.create(DIRECTORIES$OUTPUT_FIGURE, showWarnings = F)

# Load --------------------------------------------------------------------
fs <- read.flowSet(path = "data/processed/data_cytof/experiment_43739_20230721111349791_files/")

sce <- prepData(fs)
colnames(sce) <- paste0('Cell_', 1:ncol(sce))

colData(sce)$individual <- ifelse(grepl('160 ', colData(sce)$sample_id), '160', 
                                  ifelse(grepl('161 ', colData(sce)$sample_id), '161', 
                                         ifelse(grepl('162 ', colData(sce)$sample_id), '162', 
                                                ifelse(grepl('166 ', colData(sce)$sample_id), '166', 
                                                       ifelse(grepl('165 ', colData(sce)$sample_id), '165', 'TBD')))))

colData(sce)$tissue <- ifelse(grepl('PBMC ', colData(sce)$sample_id), 'PBMC', 'GALT')
colData(sce)$cell_type <- sub(".*_([A-Za-z0-9 ]+)\\.fcs", "\\1", colData(sce)$sample_id, perl=TRUE)

# Attach visne
df <- as.data.frame(t(assay(sce, 'counts'))) %>% select(tSNE1, tSNE2)
rownames(df) <- colnames(sce)

reducedDim(sce, 'visne') <- df

colours <- c('firebrick2', 'dodgerblue2', 'orange', 'seagreen', 'magenta', 'purple3', 'grey', 'brown', 'cyan')
names(colours) <- unique(sce$cell_type)

# Visne plots -------------------------------------------------------------

col <- 'viridis'
plot_list <- multi_dittoDimPlot(sce, var = c('CD45RB', 'CD20', 'CD138', 'CD38', 'CD27', 'IgG',
                                             'IgD', 'IgA', 'CD21', 'CD10', 'CD24', 'IgM'), 
                                reduction.use = "visne",
                                assay = "exprs", size = 0.1, list.out = TRUE)
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis(option=col) + 
                      theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()))
p <- cowplot::plot_grid(plotlist = plot_list, ncol = 6)
ggsave(file.path(DIRECTORIES$OUTPUT_FIGURE, 'umap_features.png'), p, width=15, height=4)

p <- (dittoDimPlot(sce[,sce$tissue == 'GALT'], var = 'cell_type', size = 0.2, color.panel = colours) + ggtitle('GALT') + 
    theme(legend.position = 'none', panel.grid = element_blank(), 
          axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())) |
  (dittoDimPlot(sce[,sce$tissue == 'PBMC'], var = 'cell_type', size = 0.2, color.panel = colours) + ggtitle('PBMC') + 
     theme(panel.grid = element_blank(),
           axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()))
ggsave(file.path(DIRECTORIES$OUTPUT_FIGURE, 'umap_split.png'), p, width=7, height=4)

df <- as.data.frame(t(assay(sce, 'exprs')))
df$tSNE1 <- assay(sce, 'counts')['tSNE1',]
df$tSNE2 <- assay(sce, 'counts')['tSNE2',]
df$tissue <- sce$tissue
df$cell_type <- sce$cell_type

p <- ggplot(df[df$tissue == 'GALT',], 
            aes(x=tSNE1, y=tSNE2, color=cell_type)) + geom_point(size=0.5) + theme_minimal() + 
  scale_color_manual(values = colours) + 
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
        axis.title = element_blank())
p2 <- ggplot(df[df$tissue == 'PBMC',], 
            aes(x=tSNE1, y=tSNE2, color=cell_type)) + geom_point(size=0.5) + theme_minimal() + 
  scale_color_manual(values = colours) + 
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
        axis.title = element_blank())
ggsave(file.path(DIRECTORIES$OUTPUT_FIGURE, 'umap_split.png'), p | p2, width=8, height=4)

dir.create(file.path(DIRECTORIES$OUTPUT_FIGURE, 'umap_features'), showWarnings = F)
for(m in c('CD45RB', 'CD20', 'CD138', 'CD38', 'CD27', 'IgG',
           'IgD', 'IgA', 'CD21', 'CD10', 'CD24', 'IgM')){
  df$marker <- df[[m]]
  p <- ggplot(df, aes(x=tSNE1, y=tSNE2, color=marker)) + geom_point(size=0.5) + theme_minimal() + 
    scale_color_gradientn(colors=c('#46095D', '#25AC82', '#FFF001', '#F79615', '#EE2C2C')) + 
    theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
          axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
          axis.title = element_blank())
  
  ggsave(file.path(DIRECTORIES$OUTPUT_FIGURE, 'umap_features', paste0(m,'.png')), p, width=4, height=4)
}

# Heatmaps ------------------------------------------------------------------

df_GALT <- average_expression(sce[,sce$tissue=='GALT'], 'cell_type', genes=rownames(sce), scale=F, method = 'mean')
rownames(df_GALT) <- paste0(rownames(df_GALT), '_GALT')

df_PBMC <- average_expression(sce[,sce$tissue=='PBMC'], 'cell_type', genes=rownames(sce), scale=F, method = 'mean')
rownames(df_PBMC) <- paste0(rownames(df_PBMC), '_PBMC')

df <- t(scale(rbind(df_GALT, df_PBMC)))
    
markers <- c('HLA-DR', 'CD40', 'CD69', 'FcRL4', 'CD138', 'CCR10', 'CCR7', 'a4', 'b7', 'b1')
pdf(file.path(DIRECTORIES$OUTPUT_FIGURE, 'heatmap_combined.pdf'), width=3.93, height=2.75)
print(pheatmap(df[markers,],
         breaks = seq(-4, 4, length.out = 101), cluster_cols = F, cluster_rows = F,
         color=colorRampPalette(c('dodgerblue3', 'white', 'firebrick2'))(100)))
dev.off()

# Pie charts
df <- as.data.frame(colData(sce)) %>% group_by(tissue, cell_type) %>% summarise(count= n())

p <- ggplot(df, aes(x='', y=count, fill=cell_type)) + geom_bar(stat='identity') + coord_polar(theta='y', direction = -1) + 
  facet_wrap(~ tissue) + scale_fill_manual(values = colours) + theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank()
  )
ggsave(file.path(DIRECTORIES$OUTPUT_FIGURE, 'pie_chart_tissue.pdf'), p)

# Proportions
df <- as.data.frame(colData(sce)) %>% group_by(tissue, cell_type, individual) %>% summarise(count= n()) %>%
  filter(cell_type == 'DN2') %>% mutate(proportion = 100*count/9214) %>% ungroup() %>%
  select(individual, tissue, proportion) %>%
  pivot_wider(names_from = 'tissue', values_from = 'proportion')
write.csv(df, file.path(DIRECTORIES$OUTPUT_FIGURE, 'DN2_data_for_prism.csv'))

df <- as.data.frame(colData(sce)) %>% group_by(tissue, cell_type, individual) %>% summarise(count= n()) %>%
  filter(cell_type == 'DN1') %>% mutate(proportion = 100*count/9214) %>% ungroup() %>%
  select(individual, tissue, proportion) %>%
  pivot_wider(names_from = 'tissue', values_from = 'proportion')
write.csv(df, file.path(DIRECTORIES$OUTPUT_FIGURE, 'DN1_data_for_prism.csv'))

df <- as.data.frame(colData(sce)) %>% group_by(tissue, cell_type, individual) %>% summarise(count= n()) %>%
  filter(cell_type %in% c('DN1', 'DN2')) %>% pivot_wider(names_from = 'cell_type', values_from = 'count') %>% 
  mutate(DNratio = DN2/DN1) %>% select(individual, tissue, DNratio) %>%
  pivot_wider(names_from = 'tissue', values_from = 'DNratio')
write.csv(df, file.path(DIRECTORIES$OUTPUT_FIGURE, 'DNratio_data_for_prism.csv'))
