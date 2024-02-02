library(cytomapper)
library(stringr)
library(imcRtools)
library(dplyr)
library(scoper)
library(scater)
library(batchelor)
library(Rphenograph)
library(igraph)
library(dittoSeq)
library(tidyr)
library(tibble)
library(viridis)
library(pheatmap)

source("config/globals.R")

markers <- c('CD11b', 'CD163', 
             'CD68', 'FITC-ITGAX',
             'CD103', 
             'FCRL4',
             'GzB', 'DNASE1L3')

DIRECTORIES$OUTPUT_CD68ITGAX <- file.path("output/CD68_ITGAX_v2")
dir.create(DIRECTORIES$OUTPUT_CD68ITGAX, showWarnings = F)

# Load ----------------------------------------------------------
spe <- readRDS("data/processed/analysis/CD69_ITGAX_segmented/spe_load_CD68_ITGAX.rds")
images <- readRDS("data/processed/analysis/images.rds")
masks <- readRDS("data/processed/analysis/CD69_ITGAX_segmented/masks_CD68_ITGAX.rds")

# Analyse -------------------------------------------------------
exp <- as.data.frame(t(assay(spe, 'exprs'))) %>% 
  select(-DNA1, -DNA2) 
exp_pl <- pivot_longer(exp, colnames(exp), names_to = 'marker', values_to = 'intensity')

p <- ggplot(exp_pl, aes(x=marker, y=intensity)) + geom_boxplot(outlier.shape=NA, fill='lightgrey') + theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) + ggtitle("Cells segmented on CD68+ and/or ITGAX+")
ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, "boxplot_all.pdf"), p)
p <- ggplot(exp_pl[exp_pl$marker %in% markers,], aes(x=marker, y=intensity)) + 
  geom_boxplot(outlier.shape=NA, fill='lightgrey') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) + ggtitle("Cells segmented on CD68+ and/or ITGAX+")
ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, "boxplot_MacDC_markers.pdf"), p)


# GC - DNASE1L3 ------------------------------------------------
dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'GC_cells'), showWarnings = F)

exp <- as.data.frame(t(assay(spe, 'exprs'))) %>% select(-DNA1, -DNA2) %>% 
  mutate(img_id = spe$img_id, is_GC = spe$is_GC)
exp_pl <- pivot_longer(exp, colnames(exp)[1:33], names_to = 'marker', values_to = 'intensity') %>% 
  group_by(img_id, is_GC, marker) %>% summarise(avg_int = mean(intensity))

DNASE1L3 <- exp_pl[exp_pl$marker == 'DNASE1L3',] %>% 
  pivot_wider(names_from = 'is_GC', values_from = 'avg_int')
write.csv(DNASE1L3, file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'GC_cells', 'DNASE1L3_GC.csv'))

exp_pl <- pivot_longer(exp, colnames(exp)[1:33], names_to = 'marker', values_to = 'intensity')
ggplot(exp_pl, aes(x=marker, y=intensity, fill=is_GC)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Batch correction ---------------------------------------------
set.seed(101)
out <- fastMNN(spe, batch = spe$mcd,
               auto.merge = TRUE,
               subset.row = markers,
               assay.type = "exprs")
reducedDim(spe, "fastMNN") <- reducedDim(out, "corrected")

p <- plotReducedDim(spe, 'fastMNN', colour_by = 'mcd', by_exprs_values = 'exprs')
ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'umap_corrected.png'), p)

dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'umap_expression'), showWarnings = F)
for(m in PANEL$Actual){
  p <- plotReducedDim(spe, 'fastMNN', colour_by = m, by_exprs_values = 'exprs')  
  ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, "umap_expression", paste0(m, ".png")), p)
}

# Cluster ------------------------------------------------------
set.seed(101)
mat <- reducedDim(spe, "fastMNN")
out <- Rphenograph(mat, k = 150)
spe$cluster_corrected <- factor(membership(out[[2]]))
spe$cluster_corrected

spe@metadata$cluster_colours <- c('firebrick2', 'dodgerblue3', 'cyan', 
                                  'darkgoldenrod', 'orange', 'springgreen4',
                                  'yellow', 'pink', 'olivedrab1',
                                  'beige', 'purple')
names(spe@metadata$cluster_colours) <- levels(spe$cluster_corrected)

p <- plotReducedDim(spe, 'fastMNN', colour_by = 'cluster_corrected', text_by = 'cluster_corrected') +
  scale_color_manual(values = spe@metadata$cluster_colours)
ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX , "umap_clusters.pdf"), p)

p <- dittoHeatmap(spe, 
             genes = markers,
             assay = "exprs", scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("cluster_corrected"),
             annot.colors = spe@metadata$cluster_colours)
ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX , "heatmap_cells_markers.pdf"), p)

p <- dittoHeatmap(spe, 
                  genes = PANEL$Actual[1:33],
                  assay = "exprs", scale = "none",
                  heatmap.colors = viridis(100), 
                  annot.by = c("cluster_corrected"),
                  annot.colors = spe@metadata$cluster_colours)
ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX , "heatmap_cells_all.pdf"), p)

exp <- as.data.frame(t(assay(spe, 'exprs'))) %>% select(-DNA1, -DNA2) %>% mutate(cluster = spe$cluster_corrected)
exp_pl <- pivot_longer(exp, colnames(exp)[1:33], names_to = 'marker', values_to = 'intensity') %>% 
  group_by(cluster, marker) %>% summarise(avg_int = mean(intensity))
exp <- pivot_wider(exp_pl, names_from = 'marker', values_from = 'avg_int') %>% column_to_rownames('cluster')

rowannot <- data.frame(cluster = 1:length(levels(spe$cluster_corrected)))
rownames(rowannot) <- rownames(exp)

rowcols <- spe@metadata$cluster_colours
names(rowcols) <- rownames(exp)

p <- pheatmap(exp[,markers], color = colorRampPalette(c('white', 'red'))(100),
              annotation_row = rowannot, 
              annotation_colors = list(cluster = rowcols))
ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'heatmap_clusters_markers.pdf'), p)

p <- pheatmap(exp, color = colorRampPalette(c('white', 'red'))(100),
              annotation_row = rowannot, 
              annotation_colors = list(cluster = rowcols))
ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'heatmap_clusters_all.pdf'), p)


# Cluster overlays -----------------------------------------------
dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'cluster_overlays'), showWarnings = F)

plotCells(masks, spe, 'cell_id', 'img_id',
          colour_by = 'cluster_corrected', display = 'single',
          colour = list(cluster_corrected = spe@metadata$cluster_colours),
          image_title = list(position = "topright",
                             margin = c(5,5),
                             font = 1,
                             cex = 1),
          save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'cluster_overlays', 'clusters.png'))
)

col <- list(cluster_corrected = c('red'))
dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'cluster_individual'), showWarnings = F)

for(i in unique(spe$cluster_corrected)){
  dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'cluster_individual', i), showWarnings = F)
  names(col$cluster_corrected) <- as.character(i)
  plotCells(masks, spe[ ,spe$cluster_corrected == i], 'cell_id', 'img_id',
            colour_by = 'cluster_corrected',
            colour = col, display = 'single',
            image_title = list(position = "topright",
                               margin = c(5,5),
                               font = 1,
                               cex = 1),
            save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'cluster_individual', i, 
                                                  paste0('cluster_', i, '.png')))
  )
}

# DNASE1L3 positive -------------------------------------------------------------------------
dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'DNASE1L3_pos'), showWarnings = F)
spe$DNASE1L3pos <- spe$cluster_corrected %in% c(5,7,8,2,6,11)

plotCells(masks, spe, 'cell_id', 'img_id',
          display = 'single',
          image_title = list(position = "topright",
                             margin = c(5,5),
                             font = 1,
                             cex = 1),
          colour_by = 'DNASE1L3pos', colour = list(DNASE1L3pos = c(`TRUE` = 'red', `FALSE` = 'lightgrey')),
          save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'DNASE1L3_pos', 
                                                paste0('DNASE1L3pos.png')))
)

df <- as.data.frame(colData(spe)) %>% group_by(img_id, DNASE1L3pos) %>% 
  summarise(avg_dist = mean(distance_to_epithelium)) %>% 
  pivot_wider(names_from = 'DNASE1L3pos', values_from = 'avg_dist')
write.csv(df, file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'DNASE1L3_pos', 
                        paste0('DNASE1L3pos_distance_epithelium.csv')))

df <- as.data.frame(colData(spe)) %>% group_by(img_id, DNASE1L3pos) %>% 
  summarise(avg_dist = mean(distance_to_epithelium_scaled)) %>% 
  pivot_wider(names_from = 'DNASE1L3pos', values_from = 'avg_dist')
write.csv(df, file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'DNASE1L3_pos', 
                        paste0('DNASE1L3pos_distance_epithelium_scaled.csv')))


# Markers between pos and neg
exp <- as.data.frame(t(assay(spe, 'exprs'))) %>% 
  select(-DNA1, -DNA2) %>% mutate(DNASE1L3pos = spe$DNASE1L3pos) %>% 
  pivot_longer(PANEL$Actual[1:33], names_to = 'marker', values_to = 'intensity')


p <- ggplot(exp, aes(x=marker, y=intensity, fill=DNASE1L3pos)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'DNASE1L3_pos', 'markers_DNASE1L3_pos_neg.pdf') ,p)


# SED ----------------------------------------------------------------
# pos_only <- spe[,spe$DNASE1L3pos]
# 
# pos_only$is_sed <- FALSE
# 
# # Process each image. Manually alter number below
# img_no <- 14
# 
# plotPixels(images[img_no], pos_only, masks[img_no], 'cell_id', 'img_id',
#            colour_by = c('E-Cadherin', 'CD20', 'CD3'), 
#            bcg = list(`E-Cadherin` = c(0,3,1),
#                       `CD20` = c(0,10,1),
#                       `CD3` = c(0,4,1)),
#            colour = list(`E-Cadherin` = c('black', 'yellow'),
#                          `CD20` = c('black', 'cyan'),
#                          `CD3` = c('black', 'magenta')),
#            scale_bar = NULL, legend=NULL, image_title = list(cex = 0.75))
# 
# selectedPoints <- c()
# 
# df <- as.data.frame(spatialCoords(pos_only[,pos_only$img_id == names(images)[img_no]]))
# rownames(df) <- colnames(pos_only[,pos_only$img_id == names(images)[img_no]])
# df$Location_Center_Y <- -1*df$Location_Center_Y
# X11()
# plot(df, pch = 16, col = "red")
# selectedPoints <- c(selectedPoints, gatepoints::fhs(df, mark = TRUE))
# 
# write.csv(data.frame(ids = selectedPoints),
#           file.path(DIRECTORIES$DATA_PROCESSED, 'analysis', 'CD69_ITGAX_segmented', paste0(img_no, '_SED_selection.csv')))
# # pos_only$is_sed[colnames(pos_only) %in% selectedPoints] <- TRUE
# 
# 
# # Load SED
# pos_only <- spe[,spe$DNASE1L3pos]
# 
# sed_cells <- bind_rows(lapply(list.files("data/processed/analysis/CD69_ITGAX_segmented/SED", full.names = T), 
#                               FUN = read.csv))$ids
# 
# pos_only$is_sed <- colnames(pos_only) %in% sed_cells
# 
# plotCells(masks[10], pos_only, 'cell_id', 'img_id',
#           colour_by = 'is_sed', colour = list(is_sed = c(`TRUE` = 'red', `FALSE` = 'dodgerblue3')))
# 
# df <- as.data.frame(t(assay(pos_only, 'exprs'))) %>% mutate(img_id = pos_only$img_id,
#                                                             is_sed = pos_only$is_sed) %>%
#   group_by(img_id, is_sed) %>% summarise(avg_lysozyme = mean(Lysozyme)) %>%
#   pivot_wider(names_from = 'is_sed', values_from = 'avg_lysozyme')
# write.csv(df, file.path(DIRECTORIES$DATA_PROCESSED, 'analysis', 'CD69_ITGAX_segmented', 'lysozyme.csv'))
# 
# df <- as.data.frame(t(assay(pos_only, 'exprs'))) %>% mutate(img_id = pos_only$img_id,
#                                                             is_sed = pos_only$is_sed) %>%
#   group_by(img_id, is_sed) %>% summarise(avg_NOX2 = mean(NOX2)) %>%
#   pivot_wider(names_from = 'is_sed', values_from = 'avg_NOX2')
# write.csv(df, file.path(DIRECTORIES$DATA_PROCESSED, 'analysis', 'CD69_ITGAX_segmented', 'NOX2.csv'))


# Save --------------------------------------------------------------------------------------
saveRDS(spe, "data/processed/analysis/spe_clustered.rds")
