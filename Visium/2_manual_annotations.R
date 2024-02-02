library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(viridis)

source("config/globals.R")
source("function/SpatialDimPlotMJP.R")

gene_annotations <- read.table("data/raw/APPGUT/A1/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                               sep = '\t') %>% 
  select('V1', 'V2') %>% rename('gene_id' = 'V1', 'gene_name' = 'V2') %>% 
  left_join(read.csv("../GeneAnnotations/gene_annotations_2022-10-31.csv") %>% 
              select('gene_id', 'description'), by='gene_id')

DIRECTORIES$OUTPUT_MANANN <- file.path(DIRECTORIES$OUTPUT, '2_manual_annotations')
dir.create(DIRECTORIES$OUTPUT_MANANN, showWarnings = F)

# Load ----
spatial_seurat <- readRDS(file.path(DIRECTORIES$DATA_PROCESSED, '1_load_seurat.rds'))

# Attach classifications & filter ----
filtered_seurat <- list()

for (id in names(spatial_seurat)){
  cl <- read.csv(file.path("data/processed", paste0(id, "_MANUAL.csv"))) %>% 
    mutate(ID = paste0(id, '_', Barcode)) %>%
    tibble::column_to_rownames('ID') %>% select(Manual)
  
  filtered_seurat[[id]] <- subset(spatial_seurat[[id]], cells = rownames(cl)) %>% AddMetaData(cl)
  
  Idents(filtered_seurat[[id]]) <- filtered_seurat[[id]]$Manual
  
  filtered_seurat[[id]] <- SCTransform(filtered_seurat[[id]], assay = 'Spatial', return.only.var.genes = F)
}

# Plots ----
crop_xmin <- 1850
crop_xmax <- 2400
crop_ymin <- 2400
crop_ymax <- 2700

dir.create(file.path(DIRECTORIES$OUTPUT_MANANN, 'annotations'), showWarnings = F)
for(id in names(filtered_seurat)){
  df <- as.data.frame(filtered_seurat[[id]]@images$slice1@coordinates)
  df$classification <- Idents(filtered_seurat[[id]])
  
  p <- SpatialPlotMJP(filtered_seurat[[id]], fill='Manual', pt_size = 3) + 
    scale_color_manual(values=COLOURS$REGION) +
    theme(legend.position = 'none')
  filename <- file.path(DIRECTORIES$OUTPUT_MANANN, 'annotations', paste0(id, '_manual_ann_fill.png'))
  ggsave(filename, width = 10, height = 10)
  
  p <- SpatialPlotMJP(filtered_seurat[[id]], border='Manual', pt_size = 2) + 
    scale_color_manual(values=COLOURS$REGION) +
    theme(legend.position = 'none')
  filename <- file.path(DIRECTORIES$OUTPUT_MANANN, 'annotations', paste0(id, '_manual_ann_border.png'))
  ggsave(filename, width = 10, height = 10)
  
  if(id == 'APPGUT_A1'){
    p <- png::readPNG(filename)
    q <- imagefx::crop.image(p, crop_ymin, crop_xmin, crop_ymax, crop_xmax)
    png::writePNG(q$img.crop, file.path(DIRECTORIES$OUTPUT_MANANN, 'annotations', 
                                        paste0(id, '_manual_ann_CROP.png')))
  }
}

# Merge ----
merged <- merge(filtered_seurat$APPGUT_A1, filtered_seurat$APPGUT_D1)

# Markers ----
merged <- PrepSCTFindMarkers(merged)
markers <- FindMarkers(merged, assay = 'SCT', ident.1 = 'SED')

markers <- markers %>% mutate(gene_name = rownames(markers)) %>% 
  left_join(gene_annotations, by='gene_name')

write.csv(markers, file.path(DIRECTORIES$OUTPUT_MANANN, 'markers_SED_vs_follicle.csv'))

# TODO - markers heatmap
markers_plot <- c('DNASE1L3', 'IL22RA2', 'CXCL14', 'CCL23', 'CCL20', 'ITGA3', 'ITGAX',
                  'CCL15', 'IL1B', 'KLF4', 'LYZ', 'TNFRSF11A', 'C1R', 'C3', 'C1QA', 'ITGB1',
                  'BCL6', 'SLBP', 'PFN1', 'TCL1A', 'SIVA1', 'SERPINA9', 'RGS13', 'CR2',
                  'CDC20')
levels(merged) <- c('SED', 'Follicle')
p <- DoHeatmap(merged, assay = 'SCT', features = markers_plot, raster = FALSE, group.bar = F, 
               group.colors = c('black', 'white')) + theme(legend.position = 'bottom')
ggsave(file.path(DIRECTORIES$OUTPUT_MANANN, 'markers_heatmap.pdf'), p, width=4.5, height=4.5)

# markers_plot <- c(rownames(markers)[markers$avg_log2FC > 0][1:15],
#                   rownames(markers)[markers$avg_log2FC < 0][1:15])
# p <- DoHeatmap(merged, assay = 'SCT', features = markers_plot)
# ggsave(file.path(DIRECTORIES$OUTPUT_MANANN, 'markers_heatmap.pdf'), p)

# Plots spatial ----
dir.create(file.path(DIRECTORIES$OUTPUT_MANANN, 'markers_spatial'), showWarnings = F)
dir.create(file.path(DIRECTORIES$OUTPUT_MANANN, 'markers_spatial', 'crops'), showWarnings = F)

for(id in names(filtered_seurat)){
  for(gene in c('DNASE1L3', 'C1QA', 'C1QB', 'C1QC', 'C1R', 'C3')){
    p <- SpatialPlotMJP_v2(filtered_seurat[[id]], gene) + theme(legend.position = 'none')
    filename <- file.path(DIRECTORIES$OUTPUT_MANANN, 'markers_spatial', 
                          paste0(id, '_region_', gene, '.png'))
    ggsave(filename, width=10, height=10)  
    
    if(id == 'APPGUT_A1'){
      p <- png::readPNG(filename)
      q <- imagefx::crop.image(p, crop_ymin, crop_xmin, crop_ymax, crop_xmax)
      png::writePNG(q$img.crop, file.path(DIRECTORIES$OUTPUT_MANANN, 'markers_spatial', 'crops',
                                          paste0(id, '_region_', gene, '_CROP.png')))
    }
  }
}

saveRDS(filtered_seurat, "data/processed/2_manann_seurat.rds")


filtered_seurat <- readRDS("data/processed/2_manann_seurat.rds")
