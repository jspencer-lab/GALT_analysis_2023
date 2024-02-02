library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(viridis)

source("config/globals.R")
source("function/SpatialDimPlotMJP.R")

DIRECTORIES$OUTPUT_LOAD <- file.path(DIRECTORIES$OUTPUT, '1_load_seurat')
dir.create(DIRECTORIES$OUTPUT_LOAD, showWarnings = F)

# Load ----
spatial_seurat <- list()
for (run_id in IDS$RUN){
  for(sample_id in IDS[[run_id]]){
    valid <- SAMPLE_METADATA$Valid[SAMPLE_METADATA$Run == run_id & SAMPLE_METADATA$ID == sample_id]
    if(valid){
      id <- paste0(run_id, '_', sample_id)

      sp <- Load10X_Spatial(file.path(DIRECTORIES$DATA_RAW, run_id, sample_id, 'outs'), )
      sp <- RenameCells(sp, add.cell.id = id)
      sp$orig.ident <- id

      # Manual alignment variables - used to line up the histology with the data co-ordinates for more
      # control over plots (SpatialDimPlot etc are limited)
      sp@misc$x_alignment_min <- SAMPLE_METADATA$Alignment_x_min[SAMPLE_METADATA$Run == run_id & 
                                                                   SAMPLE_METADATA$ID == sample_id]
      sp@misc$x_alignment_max <- SAMPLE_METADATA$Alignment_x_max[SAMPLE_METADATA$Run == run_id & 
                                                                   SAMPLE_METADATA$ID == sample_id]
      sp@misc$y_alignment_min <- SAMPLE_METADATA$Alignment_y_min[SAMPLE_METADATA$Run == run_id & 
                                                                   SAMPLE_METADATA$ID == sample_id]
      sp@misc$y_alignment_max <- SAMPLE_METADATA$Alignment_y_max[SAMPLE_METADATA$Run == run_id & 
                                                                   SAMPLE_METADATA$ID == sample_id]
      
      # Attach image
      sp@misc$image <- png::readPNG(file.path("data/raw", run_id, sample_id, 
                                              "outs/spatial/tissue_hires_image.png"))
      
      # Save to list
      spatial_seurat[[id]] <- sp
      rm(sp)
    }
  }
}

# QC ----
for(id in names(spatial_seurat)){
  p <- SpatialFeaturePlot(spatial_seurat[[id]], c('nCount_Spatial', 'nFeature_Spatial'))
  ggsave(file.path(DIRECTORIES$OUTPUT_LOAD, paste0(id, '_metrics.pdf')), p, height=5)
}

# Plots ----
dir.create(file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots'), showWarnings = F)

normalised <- list()
for(id in names(spatial_seurat)){
  normalised[[id]] <- SCTransform(spatial_seurat[[id]], assay = 'Spatial')  
}

dir.create(file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots', 'blanks'), showWarnings = F)
dir.create(file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots', 'alpha'), showWarnings = F)
dir.create(file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots', 'noalpha'), showWarnings = F)
dir.create(file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots', 'CROPS_alpha'), showWarnings = F)

crop_xmin <- 1650
crop_xmax <- 2500
crop_ymin <- 2000
crop_ymax <- 2700
for(id in names(normalised)){
  # Blanks
  p <- ggplot(data.frame(imagecol=c(0), imagerow=c(0)), aes(x=imagecol, y=imagerow)) + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
    ggpubr::background_image(normalised[[id]]@misc$image) +
    xlim(normalised[[id]]@misc$x_alignment_min, normalised[[id]]@misc$x_alignment_max) + 
    ylim(normalised[[id]]@misc$y_alignment_max, normalised[[id]]@misc$y_alignment_min)
  filename <- file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots', 'blanks',
                        paste0(id, '_H&E.png'))
  ggsave(filename, p, width=10, height=10)
  
  if(id == 'APPGUT_A1'){
    p <- png::readPNG(filename)
    q <- imagefx::crop.image(p, crop_ymin, crop_xmin, crop_ymax, crop_xmax)
    png::writePNG(q$img.crop, file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots', 'blanks',
                                        paste0(id, '_H&E_CROP.png')))
  }
  
  for(gene in c('MS4A1', 'IGHM', 'IGHD', 'MKI67', 'CD3E', 'ITGAX')){
    p <- SpatialPlotMJP(normalised[[id]], gene, alpha_fill = TRUE, pt_size=3) +
      theme(legend.position = 'none') + scale_color_viridis()

    filename <- file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots', 'alpha',
                          paste0(id, '_', gene, '.png'))
    ggsave(filename, p, width=10, height=10)
    
    if(id == 'APPGUT_A1'){
      p <- png::readPNG(filename)
      q <- imagefx::crop.image(p, crop_ymin, crop_xmin, crop_ymax, crop_xmax)
      png::writePNG(q$img.crop, file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots', 'CROPS_alpha',
                                          paste0(id, '_', gene, '_CROP.png')))
    }

    p <- SpatialPlotMJP(normalised[[id]], gene, alpha_fill = FALSE, pt_size = 3) +
      theme(legend.position = 'none') + scale_color_viridis()

    filename <- file.path(DIRECTORIES$OUTPUT_LOAD, 'spatial_feature_plots', 'noalpha',
                          paste0(id, '_', gene, '_noalpha.png'))
    ggsave(filename, p, width=10, height=10)
  }
}

# Save ----
saveRDS(spatial_seurat, file.path(DIRECTORIES$DATA_PROCESSED, '1_load_seurat.rds'))
