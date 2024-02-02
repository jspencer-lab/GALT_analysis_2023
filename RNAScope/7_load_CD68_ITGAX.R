library(cytomapper)
library(stringr)
library(imcRtools)
library(dplyr)
library(scoper)
library(scater)
library(CATALYST)
library(dittoSeq)
library(tidyr)

source("config/globals.R")

IMAGE_METADATA <- read.csv("data/processed/extract/CD68_ITGAX_cells//Image.csv")
IMAGE_METADATA$mcd_file <- stringr::str_replace(IMAGE_METADATA$Metadata_ID, '_s0_a[0-9]', '')

DIRECTORIES$OUTPUT_CD68ITGAX <- file.path("output/CD68_ITGAX")
dir.create(DIRECTORIES$OUTPUT_CD68ITGAX, showWarnings = F)

# Load masks -----------
masks <- loadImages(file.path("data/processed/extract/CD68_ITGAX_cells/masks"), pattern='_CD68_ITGAX_cells.tiff')
mcols(masks)$img_id <- names(masks) <- str_match(names(masks), "(.*?)\\s*_ac_full_CD68_ITGAX_cells")[,2]
masks <- scaleImages(masks, 2^16-1)
# Reorder
masks <- masks[IMAGE_METADATA$Metadata_ID]

# Create spe ----------------
cell_data <- read.csv("data/processed/extract/CD68_ITGAX_cells/CD68_ITGAX_cells.csv") %>% 
  left_join(IMAGE_METADATA %>% select(ImageNumber, Metadata_ID, mcd_file), by='ImageNumber') %>% 
  rename(img_id = 'Metadata_ID', cell_id = 'ObjectNumber') %>% 
  left_join(SAMPLE_METADATA %>% select(ID, tissue, patient_id, block_id, slide), by = c('img_id' = 'ID'))

int_data <- t(as.matrix(cell_data[, paste0('Intensity_MeanIntensity_Stack_', paste0('c', 1:nrow(PANEL)))]))*(2**16)
rownames(int_data) <- PANEL$Actual
colnames(int_data) <- cell_data$cell_id
spatial_data <- DataFrame(cell_data[,c("Location_Center_X", "Location_Center_Y")])

# Create SpatialExperiment object
spe <- SpatialExperiment(list(counts=int_data,
                              exprs=asinh(int_data/1)),
                         spatialData = spatial_data,
                         spatialCoordsNames = c("Location_Center_X", "Location_Center_Y"),)
colData(spe)$img_id <- cell_data$img_id
colData(spe)$cell_id <- cell_data$cell_id
colData(spe)$tissue <- cell_data$tissue
colData(spe)$mcd <- cell_data$mcd_file

colnames(spe) <- paste0(spe$img_id, '_', spe$cell_id)

rowData(spe)$Metal_tag <- PANEL$Metal_tag

# Spillover correct -----------------------------------------------
sm <- readRDS("reference/spillover_matrix.rds")

rowData(spe)$channel_name <- paste0(PANEL$Metal_tag, 'Di')

spe <- compCytof(spe, sm, 
                 transform = TRUE, cofactor = 1,
                 isotope_list = isotope_list, 
                 overwrite = FALSE)

assay(spe, 'exprs') <- assay(spe, 'compexprs')
assay(spe, 'counts') <- assay(spe, 'compcounts')
assay(spe, 'compexprs') <- assay(spe, 'compcounts') <- NULL

# Output ----------------------------------------------------------
images <- readRDS("data/processed/analysis/images.rds")

dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'cells'), showWarnings = F)

plotPixels(images, spe, masks, 'cell_id', 'img_id',
           colour_by = c('FITC-ITGAX', 'CD68'),
           bcg = list(`FITC-ITGAX` = c(0,5,1),
                      CD68 = c(0,5,1)),
           display = 'single',
           image_title = list(position = "topright",
                              margin = c(5,5),
                              font = 1,
                              cex = 1),
           legend = list(colour_by.title.font=4, colour_by.title.cex=4),
           save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'cells', 'clusters.png')))

# Distance to epithelium ------------------------------------------
epi_images <- loadImages(file.path("data/processed/extract/Epithelium/"), pattern='.tiff')

epi_coords <- list()
for(n in names(epi_images)){
  id <- str_replace(n, '_epithelium', '')
  
  imgdata <- imageData(epi_images[[n]])
  colnames(imgdata) <- 1:dim(imgdata)[2]
  
  epi_coords[[id]] <- as.data.frame(imgdata) %>% 
    mutate(x = 1:dim(imgdata)[1]) %>% 
    pivot_longer(1:dim(imgdata)[2], names_to = 'y', values_to = 'val') %>%
    mutate(y = as.numeric(y)) %>% filter(val == 1)
}

spe$distance_to_epithelium <- NA
spe$distance_to_epithelium_scaled <- NA

for(id in unique(spe$img_id)){
  epi <- t(epi_coords[[id]] %>% select(x, y))
  cell <- t(spatialCoords(spe)[spe$img_id == id,])
  
  dists <- sqrt(apply(cell, 2, function(center) {
    colSums((epi - center)^2)
  }))

  dists <- apply(dists, 2, min)
  spe$distance_to_epithelium[spe$img_id == id] <- dists
  spe$distance_to_epithelium_scaled[spe$img_id == id] <- dists/max(dists)
}

dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'distance_to_epithelium'), showWarnings = F)
plotCells(mask = masks, object = spe, cell_id = 'cell_id', img_id = 'img_id', 
          image_title = list(position = "topright",
                             margin = c(5,5),
                             font = 1,
                             cex = 1),
          display = 'single',
          colour_by = 'distance_to_epithelium',
          colour = list('distance_to_epithelium' = rev(viridis::viridis(100))), # Reversed - closer are yellow
          save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'distance_to_epithelium',
                                                'distance_to_epithelium.png')))

dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'distance_to_epithelium_scaled'), showWarnings = F)
plotCells(mask = masks, object = spe, cell_id = 'cell_id', img_id = 'img_id', 
          image_title = list(position = "topright",
                             margin = c(5,5),
                             font = 1,
                             cex = 1),
          display = 'single',
          colour_by = 'distance_to_epithelium_scaled',
          colour = list('distance_to_epithelium_scaled' = rev(viridis::viridis(100))), # Reversed - closer are yellow
          save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'distance_to_epithelium_scaled',
                                                'distance_to_epithelium_scaled.png')))

# GC --------------------------------------
gc_cells <- read.csv("data/processed/extract/CD68_ITGAX_cells/GC_Cell_objects.csv") %>% 
  filter(Location_Center_X != 'NaN') %>%
  left_join(IMAGE_METADATA, by='ImageNumber') %>%
  mutate(cell_id = paste0(Metadata_ID, '_', ObjectNumber))
  
spe$is_GC <- colnames(spe) %in% gc_cells$cell_id

dir.create(file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'GC_cells'), showWarnings = F)
plotCells(masks, spe, 'cell_id', 'img_id',
          display = 'single',
          colour_by = 'is_GC', colour = list(is_GC = c(`TRUE` = 'red', `FALSE` = 'lightgrey')),
          image_title = list(cex=5),
          save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_CD68ITGAX, 'GC_cells', 'gc_cells.png')))

# Save ------------------------------------
saveRDS(spe, "data/processed/analysis/CD69_ITGAX_segmented/spe_load_CD68_ITGAX.rds")
saveRDS(masks, "data/processed/analysis/CD69_ITGAX_segmented/masks_CD68_ITGAX.rds")
