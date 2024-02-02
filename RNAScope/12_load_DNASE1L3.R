library(cytomapper)
library(stringr)
library(imcRtools)
library(dplyr)
library(scoper)
library(scater)
library(CATALYST)
library(dittoSeq)

source("config/globals.R")

metadata <- read.csv("data/processed/extract/DNASE1L3_cells/Image.csv")


# Load images and masks -----------
masks <- loadImages(file.path("data/processed/extract/DNASE1L3_cells/masks"), pattern='_DNASE1L3_cells.tiff')
mcols(masks)$img_id <- names(masks) <- str_match(names(masks), "(.*?)\\s*_ac_full_DNASE1L3_cells")[,2]
masks <- scaleImages(masks, 2^16-1)
# Reorder
masks <- masks[metadata$Metadata_ID]

# Create spe ----------------
cell_data <- read.csv("data/processed/extract/DNASE1L3_cells/DNASE1L3_cells.csv") %>% 
  left_join(metadata %>% select(ImageNumber, Metadata_ID), by='ImageNumber') %>% 
  rename(img_id = 'Metadata_ID', cell_id = 'ObjectNumber') %>% 
  left_join(SAMPLE_METADATA %>% select(ID, tissue, patient_id, block_id, slide), by = c('img_id' = 'ID'))

int_data <- t(as.matrix(cell_data[, paste0('Intensity_MeanIntensity_Stack_', paste0('c', 1:35))]))*(2**16)
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

# Save ------------------------------------
saveRDS(spe, "data/processed/analysis/DNASE1L3_segmented/spe_DNASE1L3.rds")
saveRDS(masks, "data/processed/analysis/DNASE1L3_segmented/masks_DNASE1L3.rds")
