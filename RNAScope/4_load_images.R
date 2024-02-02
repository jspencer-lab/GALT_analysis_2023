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
images <- loadImages(file.path("data/processed/extract/cpout/images/"), pattern='_full.tiff')
mcols(images)$img_id <- names(images) <- str_match(names(images), "(.*?)\\s*_ac_full")[,2]
channelNames(images) <- PANEL$Actual
# Reorder
images <- images[metadata$Metadata_ID]


# Files for GC masking -----------
dir.create(file.path(DIRECTORIES$DATA_EXTRACT, 'GC'), showWarnings = F)
dir.create(file.path(DIRECTORIES$DATA_EXTRACT, 'GC', 'reference'), showWarnings = F)

plotPixels(images, colour_by = c('E-Cadherin', 'CD45RB', 'Ki-67'),
           colour = list(`E-Cadherin` = c('black', 'orange'),
                         `CD45RB` = c('black', 'magenta'),
                         `Ki-67` = c('black', 'white')),
           bcg = list(`E-Cadherin` = c(0,5,1),
                      `CD45RB` = c(0,20,1),
                      `Ki-67` = c(0,20,1)),
           display = 'single',
           scale_bar = NULL,
           image_title = NULL,
           save_plot = list(filename = file.path(DIRECTORIES$DATA_EXTRACT, 'GC', 'reference', 'GC_reference.png')))

# Save -----------------------------
saveRDS(images, "data/processed/analysis/images.rds")
