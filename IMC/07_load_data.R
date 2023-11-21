# Overview ----------------------------------------------------------------

# 7_load_data
# Load in cell data, images and masks outputted from 
# IMCSegmentationPipeline

library(imcRtools)
library(dplyr)
library(cytomapper)
library(stringr)
library(dittoSeq)
library(patchwork)

source("config/globals.R")

DIRECTORIES$OUTPUT_LOAD <- file.path(DIRECTORIES$OUTPUT, '07_load')
dir.create(file.path(DIRECTORIES$OUTPUT_LOAD), showWarnings = F)

# Cell data ---------------------------------------------------------------

# Read in data
spe <- read_cpout(file.path(DIRECTORIES$CPOUT),
                  extract_metal_from = "Metal_tag")

# Add arcsinh transformed assay
assay(spe, "exprs") <- asinh(counts(spe)/1)

# Rename rows
rownames(spe) <- rowData(spe)$Actual

# Don't use DNA channels
rowData(spe)$use_channel <- c(rep(T, 35), rep(F,2))

# Attach metadata
spe$cell_id <- spe$ObjectNumber
spe$img_id <- paste0(spe$Metadata_acname, '_s0_a', spe$Metadata_acid)

colnames(spe) <- paste0(spe$img_id, '_', spe$cell_id)

df <- as.data.frame(colData(spe)) %>% left_join(SAMPLE_METADATA, by='img_id')
spe$pt_id <- df$PT_ID
spe$block_id <- df$Block_ID
spe$tissue <- df$Tissue

dir.create(file.path(DIRECTORIES$OUTPUT_LOAD, 'arcsinh_transform'), showWarnings = F)

for(m in PANEL$Actual){
  p <- (dittoRidgePlot(spe, var = m, group.by = 'Metadata_acname', assay='counts') + theme(legend.position = 'none')) |
  (dittoRidgePlot(spe, var = m, group.by = 'Metadata_acname', assay='exprs') + 
     theme(legend.position = 'none', axis.text.y = element_blank()) + 
     ggtitle(paste0(m, ' - arcsinh transformed')))
  ggsave(file.path(DIRECTORIES$OUTPUT_LOAD, 'arcsinh_transform', paste0(m, '.pdf')), p, width = 10, height = 5)
}

# Images ------------------------------------------------------------------
images <- loadImages(file.path(DIRECTORIES$IMAGES), pattern = '.tiff')

channelNames(images) <- rownames(spe)
names(images) <- mcols(images)$img_id <- str_replace(names(images), '_ac_full', '')

# Masks -------------------------------------------------------------------
masks <- loadImages(file.path(DIRECTORIES$MASKS), as.is = TRUE)

names(masks) <- mcols(masks)$img_id <- str_replace(names(masks), '_ac_ilastik_s2_Probabilities_mask', '')

all.equal(names(images), names(masks))

# Save --------------------------------------------------------------------
saveRDS(spe, file.path(DIRECTORIES$DATA_ANALYSIS, 'spe_load.rds'))
saveRDS(images, file.path(DIRECTORIES$DATA_ANALYSIS, 'images.rds'))
saveRDS(masks, file.path(DIRECTORIES$DATA_ANALYSIS, 'masks.rds'))

rowData(spe)
