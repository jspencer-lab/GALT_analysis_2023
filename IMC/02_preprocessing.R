
# Overview ----------------------------------------------------------------

# 2_preprocessing
# Filtering the image output from 1_imc_image_extract.ipynb to only
# those in the sample metadata

source("config/globals.R")

# Libraries ---------------------------------------------------------------


# Remove unneeded files ---------------------------------------------------

# Images
valid_images <- file.path(DIRECTORIES$IMAGES, c(paste0(SAMPLE_METADATA$ID, '_ac_full.csv'), 
                                                paste0(SAMPLE_METADATA$ID, '_ac_full.tiff')))
invalid_images <- setdiff(list.files(DIRECTORIES$IMAGES, full.names = T), valid_images)
unlink(invalid_images)

# Ilastik
valid_ilastik <- file.path(DIRECTORIES$ILASTIK, c(paste0(SAMPLE_METADATA$ID, '_ac_ilastik.csv'), 
                                                  paste0(SAMPLE_METADATA$ID, '_ac_ilastik.tiff')))
invalid_ilastik <- setdiff(list.files(DIRECTORIES$ILASTIK, full.names = T), valid_ilastik)
unlink(invalid_ilastik)

