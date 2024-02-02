source("config/globals.R")

# Remove unneeded files ----
files_to_remove <- setdiff(list.files(DIRECTORIES$DATA_IMAGES,  full.names = T), 
        file.path(DIRECTORIES$DATA_IMAGES, c(paste0(SAMPLE_METADATA$ID, '_full.csv'), 
                                             paste0(SAMPLE_METADATA$ID, '_full.tiff'))))

file.remove(files_to_remove)

files_to_remove <- setdiff(list.files(DIRECTORIES$DATA_ILASTIK,  full.names = T), 
                           file.path(DIRECTORIES$DATA_ILASTIK, 
                                     c(paste0(SAMPLE_METADATA$ID, '_ilastik.csv'), 
                                       paste0(SAMPLE_METADATA$ID, '_ilastik.tiff'))))

file.remove(files_to_remove)
