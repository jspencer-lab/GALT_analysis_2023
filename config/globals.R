# Directories -------------------------------------------------------------------

DIRECTORIES <- list()
DIRECTORIES$DATA <- file.path("data")
DIRECTORIES$DATA_RAW <- file.path("data", "raw")
DIRECTORIES$DATA_PROCESSED <- file.path("data", "processed")
DIRECTORIES$DATA_SPILLOVER <- file.path("data", "spillover")
DIRECTORIES$OUTPUT <- file.path("output")

DIRECTORIES$CPOUT <- file.path(DIRECTORIES$DATA_PROCESSED, 'data_extract', 'cpout')
DIRECTORIES$IMAGES <- file.path(DIRECTORIES$CPOUT, 'images')
DIRECTORIES$MASKS <- file.path(DIRECTORIES$CPOUT, 'masks')
DIRECTORIES$ILASTIK <- file.path(DIRECTORIES$DATA_PROCESSED, 'data_extract', 'ilastik')

DIRECTORIES$DATA_ANALYSIS <- file.path(DIRECTORIES$DATA_PROCESSED, 'data_analysis')
dir.create(DIRECTORIES$DATA_ANALYSIS, showWarnings = F)

# Panel ------------------------------------------------------------------------
PANEL <- read.csv("config/IMC_panel.csv")

# Samples ----------------------------------------------------------------------
SAMPLE_METADATA <- read.csv("config/IMC_sample_metadata.csv")
# Get sample ID to match output from IMCSegmentationPipeline format (mcd_name _s0_ mcd_acID)
SAMPLE_METADATA$img_id <- paste0(SAMPLE_METADATA$MCD_name, '_s0_', SAMPLE_METADATA$AC_ID)

# Cell types -------------------------------------------------------------------
T_CELLS <- "T_Cell"
B_CELLS <- "B_Cell"
B_OR_T_CELLS <- "B_or_T_Cell"
EPITHELIAL_CELLS <- "Epithelial_Cell"
ENDOTHELIAL_CELLS <- "Endothelial_Cell"
MACDC <- "MacDC"
UNLABELLED <- 'Unlabelled'

# Colours -----------------------------------------------------------------------
COLOURS <- list()
COLOURS$CELL_TYPES <- as.character(c('seagreen3', 'dodgerblue2', 'cyan', 'orange', 'tomato', 'yellow', 'lightgrey'))
names(COLOURS$CELL_TYPES) <- c("T_Cell",
                               "B_Cell",
                               "B_or_T_Cell",
                               "Epithelial_Cell",
                               "Endothelial_Cell",
                               "MacDC",
                               "Unlabelled")
COLOURS$IMG_ID <- colorRampPalette(RColorBrewer::brewer.pal(5, 'Set3'))(16)
COLOURS$TISSUE <- RColorBrewer::brewer.pal(3, 'Set2')
COLOURS$MCD <- RColorBrewer::brewer.pal(3, 'Set1')

BLANK_THEME <- ggplot2::theme_minimal()+
  ggplot2::theme(
    axis.title = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    axis.text = element_blank()
  )
