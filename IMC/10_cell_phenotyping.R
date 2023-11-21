# Overview ----------------------------------------------------------------

# 10_cell_phenotyping
# TODO

library(imcRtools)
library(dplyr)
library(cytomapper)
library(stringr)
library(ggplot2)
library(CATALYST)
library(dittoSeq)
library(viridis)
library(patchwork)
library(corrplot)
library(tidyverse)
library(ggrepel)
library(EBImage)
library(scuttle)
library(mclust)
library(scater)
library(scran)
library(dichromat)
library(caret)
library(doParallel)
library(ggridges)

source("config/globals.R")

DIRECTORIES$OUTPUT_PHENOTYPING <- file.path(DIRECTORIES$OUTPUT, '10_cell_phenotyping')
dir.create(DIRECTORIES$OUTPUT_PHENOTYPING, showWarnings = F)

# Load ----------------------------------------------------------------------
spe <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'spe_QC.rds'))
images <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'images.rds'))
masks <- readRDS(file.path(DIRECTORIES$DATA_ANALYSIS, 'masks.rds'))

# Gating --------------------------------------------------------------------
# cytomapperShiny(spe, masks, images, 'cell_id', 'img_id')

# Load gated cells ----------------------------------------------------------
gated_spes <- sapply(list.files("output/10_cell_phenotyping/1_gating/manual_gates", full.names = F), function(x){
  readRDS(file.path("output/10_cell_phenotyping/1_gating/manual_gates", x))}
)

# Get cells and their labels
df <- data.frame(filename = names(gated_spes),
                 img = as.character(lapply(gated_spes, function(x){ return(x$img_id[1])})),
                 cell = as.character(lapply(gated_spes, function(x){ return(x$cytomapper_CellLabel[1])})))
table(table(df$img) == 8) # Ensure all images appear 8 times (1 per cell type)
table(table(df$cell) == 16) # Ensure all cell types appear 16 times (1 per image)

# Merge gated_spes
cols <- names(colData(gated_spes[[1]]))
# Ensure all have the same colData columns and wipe rowData
gated_spes <- lapply(gated_spes, function(x){
   colData(x) <- colData(x)[,cols]
   rowData(x) <- NA
   return(x)
})
# Merge
labelled_spe <- do.call("cbind", gated_spes)
# Apply the rowData from the original spe
rowData(labelled_spe) <- rowData(spe)

# Construct dataframe and relabel - all MacDC subtypes as just MacDC
df <- data.frame(cell = colnames(labelled_spe),
                 label.1 = labelled_spe$cytomapper_CellLabel) %>%
  mutate(label.2 = ifelse(label.1 %in% c("MacDC_CD11b", "MacDC_CD11c", "MacDC_CD68"), MACDC, label.1))
 
# Calculate dupes
dupes <- names(which(table(df$cell) > 1))
true_dupes <- unique(df$cell[df$cell %in% dupes & df$label.2 != MACDC])

# Split df
df_dupes <- split(df, df$cell %in% true_dupes)$`TRUE`
df_nodupes <- split(df, df$cell %in% true_dupes)$`FALSE`

df_nodupes <- df_nodupes[!duplicated(df_nodupes$cell),]

# How many cells got labelled?
100*length(unique(df$cell)) / ncol(spe) # 
 
# Of the labelled cells, how many duplicates do we have?
100*length(unique(df_dupes$cell)) / length(unique(df$cell)) # <0.003% of labelled cells are duplicates
 
# How many cells get a label?
100* nrow(df_nodupes) / ncol(spe) # 
 
# Add labels to spe object
label_vector <- rep(UNLABELLED, ncol(spe))
names(label_vector) <- colnames(spe)
label_vector[df_nodupes$cell] <- df_nodupes$label.2

colData(spe)$classification_gated <- label_vector
 
# Gating QC ------------------------------------------------------------------
spe$classification_gated[spe$classification_gated == B_OR_T_CELLS] <- UNLABELLED

df <- as.data.frame(table(spe$classification_gated)) %>% rename(Cell = 'Var1')
p <- ggplot(df, aes(x='', y=Freq, fill=Cell)) + geom_bar(stat='identity') + 
  coord_polar(theta='y') + scale_fill_manual(values=COLOURS$CELL_TYPES) + 
  BLANK_THEME
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating', 'QC', 'pie_labelled.pdf'), p, height=2.5, width=3.3)

# Output
dir.create(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating'), showWarnings = F)
dir.create(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating', 'classification_gated'), showWarnings = F)
plotCells(masks, spe, 'cell_id', 'img_id',
          colour_by = 'classification_gated', colour = list(classification_gated = COLOURS$CELL_TYPES),
          display = 'single',
          image_title = list(position = "topright",
                             margin = c(5,5),
                             font = 1,
                             cex = 1),
          save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating', 'classification_gated',
                                                'classification_gated.png')))

dir.create(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating', 'QC'), showWarnings = F)
p <- as_tibble(colData(spe)) %>%
  group_by(img_id) %>%
  summarise(labelled_cells = 100*sum(classification_gated != UNLABELLED)/n(),
            number_cells = n()) %>%
  as.data.frame() %>% 
  ggplot(aes(x=img_id, y=labelled_cells, fill=img_id)) + geom_bar(stat = 'identity') + ylim(0,100) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = 'none') +
  scale_fill_manual(values = COLOURS$IMG_ID)
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating', 'QC', 'labelled_per_image.pdf'), p)

p <- as_tibble(colData(spe)) %>%
  group_by(Metadata_acname, img_id) %>%
  summarise(labelled_cells = 100*sum(classification_gated != UNLABELLED)/n(),
            number_cells = n()) %>%
  as.data.frame() %>% ggplot(aes(x=Metadata_acname, y=labelled_cells, fill=Metadata_acname)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = 'none')  +
  scale_fill_manual(values = COLOURS$MCD)
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating', 'QC', 'labelled_per_mcd.pdf'), p)

# Class per image
p <- as_tibble(colData(spe)) %>%
  group_by(classification_gated, img_id) %>%
  summarise(number_cells = n()) %>%
  as.data.frame() %>% 
  ggplot(aes(x=img_id, y=number_cells, color= classification_gated)) + geom_jitter() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = 'none')  +
  scale_color_manual(values = COLOURS$CELL_TYPES)
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating', 'QC', 'celltypes_per_image.pdf'), p)


# Heatmaps - sampled down to avoid too many cells
names(COLOURS$IMG_ID) <- unique(spe$img_id)
markers <- c('CD20', 'CD31', 'E_Cadherin', 'CD68', 'CD11c', 'CD11b', 'CD3')

set.seed(101)
scaled_spe <- spe[markers, sample(colnames(spe), 2000)]
assay(scaled_spe, "scaled_exprs") <- t(scale(t(assay(scaled_spe, 'exprs'))))

p <- dittoHeatmap(scaled_spe, assay = "exprs",
             annot.by = c("classification_gated"), 
             order.by = "classification_gated", cluster_rows = FALSE,
             scale = "none", heatmap.colors = viridis(100), 
             annotation_colors = list(classification_gated = COLOURS$CELL_TYPES,
                                      img_id = COLOURS$IMG_ID))
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating', 'QC', 'celltypes_heatmap.pdf'), p)

p <- dittoHeatmap(scaled_spe, assay = "scaled_exprs",
             annot.by = c("classification_gated"), 
             order.by = c("classification_gated"), cluster_rows = FALSE,
             annotation_colors = list(classification_gated = COLOURS$CELL_TYPES,
                                      img_id = COLOURS$IMG_ID),
             heatmap.colors = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-3, 3, length.out = 101))
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '1_gating', 'QC', 'celltypes_heatmap_scaled.pdf'), p)

# Training  ------------------------------------------------------------------

dir.create(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '2_training'), showWarnings = F)

labelled_spe <- spe[,spe$classification_gated != UNLABELLED]
unlabelled_spe <- spe[,spe$classification_gated == UNLABELLED]

# TODO: currently downsampled
set.seed(101)
labelled_spe <- labelled_spe[,sample(seq_len(ncol(labelled_spe)), 50000)]

set.seed(101)
trainIndex <- createDataPartition(factor(labelled_spe$classification_gated), p = 0.70)
train_spe <- labelled_spe[,trainIndex$Resample1]
test_spe <- labelled_spe[,-trainIndex$Resample1]

# Define seeds for parallel processing
# Per iteration, we evaluate 10 models while tuning mtry
set.seed(101)
seeds <- vector(mode = "list", length = 11)
for (i in 1:10) {
  seeds[[i]] <- sample.int(5000, 10)
}

seeds[[11]] <- sample.int(5000, 1)

fitControl <- trainControl(method = "repeatedcv",
                           repeats = 1,
                           number = 10,
                           seeds = seeds)

# TODO reduce num of cores if needed - difference between virtual (16) and physical (8)?
cl <- makePSOCKcluster(round(detectCores()/2,0), outfile = "")
registerDoParallel(cl)

set.seed(101)
start = Sys.time()
rffit <- train(x = t(assay(train_spe, "exprs")[rowData(spe)$use_channel,]), 
               y = factor(train_spe$classification_gated),
               method = "rf", ntree = 500,
               tuneLength = 10,
               trControl = fitControl,
               allowParallel = TRUE)
stopCluster(cl)
end = Sys.time()
print(end-start)

saveRDS(rffit, 'rffit.rds')

# Plots 
p <- ggplot(rffit) + 
  geom_errorbar(data = rffit$results,
                aes(ymin = Accuracy - AccuracySD,
                    ymax = Accuracy + AccuracySD),
                width = 0.4)
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '2_training', 'rffit_error.pdf'), p)

pdf(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '2_training', 'rffit_variable_importance.pdf'))
plot(varImp(rffit))
dev.off()

confusionMatrix(rffit)

# Testing -------------------------------------------------------------------------------
dir.create(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '3_testing'), showWarnings = F)

cur_mat <- t(assay(test_spe, "exprs")[rowData(test_spe)$use_channel,])

# Predict cell phenotypes in test data
cur_pred <- predict(rffit, 
                    newdata = cur_mat)

cm <- confusionMatrix(data = cur_pred, 
                      reference = factor(test_spe$classification_gated), 
                      mode = "everything")

cm

p <- data.frame(cm$byClass) %>%
  mutate(class = sub("Class: ", "", rownames(cm$byClass))) %>%
  ggplot() + 
  geom_point(aes(1 - Specificity, Sensitivity, 
                 size = Detection.Rate,
                 fill = class),
             shape = 21) + 
  theme_classic(base_size = 15) + 
  ylab("Sensitivity (TPR)") +
  xlab("1 - Specificity (FPR)")
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '3_testing', 'sensivity_vs_specificity.pdf'), p, width=8, height=8)

cur_pred <- predict(rffit, 
                    newdata = cur_mat, 
                    type = "prob")
cur_pred$truth <- factor(test_spe$classification_gated)

p <- cur_pred %>%
  pivot_longer(cols = B_Cell:T_Cell) %>%
  ggplot() +
  geom_boxplot(aes(x = name, y = value, fill = name), outlier.size = 0.5) +
  facet_wrap(. ~ truth, ncol = 1) + 
  theme(panel.background = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '3_testing', 'probs_per_class.pdf'), p, width=5, height=10)

# Classify ---------------------------------------------------------------------------
dir.create(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '4_classify'), showWarnings = F)

# Select unlabeled data
cur_mat <- t(assay(unlabelled_spe, "exprs")[rowData(unlabelled_spe)$use_channel,])

# Predict cell phenotypes
cell_class <- as.character(predict.train(rffit, 
                                         newdata = cur_mat, 
                                         type = "raw"))
names(cell_class) <- rownames(cur_mat)

table(cell_class)

# Extract classification probabilities
cell_prob <- predict.train(rffit, 
                           newdata = cur_mat, 
                           type = "prob")



# Distribution of maximum probabilities
tibble(max_prob = rowMax(as.matrix(cell_prob)),
       type = cell_class) %>%
  ggplot() +
  geom_density_ridges(aes(x = max_prob, y = cell_class, fill = cell_class)) +
  scale_fill_manual(values = COLOURS$CELL_TYPES) +
  theme_classic(base_size = 15) +
  xlab("Maximum probability") +
  ylab("Cell type") + 
  xlim(c(0,1.2))

cell_class[rowMax(as.matrix(cell_prob)) < 0.4] <- UNLABELLED

# Store labels in SpatialExperiment onject
cell_labels <- spe$classification_gated
cell_labels[colnames(unlabelled_spe)] <- cell_class
spe$cell_type <- cell_labels 

m <- table(spe$cell_type, spe$img_id)
p <- pheatmap::pheatmap(scale(m), cluster_cols = F, cluster_rows = F)
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '4_classify', 'class_per_image.pdf'), p)

set.seed(101)
scaled_spe <- spe[rowData(spe)$use_channel, sample(colnames(spe), 2000)]
assay(scaled_spe, "scaled_exprs") <- t(scale(t(assay(scaled_spe, 'exprs'))))

# dittoHeatmap(scaled_spe, assay = "exprs",
#                   annot.by = c("cell_type"), 
#                   order.by = "cell_type", cluster_rows = FALSE,
#                   scale = "none", heatmap.colors = viridis(100), 
#                   annotation_colors = list(cell_type = COLOURS$CELL_TYPES,
#                                            img_id = COLOURS$IMG_ID))
p <- dittoHeatmap(scaled_spe, assay = "scaled_exprs",
                  annot.by = c("cell_type"), 
                  order.by = c("cell_type"), cluster_rows = T,
                  annotation_colors = list(cell_type = COLOURS$CELL_TYPES,
                                           img_id = COLOURS$IMG_ID),
                  heatmap.colors = colorRampPalette(c("dodgerblue3", "white", "red"))(100),
                  breaks = seq(-3, 3, length.out = 101))
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '4_classify', 'heatmap_per_cell_type.pdf'), p)


dir.create(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '4_classify',  "cell_type_overlay"), showWarnings = F)
plotCells(masks, spe, 'cell_id', 'img_id',
          colour_by = 'cell_type', colour = list(cell_type = COLOURS$CELL_TYPES),
          display = 'single',
          image_title = list(position = "topright",
                             margin = c(5,5),
                             font = 1,
                             cex = 1),
          scale_bar = list(length=500, label="", lwidth=10),
          save_plot = list(filename = file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '4_classify', 'cell_type_overlay',
                                                'cell_type.png')))

df <- as.data.frame(table(spe$cell_type)) %>% rename(Cell = 'Var1')
p <- ggplot(df, aes(x='', y=Freq, fill=Cell)) + geom_bar(stat='identity') + 
  coord_polar(theta='y') + scale_fill_manual(values=COLOURS$CELL_TYPES) + 
  BLANK_THEME
ggsave(file.path(DIRECTORIES$OUTPUT_PHENOTYPING, '4_classify', 'pie_classified.pdf'), p, height=2.5, width=3.3)

saveRDS(spe, file.path(DIRECTORIES$DATA_ANALYSIS, 'spe_classified_v2.rds'))