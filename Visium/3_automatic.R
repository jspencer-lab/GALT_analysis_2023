library(BayesSpace)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(harmony)
library(ggpubr)
library(RcmdrMisc)
library(corrplot)
library(viridis)

source("config/globals.R")
source("function/SpatialDimPlotMJP.R")

DIRECTORIES$OUTPUT_AUTOMATIC <- file.path(DIRECTORIES$OUTPUT, '3_automatic_clustering')
dir.create(DIRECTORIES$OUTPUT_AUTOMATIC, showWarnings = F)

crop_xmin <- 1850
crop_xmax <- 2400
crop_ymin <- 2400
crop_ymax <- 2700

# Functions ----
convert_seurat_ST_to_sce <- function(seu){
  sce <- as.SingleCellExperiment(seu, assay='Spatial')
  # Add coordinates, correct imagerow and imagecol for issue where images have been re-oriented in SpaceRanger
  SummarizedExperiment::colData(sce)$row <- seu@images[[names(seu@images)[1]]]@coordinates$row
  SummarizedExperiment::colData(sce)$col <- seu@images[[names(seu@images)[1]]]@coordinates$col
  SummarizedExperiment::colData(sce)$imagecol <- seu@images[[names(seu@images)[1]]]@coordinates$imagerow
  SummarizedExperiment::colData(sce)$imagerow <- seu@images[[names(seu@images)[1]]]@coordinates$imagecol
  
  return(sce)
}

# Load ----
spatial_seurat <- readRDS(file.path(DIRECTORIES$DATA_PROCESSED, '1_load_seurat.rds'))

# Convert to SPE ----
spes <- sapply(names(spatial_seurat), function(x){
  convert_seurat_ST_to_sce(spatial_seurat[[x]])
})

# Pre-processing for BayesSpace ----
for(id in names(spes)){
  set.seed(101)
  spes[[id]] <- spatialPreprocess(spes[[id]], platform="Visium", n.PCs=15, n.HVGs=2000, 
                                  log.normalize=T)
  
  spes[[id]] <- qTune(spes[[id]], qs=seq(2, 10), platform="Visium")
  
  p <- qPlot(spes[[id]]) +  ggtitle(id)
  ggsave(file.path(DIRECTORIES$OUTPUT_AUTOMATIC, paste0('q_plot_', id, '.pdf')), p)  
}

# Spatial cluster ----
# Number of clusters based on q plots
qs <- list('APPGUT_A1'=7, 'APPGUT_D1'=8)

for (id in names(spes)){
  set.seed(101)
  spes[[id]] <- spatialCluster(spes[[id]], q=qs[[id]], platform="Visium",
                               init.method="mclust", model="t",
                               nrep=10000, burn.in=500)
  
  Idents(spatial_seurat[[id]]) <- spatial_seurat[[id]]$spatial_cluster <- 
    factor(spes[[id]]$spatial.cluster, levels=1:max(spes[[id]]$spatial.cluster))
  
}

# Plots ----
dir.create(file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'spatial_clusters'), showWarnings = F)
cols <- c('brown', 'cyan', 'limegreen', 'dodgerblue2', 'pink', 'red', 'gold', 'black')
for(id in names(spatial_seurat)){
  p <- SpatialPlotMJP(spatial_seurat[[id]], fill = 'spatial_cluster', pt_size = 3) + 
    scale_color_manual(values = cols) + theme(legend.position = 'none')
  filename <- file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'spatial_clusters', paste0(id, "_spatial_cluster.png"))
  ggsave(filename, p, width = 10, height=10)
  
  if(id == 'APPGUT_A1'){
    p <- png::readPNG(filename)
    q <- imagefx::crop.image(p, crop_ymin, crop_xmin, crop_ymax, crop_xmax)
    png::writePNG(q$img.crop, file.path(DIRECTORIES$OUTPUT_AUTOMATIC, "spatial_clusters",
                                        "APPGUT_A1_spatial_cluster_CROP.png"))
  }
  
  p <- SpatialPlotMJP(spatial_seurat[[id]], fill = 'spatial_cluster', pt_size = 3) + 
    scale_color_manual(values = cols, name= 'Spatial \ncluster') + theme(legend.position = 'bottom')
  filename <- file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'spatial_clusters', 
                        paste0(id, "_spatial_cluster_legend.png"))
  ggsave(filename, p, width = 5, height=5)
}

# Lymphoid assignment ----
dir.create(file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'is_lymphoid'), showWarnings = F)
lymphoid <- list('APPGUT_A1'=c(7), 'APPGUT_D1'=c(6))


for(id in names(spes)){
  spes[[id]]$lymphoid <- spes[[id]]$spatial.cluster %in% lymphoid[[id]]
  spatial_seurat[[id]]$lymphoid <- ifelse(spatial_seurat[[id]]$spatial_cluster %in% lymphoid[[id]],
                                          'lymphoid_tissue',
                                          'other')
  
  filename <- file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'is_lymphoid', paste0(id, '_is_lymphoid.png'))
  p <- SpatialPlotMJP(spatial_seurat[[id]], 'lymphoid', pt_size=3) + 
    theme(legend.position = 'none') +
    scale_color_manual(values = c('firebrick2', 'darkgrey'))
  ggsave(filename, p, width=10, height=10)
  
  if(id == 'APPGUT_A1'){
    p <- png::readPNG(filename)
    q <- imagefx::crop.image(p, crop_ymin, crop_xmin, crop_ymax, crop_xmax)
    png::writePNG(q$img.crop, file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'is_lymphoid',
                                        paste0(id, '_is_lymphoid_CROP.png')))
  }
  
  filename <- file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'is_lymphoid', paste0(id, '_is_lymphoid_legend.png'))
  p <- SpatialPlotMJP(spatial_seurat[[id]], 'lymphoid', pt_size=3) + 
    theme(legend.position = 'bottom') +
    scale_color_manual(values = c('firebrick2', 'darkgrey'))
  ggsave(filename, p, width=10, height=10)
}


# Filter ----
lymphoid_seurat <- sapply(names(spatial_seurat), function(id){
  subset(spatial_seurat[[id]], subset = lymphoid == 'lymphoid_tissue')
})
lymphoid_spes <- sapply(names(spes), function(id){
  spes[[id]][,spes[[id]]$lymphoid]
})

# Checkpoint ----
saveRDS(spatial_seurat, file.path(DIRECTORIES$DATA_PROCESSED, '3_seurat.rds'))
saveRDS(lymphoid_seurat, file.path(DIRECTORIES$DATA_PROCESSED, '3_lymphoid_seurat.rds'))
saveRDS(spes, file.path(DIRECTORIES$DATA_PROCESSED, '3_spes.rds'))
saveRDS(lymphoid_spes, file.path(DIRECTORIES$DATA_PROCESSED, '3_lymphoid_spes.rds'))

# Lymphoid subcluster ----
spatial_seurat <- readRDS(file.path(DIRECTORIES$DATA_PROCESSED, '3_seurat.rds'))
lymphoid_seurat <- readRDS(file.path(DIRECTORIES$DATA_PROCESSED, '3_lymphoid_seurat.rds'))
spes <- readRDS(file.path(DIRECTORIES$DATA_PROCESSED, '3_spes.rds'))
lymphoid_spes <- readRDS(file.path(DIRECTORIES$DATA_PROCESSED, '3_lymphoid_spes.rds'))


lymphoid_seurat <- sapply(names(spatial_seurat), function(x){
  q <- SCTransform(lymphoid_seurat[[x]], assay = 'Spatial', method = "glmGamPoi", 
                   return.only.var.genes = F, min_cells=1)
  return(q)
})

# Merge
var_feat <- SelectIntegrationFeatures(lymphoid_seurat, nfeatures=3000)
lymphoid_data_merged <- merge(lymphoid_seurat$APPGUT_A1, lymphoid_seurat$APPGUT_D1)
VariableFeatures(lymphoid_data_merged) <- var_feat

# Harmonise
lymphoid_data_merged <- RunPCA(lymphoid_data_merged)
lymphoid_data_merged <- RunHarmony(lymphoid_data_merged, group.by.vars = 'orig.ident', assay.use = 'SCT')

p <- DimPlot(lymphoid_data_merged, reduction = 'harmony', group.by = 'orig.ident') + ggtitle("") + 
  scale_color_manual(values = c('purple', 'seagreen2'))
ggsave("output/Supplementary/harmony.pdf", p, height = 5, width = 7)

# Spatial clustering ----
DefaultAssay(lymphoid_data_merged) <- 'SCT'
lymphoid_data_merged <- PrepSCTFindMarkers(lymphoid_data_merged)

lymphoid_data_merged <- FindNeighbors(lymphoid_data_merged, reduction='harmony', dims=1:15) %>% FindClusters()
lymphoid_data_merged <- RenameIdents(lymphoid_data_merged, 
                                     '0' = 'B cell/T cell', 
                                     '1' = 'GC 1', 
                                     '2' = 'GC 2', 
                                     '3' = 'T cell zone', 
                                     '4' = 'SED/Epithelium')

lymphoid_data_merged$seurat_clusters <- Idents(lymphoid_data_merged)


SpatialDimPlot(lymphoid_data_merged)

# Reattach to individual objects for visuals
dir.create(file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'subcluster'), showWarnings = F)
for(id in names(lymphoid_seurat)){
  # df <- lymphoid_data_merged@meta.data[lymphoid_data_merged$orig.ident == id, 'seurat_clusters']
  # lymphoid_seurat[[id]]$subcluster <- df
  
  p <- SpatialPlotMJP(lymphoid_seurat[[id]], 'subcluster', pt_size = 3) + 
    scale_color_manual(values = COLOURS$SUBCLUSTER) + theme(legend.position = 'none')
  
  filename <- file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'subcluster', paste0(id, '_subclusters.png'))
  ggsave(filename, p, width=10, height=10)

  if(id == 'APPGUT_A1'){
    p <- png::readPNG(filename)
    q <- imagefx::crop.image(p, crop_ymin, crop_xmin, crop_ymax, crop_xmax)
    png::writePNG(q$img.crop, file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'subcluster',
                                        paste0(id, '_subclusters_CROP.png')))
  }

  p <- SpatialPlotMJP(lymphoid_seurat[[id]], 'subcluster', pt_size = 3) +
    scale_color_manual(values = COLOURS$SUBCLUSTER) +
    theme(legend.position = 'bottom')

  filename <- file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'subcluster', paste0(id, '_subclusters_legend.png'))
  ggsave(filename, p, width=10, height=10)
}

p <- SpatialPlotMJP(lymphoid_seurat$APPGUT_A1, 'subcluster', pt_size = 3) + 
  scale_color_manual(values = COLOURS$SUBCLUSTER) + 
  theme(legend.position = 'bottom')
ggsave('legend.pdf', p, height = 5, width=10)


# Markers
m <- FindAllMarkers(lymphoid_data_merged, assay = 'SCT', only.pos = T)
rownames(m) <- 1:dim(m)[1]
write.csv(m, file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'subcluster', 'markers.csv'))

DotPlot(lymphoid_data_merged, 'SCT', list(GC = c('MS4A1', 'TCL1A', 'SERPINA9', 'BCL6', 'BCL7A'),
                                       GC_diff = c('ANP32E', 'FTL', 'HIST1H4A', 'UCP2', 'LTF'), # GC diff genes
                                       T_cell = c('CD3E', 'CD4', 'IL7R', 'CD69', 'CCL21'),
                                       others = c('CYBB', 'IGHG1'),
                                       SED = c('CDH1', 'CCL20', 'CCL23', 'ITGAX'),
                                       lupus = c('DNASE1L3', 'C1QA', "C3", "C1R"))) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, hjust=1))

levels(lymphoid_data_merged) <- levels(lymphoid_data_merged)[5:1]
p <- DotPlot(lymphoid_data_merged, 'SCT', c('MS4A1', 'TCL1A', 'SERPINA9', 'BCL6', 'BCL7A',
                                       'CD3E', 'CD4', 'IL7R', 'CD69', 'CCL21',
                                       'CDH1', 'CCL20', 'CCL23', 'ITGAX',
                                       'DNASE1L3', 'C1QA', "C3", "C1R")) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'subcluster', 'dotplot.pdf'),p, width=8.5, height=6)

# 
# m_gc <- FindMarkers(lymphoid_data_merged, ident.1 = 'GC 1', ident.2 = 'GC 2')
# 'MS4A1' %in% rownames(m_gc)
# 
# m_unknown <- FindMarkers(lymphoid_data_merged, ident.1 = 'Unknown', ident.2 = c('GC 1', 'GC 2'), 
#                          only.pos = T)

# Lupus ----
dir.create(file.path(DIRECTORIES$OUTPUT_AUTOMATIC, 'lupus'),  showWarnings = F)
lupus_genes <- read.csv("reference/lupus_related_genes.csv")

lupus_genes_compliment <- lupus_genes$gene[lupus_genes$type == 'compliment']
lupus_genes_celldeath <- lupus_genes$gene[lupus_genes$type == 'cell_death']

# Correlations ----
df <- as.data.frame(t(lymphoid_data_merged@assays$SCT@data))
df <- df[,c(lupus_genes_compliment, lupus_genes_celldeath)]

corr <- rcorr.adjust(df, type="spearman")
corr$P[corr$P == "<.0001"] <- 0.00009
corr$P<- data.matrix(corr$P)
class(corr$P) <- "numeric"
pdf(file.path(DIRECTORIES$OUTPUT_AUTOMATIC, "lupus", "spatial_correlations.pdf"))
corrplot(corr$R$r, type="lower", order="original", method = "color", 
         p.mat = corr$P, sig.level = c(0.0001, .001, .01, .05), 
         pch.cex = 1, tl.cex = 1.2, cl.cex = 1.1,
         col = colorRampPalette(c('dodgerblue3', 'white', 'firebrick2'))(100),
         insig = "label_sig", tl.col = "black", tl.srt = 90, outline = "gray", mar=c(0,0,1,0), diag = FALSE)
dev.off()

write.csv(corr$R$r, file.path(DIRECTORIES$OUTPUT_AUTOMATIC, "lupus", "corr_R.csv"))
write.csv(corr$P, file.path(DIRECTORIES$OUTPUT_AUTOMATIC, "lupus", "corr_P.csv"))

pdf(file.path(DIRECTORIES$OUTPUT_AUTOMATIC, "lupus", "spatial_correlations_R.pdf"))
corrplot(corr$R$r, type="lower", order="original", method = "number", 
         pch.cex = 1, tl.cex = 1.2, cl.cex = 1, number.cex=0.75,
         col = colorRampPalette(c('black', 'black', 'black'))(100),
         tl.col = "black", outline = "black", diag = FALSE)
dev.off()

pdf(file.path(DIRECTORIES$OUTPUT_AUTOMATIC, "lupus", "spatial_correlations_P.pdf"))
corrplot(corr$P, type="lower", order="original", method = "number", 
         pch.cex = 1, tl.cex = 1.2, cl.cex = 1, number.cex=0.65, number.digits=4,
         col = colorRampPalette(c('black', 'black', 'black'))(100),
         tl.col = "black", outline = "black", diag = FALSE)
dev.off()

#Supp
dir.create("output/Supplementary",showWarnings = F)

p <- qPlot(spes$APPGUT_A1) + ggtitle("APP1") + theme_classic() 
ggsave("output/Supplementary/app1_qplot.pdf", width=3, height = 3)
p <- qPlot(spes$APPGUT_D1) + ggtitle("APP2") + theme_classic() 
ggsave("output/Supplementary/app2_qplot.pdf", width=3, height = 3)


lymphoid_data_merged <- RunUMAP(lymphoid_data_merged, reduction = 'harmony', dims=1:15)
p1 <- ggplot(data.frame(UMAP_1=lymphoid_data_merged@reductions$umap@cell.embeddings[,1],
                  UMAP_2=lymphoid_data_merged@reductions$umap@cell.embeddings[,2],
                  subcluster = lymphoid_data_merged$seurat_clusters),
       aes(x=UMAP_1, y=UMAP_2, fill=subcluster)) + geom_point(shape=21, size=3) + 
  scale_fill_manual(values=COLOURS$SUBCLUSTER) + theme_classic() + NoLegend()
ggsave(file.path("output/Supplementary", 'umap.pdf'), p1, width = 4, height=4)

plots <- list()
for(gene in c('DNASE1L3', 'C1QA', 'C1R', 'C3', 'MS4A1', 'CD3E', 'CCL20')){
  plots[[gene]] <- ggplot(data.frame(UMAP_1=lymphoid_data_merged@reductions$umap@cell.embeddings[,1],
                                     UMAP_2=lymphoid_data_merged@reductions$umap@cell.embeddings[,2],
                                     g = GetAssayData(lymphoid_data_merged, 'data')[gene,]),
               aes(x=UMAP_1, y=UMAP_2, col=g)) + scale_color_viridis(name=gene) + 
    geom_point(size=3) + theme_classic()
  ggsave(file.path("output/Supplementary", paste0('umap_', gene, '.pdf')), plots[[gene]], 
         width = 4, height=4)
}

saveRDS(lymphoid_seurat, "data/processed/3_lymphoid_seurat_clustered.rds")
saveRDS(lymphoid_data_merged, "data/processed/3_lymphoid_data_emerged_clustered.rds")

lymphoid_seurat <- readRDS("data/processed/3_lymphoid_seurat_clustered.rds")
lymphoid_data_merged <- readRDS("data/processed/3_lymphoid_data_emerged_clustered.rds")

SpatialPlotMJP(lymphoid_seurat$APPGUT_A1, 'subcluster', pt_size = 3) + 
  scale_color_manual(values = COLOURS$SUBCLUSTER)

GC <- subset(lymphoid_data_merged, idents = c('GC 1', 'GC 2'))

GC <- RunUMAP(GC, reduction = 'harmony', dims=1:15)
DimPlot(GC)

FeaturePlot(lymphoid_data_merged, 'CXCL13')
