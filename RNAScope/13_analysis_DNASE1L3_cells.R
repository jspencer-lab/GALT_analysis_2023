library(cytomapper)
library(stringr)
library(imcRtools)
library(dplyr)
library(scoper)
library(scater)
library(tidyr)

source("config/globals.R")

DIRECTORIES$OUTPUT_DNASE1L3 <- file.path("output/DNASE1L3_segmented")
dir.create(DIRECTORIES$OUTPUT_DNASE1L3, showWarnings = F)

# Load ------------------------------------------
spe <- readRDS("data/processed/analysis/DNASE1L3_segmented/spe_DNASE1L3.rds")
images <- readRDS("data/processed/analysis/images.rds")
masks <- readRDS("data/processed/analysis/DNASE1L3_segmented/masks_DNASE1L3.rds")

# # Analysis --------------------------------------
# exp <- as.data.frame(t(assay(spe, 'exprs')))
# exp_pl <- pivot_longer(exp, colnames(exp), names_to = 'marker', values_to = 'intensity')
# 
# ggplot(exp_pl, aes(x=marker, y=intensity)) + geom_boxplot(fill='lightgrey') + theme_classic() +
#   theme(axis.text.x = element_text(angle=45, hjust=1))
# 
# ggplot(exp_pl %>% filter(marker %in% c('CD11b', 'CD20', 'CD3', 'CD31', 'CD68', 'FITC-ITGAX', 'E-Cadherin')), 
#        aes(x=marker, y=intensity)) + geom_boxplot(fill='lightgrey', outlier.shape = NA) + theme_classic() +
#   theme(axis.text.x = element_text(angle=45, hjust=1))
# 
# 
# exp2 <- exp %>% mutate(img_id = spe$img_id) %>% pivot_longer(colnames(exp), names_to = 'marker', values_to = 'intensity') %>%
#   group_by(img_id, marker) %>% summarise(mean_int_per_image = mean(intensity))
# ggplot(exp2 %>% filter(marker %in% c('CD11b', 'CD20', 'CD3', 'CD31', 'CD68', 'FITC-ITGAX', 'E-Cadherin')), 
#                        aes(x=marker, y=mean_int_per_image)) + geom_boxplot(fill='lightgrey') + theme_classic() +
#   theme(axis.text.x = element_text(angle=45, hjust=1))

# SED ---------------------------------------------
sed_cells <- bind_rows(lapply(list.files("data/processed/analysis/DNASE1L3_segmented/SED", full.names = T), 
                              FUN = read.csv))$ids

spe$is_sed <- colnames(spe) %in% sed_cells

df_lys <- data.frame(img_id = spe$img_id,
                     is_sed = ifelse(spe$is_sed, 'is_sed', 'not_sed'),
                     lysozyme = assay(spe, 'exprs')['Lysozyme', ]) %>% 
  group_by(img_id, is_sed) %>% 
  summarise(avg_lyz = mean(lysozyme)) %>% 
  tidyr::pivot_wider(names_from = 'is_sed', values_from = 'avg_lyz')
write.csv(df_lys, "output/DNASE1L3_segmented/lysozyme.csv")

df_NOX2 <- data.frame(img_id = spe$img_id,
                     is_sed = ifelse(spe$is_sed, 'is_sed', 'not_sed'),
                     NOX2 = assay(spe, 'exprs')['NOX2', ]) %>% 
  group_by(img_id, is_sed) %>% 
  summarise(avg_NOX2 = mean(NOX2)) %>% 
  tidyr::pivot_wider(names_from = 'is_sed', values_from = 'avg_NOX2')
write.csv(df_NOX2, "output/DNASE1L3_segmented/NOX2.csv")
