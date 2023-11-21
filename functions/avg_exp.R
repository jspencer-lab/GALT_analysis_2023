average_expression <- function(spe, grouping, assay = 'exprs', genes = rownames(spe), scale=F, method='mean'){
  mat <- assay(spe, assay)
  
  mat <- as.data.frame(t(assay(spe, assay)))[,genes]
  mat$grouping <- spe[[grouping]]
  if(method == 'mean'){
    mat <- mat %>% group_by(grouping) %>% summarise(across(everything(), mean))  
  } else if (method == 'median'){
    mat <- mat %>% group_by(grouping) %>% summarise(across(everything(), median))  
  }
  
  mat <- column_to_rownames(mat, 'grouping')
  if(scale){
    mat <- scale(mat)
  }
  return(mat)  
}
# 
# pheatmap::pheatmap(average_expression(bcells, 'clusters', genes=PANEL$Actual[PANEL$B_cell_marker], scale=T),
#                    breaks = seq(-3, 3, length.out = 101), 
#                    color = colorRampPalette(c('dodgerblue3', 'white', 'firebrick2'))(100))
