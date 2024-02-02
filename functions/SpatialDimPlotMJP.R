SpatialPlotMJP <- function(seurat_obj, fill=NULL, border=NULL, 
                        alpha_fill = FALSE, forced_alpha='none', assay='SCT', 
                        pt_size=2, xlim_min = NULL, xlim_max = NULL, ylim_min = NULL, 
                        ylim_max = NULL){
  
  df <- as.data.frame(seurat_obj@images$slice1@coordinates)
  if(!is.null(fill)){
    if(fill %in% rownames(seurat_obj)){
      df[[fill]] <- GetAssayData(seurat_obj, assay = assay)[fill,]
    } else {
      df[[fill]] <- seurat_obj@meta.data[[fill]]  
    }
  }

  if (!is.null(border)){
    if(border != 'black'){
      df[[border]] <- seurat_obj@meta.data[[border]]  
    }
  }
  
  if(is.null(xlim_min)){
    xlim_min <- seurat_obj@misc$x_alignment_min
  }
  if(is.null(xlim_max)){
    xlim_max <- seurat_obj@misc$x_alignment_max
  }
  if(is.null(ylim_min)){
    ylim_min <- seurat_obj@misc$y_alignment_min
  }
  if(is.null(ylim_max)){
    ylim_max <- seurat_obj@misc$y_alignment_max
  }
  
  p <- ggplot(df, aes(x=imagecol, y=imagerow)) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
    ggpubr::background_image(seurat_obj@misc$image) 

  # TODO - better way of doing this
  if(is.null(border)){
    p <- p + geom_point(aes(color=.data[[fill]]), size=pt_size)
  } else if(border == 'black'){
    p <- p + geom_point(aes(fill=.data[[fill]]), size=pt_size, shape=21, stroke=1)
  } else {
    if(!is.null(fill)){
      p <- p + geom_point(aes(fill=.data[[fill]], color=.data[[border]]), size=pt_size, shape=21, stroke=1)
    } else {
      p <- p + geom_point(aes(color=.data[[border]]), size=pt_size, shape=21, stroke=1)  
    }
  }
  
  if (alpha_fill){
    p <- p + aes(alpha = .data[[fill]])
  } else if (forced_alpha != 'none'){
    p <- p + aes(alpha = forced_alpha)
  }
  
  p <- p + xlim(xlim_min, xlim_max) + ylim(ylim_max, ylim_min)
  
  return(p)
}

SpatialPlotMJP_v2 <- function(seurat_obj, fill=NULL, 
                              xlim_min = NULL, xlim_max = NULL, ylim_min = NULL, ylim_max = NULL){
  
  df <- as.data.frame(seurat_obj@images$slice1@coordinates)
  df[[fill]] <- GetAssayData(seurat_obj, assay = 'SCT')[fill,]
  df[['Manual']] <- seurat_obj@meta.data[['Manual']]  
  
  if(is.null(xlim_min)){
    xlim_min <- seurat_obj@misc$x_alignment_min
  }
  if(is.null(xlim_max)){
    xlim_max <- seurat_obj@misc$x_alignment_max
  }
  if(is.null(ylim_min)){
    ylim_min <- seurat_obj@misc$y_alignment_min
  }
  if(is.null(ylim_max)){
    ylim_max <- seurat_obj@misc$y_alignment_max
  }
  
  p <- ggplot(df, aes(x=imagecol, y=imagerow)) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
    ggpubr::background_image(seurat_obj@misc$image) 
  
  p <- p + geom_point(aes(color = .data[[fill]], alpha = .data[[fill]]), size=2) + viridis::scale_color_viridis() +
    ggnewscale::new_scale_color()
  p <- p + geom_point(data=df, aes(color = .data[['Manual']]), shape=21, size=2, stroke=1) + 
    scale_color_manual(values=c('white', 'black'))
  
  p <- p + xlim(xlim_min, xlim_max) + ylim(ylim_max, ylim_min)
  
  return(p)
}
