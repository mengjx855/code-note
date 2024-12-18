#### Jinxin Meng, 20241204, 20241204, v0.1 ####

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(umap)

# 对高维单细胞数据的可视化展示，以t-SNE为代表的非线性降维技术，由于其能够避免集群表示的过度拥挤，
# 在重叠区域上能表示出不同的集群而被广泛运用。然而，任何技术方法都不是完美的，t-SNE也一样，
# 它的局限性体现在丢失大规模信息(集群间关系)、计算时间较慢以及无法有效地表示非常大的数据集]等方面。
# UMAP就是这样一个能解决这些问题的降维和可视化的工具。
# UMAP：统一流形逼近与投影(UMAP，Uniform Manifold Approximation and Projection)是一种新的降维流形学习技术。
# UMAP是建立在黎曼几何和代数拓扑理论框架上的。UMAP是一种非常有效的可视化和可伸缩降维算法。
# 在可视化质量方面，UMAP算法与t-SNE具有竞争优势，但是它保留了更多全局结构、具有优越的运行性能、更好的可扩展性。
# 此外，UMAP对嵌入维数没有计算限制，这使得它可以作为机器学习的通用维数约简技术。
# 小数据集中，t-SNE和UMAP差别不是很大
# 大数据集中，UMAP优势明显（30多万个细胞的降维可视化分析）
# 通过数据降维和可视化展示的比较显示，PCA分群效果最差，UMAP和t-SNE都成功将与相似细胞群相对应的簇聚集在一起。
# 与t-SNE相比，UMAP还提供了有用的和直观的特性、保留了更多的全局结构，特别是细胞子集的连续性。

#### plot_umap ####
# group_order: group level
# group_color: group color
# profile
#         s1  s2 ...
# gene1   2   3  
# gene2   12  2  
# ...
# group
# sampple group        
# gene1   ctr
# gene2   case
# ...
plot_umap <- function(profile, group, group_colnames = NULL, group_order = NULL, group_color = NULL, 
                      display_type = "line", title = NULL, ellipse_level = .75, aspect_ratio = 3/4,
                      theme = "default") {
  
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  if(is.null(group_order)) group_order <- unique(group$group)
  if(is.null(title)) title <- 'Uniform Manifold Approximation and Projection Analysis'
  
  if(is.null(group_color)){
    color <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3")
    group_color <- rep(color, times = ceiling(length(group_order)/length(color)))[seq_len(length(group_order))]
  }
  
  umap_out <- umap::umap(data.frame(t(profile), check.names = F))
  
  plot_data <- umap_out$layout %>% 
    data.frame() %>% 
    dplyr::rename_with(~ c('X1', 'X2')) %>% 
    rownames_to_column('sample') %>% 
    left_join(group, by = 'sample') %>% 
    mutate(group = factor(group, levels = group_order))
  
  if (display_type == "line") {
    
    means <- plot_data %>% 
      dplyr::select(-sample) %>% 
      group_by(group) %>% 
      summarise_all(mean) %>% 
      dplyr::rename(X1mean = X1, X2mean = X2) %>% 
      merge(x = plot_data %>% 
              dplyr::select(-sample) %>% 
              relocate(group), 
            y = ., by = "group") %>% 
      mutate(color = group_color[match(group, group_order)])
    
    p <- ggplot(means, aes(x = X1, y = X2, fill = group)) +
      geom_vline(xintercept = 0, color = "gray50", linewidth = .3, lty = "dashed") + 
      geom_hline(yintercept = 0, color = "gray50", linewidth = .3, lty = "dashed")
    
    for (i in 1:nrow(means)) {
      path <-  means[i, 2:5] %>% 
        unlist() %>% 
        as.numeric()
      var_color <-  means[i, 6] %>%
        unlist() %>% 
        as.character()
      p <- p + annotate(geom = "segment", x = path[1], y = path[2], 
                        xend = path[3], yend = path[4], color = var_color, lwd = .4)
    }

    p <- p + geom_point(aes(color = group), size = 1, show.legend = F) +
      geom_point(data = means %>% 
                   dplyr::select(group, X1mean, X2mean) %>% 
                   unique(), 
                 aes(x = X1mean, y = X2mean, color = group), 
                 size = 2, inherit.aes = F) +
      stat_ellipse(aes(color = group, fill = group), geom = 'polygon', 
                   level = ellipse_level, alpha = .05, lty = 2, lwd = .3, show.legend = F) +
      scale_color_manual(values = group_color) +
      scale_fill_manual(values = group_color) +
      labs(x = 'UMAP_1', y = 'UMAP_2', title = title, color = "group")
  
  } else if (display_type == "dot") {
    
    p <- ggplot(plot_data, aes(x = X1, y = X2, color = group)) +
      geom_vline(xintercept = 0, lty = 2, lwd = .4) +
      geom_hline(yintercept = 0, lty = 2, lwd = .4) +
      geom_point(size = 1.5) +
      stat_ellipse(aes(fill = group), geom = 'polygon', level = ellipse_level, 
                   alpha = .05, lty = 2, lwd = .3, show.legend = F) +
      scale_color_manual(values = group_color) +
      scale_fill_manual(values = group_color) +
      labs(x = 'UMAP_1', y = 'UMAP_2', title = title, color = "group")
    
  } else {
    
    stop("ERROR in show_type parameter .. dot or line.") 
  
  }
  
  if (theme == "pubr") {
    
    p <- p + 
      ggpubr::theme_pubr() + 
      theme(aspect.ratio = aspect_ratio, 
            legend.position = 'right')
    
  } else {
    
    p <- p +
      theme_bw() +
      theme(axis.ticks = element_line(linewidth = .5, color = "black"),
            axis.title = element_text(size = 12, color = "black"),
            axis.text = element_text(size = 12, color = "black"),
            axis.line = element_blank(),
            plot.title = element_text(size = 12, color = "black"),
            plot.subtitle = element_text(size = 12, color = "black"),
            panel.background = element_blank(),
            panel.border = element_rect(linewidth = .5, color = "black"),
            panel.grid = element_blank(),
            legend.text = element_text(size = 12, color = "black"),
            legend.title = element_text(size = 12, color = "black"),
            aspect.ratio = aspect_ratio)
  }
  
  cat("  ggsave(file = \"umap.pdf\", width = 6, height = 4.5)\n")
  return(p)
}
