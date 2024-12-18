#### Jinxin Meng, 20241204, 20241204, v0.1 ####

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(Rtsne)

# t-SNE，全称为 t-distributed Stochastic Neighbor Embedding，翻译为 t分布-随机邻近嵌入。
# 是通过将数据点之间的相似度转化为条件概率，原始空间中数据点的相似度由高斯联合分布表示，
# 嵌入空间中数据点的相似度由学生t分布表示 能够将高维空间中的数据映射到低维空间中，并保留数据集的局部特性。
# t-SNE本质是一种嵌入模型，主要用于高维数据的降维和可视化。
# 如果想象在一个三维的球里面有均匀分布的点，不难想象，如果把这些点投影到一个二维的圆上一定会有很多点是重合的。
# 所以，为了在二维的圆上想尽可能表达出三维里的点的信息，大神Hinton采取的方法：
# 把由于投影所重合的点用不同的距离（差别很小）表示。
# 这样就会占用原来在那些距离上的点，原来那些点会被赶到更远一点的地方。
# t分布是长尾的，意味着距离更远的点依然能给出和高斯分布下距离小的点相同的概率值。
# 从而达到高维空间和低维空间对应的点概率相同的目的。
profile <- read.delim("test/microeco/taxa_abund/mpa_abund.tsv", row.names = 1)
group <- read.delim("test/microeco/sample_table.tsv") %>%
  select(sample=  ID, group = Group)


#### plot_tsne ####
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
plot_tsne <- function(profile, group, group_colnames = NULL, group_order = NULL, group_color = NULL, 
                      display_type = "line", seed = 2024, theta = .5, perplexity = NULL,
                      verbose = F, title = NULL, ellipse_level = .75, aspect_ratio = 3/4,
                      theme = "default") {
  
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  if(is.null(group_order)) group_order <- unique(group$group)
  if(is.null(title)) title <- 't-distributed Stochastic Neighbor Embedding Analysis'
  if(is.null(perplexity)) perplexity <- floor(ncol(profile) - 1) / 3
  
  if(is.null(group_color)){
    color <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3")
    group_color <- rep(color, times = ceiling(length(group_order)/length(color)))[seq_len(length(group_order))]
  }
  
  set.seed(seed)
  tsne_out <- Rtsne::Rtsne(data.frame(t(profile), check.names = F), dims = 2, pca = T, 
                           max_iter = 1000, theta = theta, perplexity = perplexity, verbose = verbose)
  # data 用于降维的原始数据，其中行代表特征，列代表样本，与我们生信分析中常用到的表达谱矩阵相反，需要利用 t() 函数进行转置。
  # dims 降维后的维度数，默认为2，这样降维后的数据可以用平面直角坐标系的散点图进行表示，如果设置为3，则会得到一个3维的降维结果。
  # pca 逻辑型变量，规定是否在t-SNE前预先进行PCA分析，默认为True。
  # max_iter 迭代次数，默认为1000，过小则结果未完全收敛，过大则浪费计算量。
  # theta 计算速度与精确度之间的权衡，范围在0~1之间，越接近0越精确，默认0.5。这个参数影响最终计算结果，可以根据图像效果进行选取。
  # perplexity 这个值是一个正整数，且满足 3*perplexity < nrow(data) - 1 。作为计算数据点相似度的参数， perplexity 可以简单理解为对每个点具有的近邻数量的猜测，代表了平衡数据的局部和全局方面之间的程度，对生成的图像有复杂的影响。一般来说，随着 perplexity 值的增加，形状越来越清晰。
  # verbose 是否输出计算进度，在数据集较大的时候比较实用。
  
  plot_data <- tsne_out$Y %>% 
    data.frame() %>% 
    dplyr::rename_with(~ c('X1', 'X2')) %>% 
    add_column(sample = colnames(profile), .before = 1) %>% 
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
      labs(x = 'tSNE_1', y = 'tSNE_2', title = title, color = "group")
    
  } else if (display_type == "dot") {
    
    p <- ggplot(plot_data, aes(x = X1, y = X2, color = group)) +
      geom_vline(xintercept = 0, lty = 2, lwd = .4) +
      geom_hline(yintercept = 0, lty = 2, lwd = .4) +
      geom_point(size = 1.5) +
      stat_ellipse(aes(fill = group), geom = 'polygon', level = ellipse_level, 
                   alpha = .05, lty = 2, lwd = .3, show.legend = F) +
      scale_color_manual(values = group_color) +
      scale_fill_manual(values = group_color) +
      labs(x = 'tSNE_1', y = 'tSNE_2', title = title, color = "group")
    
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
  
  cat("  ggsave(file = \"tSNE.pdf\", width = 6, height = 4.5)\n")
  return(p)
}
