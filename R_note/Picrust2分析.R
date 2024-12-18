#####Info####
# Encoding: utf-8
# Creator: network
# Modifier: 小梦、李铭晗
# Created Date: 2021-10-15
# Modified Date：2022-03-13

####option####
setwd("F:/Proj_Other/07.chicken_16s_wangxy/Picrust2/")
pacman::p_load(tidyverse,dplyr,tibble,tidyr,ggplot2,patchwork)

####metabolism KEGG####
load("F:/Database/KEGG/kegg_pathway_database-2021.Rdata")
kegg_database_metabolism <- kegg_database %>% subset(level1des%in%"Metabolism")
group <- read.delim("../Data/sample_group.txt", row.names = NULL, sep = '\t', check.names = F) # 
profile <- read.delim("pred_metagenome_unstrat.tsv", sep = "\t", row.names = 1) %>% select(group$sample)
profile$level2 <- kegg_database_metabolism$level2[match(rownames(profile), kegg_database_metabolism$level4)]
profile <- profile[!is.na(profile$level2),]
profile <- profile %>% group_by(level2) %>% summarise_all(sum) %>% column_to_rownames(var = "level2")
profile_relative <- data.frame(apply(profile, 2, function(x) x/sum(x)) * 100)  # 计算相对丰度
group_pair <- data.frame(pair1 = c("HP","YP","HP","HE"), pair2 = c("HE","YE","YP","YE"))  # 设置分组对
error_pair <- c()
abun_bar_dt <- rbind()
diff_mean_dt <- rbind()
for (j in 1:nrow(group_pair)) {  # 生成第i对比较数据集
  pair <- unlist(group_pair[j,])
  filename <- paste0(pair[1], "_vs_", pair[2])
  sample_list_j <- group$sample[group$group%in%pair]
  group_list_j <- group$group[group$group%in%pair]
  dt_j <- data.frame(t(profile_relative[,sample_list_j])) %>% filter(apply(., 1, mean) > 1)  # 代谢平均丰度要大于1% 
  dt_j <- cbind(dt_j, group = as.factor(group_list_j))
  diff <- dt_j %>% select_if(is.numeric) %>% map_df(~ broom::tidy(t.test(. ~ group, data = dt_j)), .id = 'var')   # t-test
  # diff$p.value <- p.adjust(diff$p.value,"bonferroni")  # FDR 
  diff <- diff %>% filter(p.value < 0.05)  # 过滤
  
  # diff <- dt_j %>% select_if(is.numeric) %>% map_df(~ broom::tidy(wilcox.test(. ~ group, data = dt_j)), .id = 'var')   # wilcox-test
  # 
  # diff$p.value <- p.adjust(diff$p.value,"bonferroni")  # FDR
  # 
  # diff <- diff %>% filter(p.value < 0.05)  # 过滤
  
  if (nrow(diff) == 0) {  # 判断是否有差异结果，没有结果则保存item
    error_pair <- c(error_pair, filename)
    next
  } else {  # 生成绘图数据
    dt_level2 <- unique(data.frame(kegg_database_metabolism$level2,kegg_database_metabolism$level2des))  # 提取level2编号和描述信息
    abun_bar <- dt_j[,c(diff$var,"group")] %>% gather(key = variable, value = value, -group) %>% group_by(variable, group) %>% summarise_at(.vars = "value", list(Mean = mean))  # 左侧条形囿数据
    abun_bar$variable <- dt_level2$kegg_database_metabolism.level2des[match(abun_bar$variable,dt_level2$kegg_database_metabolism.level2)]  # 左侧条形囿数据
    abun_bar_j <- data.frame(abun_bar)
    abun_bar_j$class <- rep(filename, each = nrow(abun_bar_j))
    abun_bar_dt <- rbind(abun_bar_dt, abun_bar_j)
    diff_mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]  # 右侧散点囿数据
    diff_mean$group <- c(ifelse(diff_mean$estimate > 0, levels(dt_j$group)[1], levels(dt_j$group)[2]))
    diff_mean <- diff_mean[order(diff_mean$estimate, decreasing = T),]
    diff_mean$var <- dt_level2$kegg_database_metabolism.level2des[match(diff_mean$var,dt_level2$kegg_database_metabolism.level2)]
    diff_mean_j <- data.frame(diff_mean)   # 右侧图保存数据
    diff_mean_j$class <- rep(filename, each = nrow(diff_mean_j))
    diff_mean_dt <- rbind(diff_mean_dt, diff_mean_j)
    # 左侧条形囿
    cbbPalette <- c("#E69F00", "#56B4E9")
    abun_bar$variable <- factor(abun_bar$variable, levels = rev(diff_mean$var))
    p1 <- ggplot(abun_bar, aes(variable, Mean, fill = group)) +
      scale_x_discrete(limits = levels(diff_mean$var), labels = function(x) str_wrap(x, width = 20)) +
      labs(x = "", y = "Mean proportion (%)") +
      coord_flip() + 
      theme(panel.background = element_rect(fill = 'transparent'),
            panel.grid = element_blank(),
            axis.ticks.length = unit(0.4, "lines"),
            axis.ticks = element_line(color = 'black'),
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(colour = 'black', size = 12, face = "bold"),
            axis.text = element_text(colour = 'black',size = 12, face = "bold"),
            legend.title = element_blank(),
            legend.text = element_text(size = 10, face = "bold",colour = "black", margin = margin(r = 20)),
            legend.position = c(-0.12,-0.07),
            legend.background = element_rect(fill = "transparent"),
            legend.direction = "horizontal",
            legend.key.width = unit(0.8,"cm"),
            legend.key.height = unit(0.5,"cm"))
    if (nrow(diff_mean) > 1) {
      for (i in 1:(nrow(diff_mean)-1)) {
      p1 <- p1 + 
        annotate('rect', xmin = i + 0.5, xmax = i + 1.5, ymin = -Inf, ymax = Inf, fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
      p1 <- p1 + 
        geom_bar(stat = "identity", position = "dodge", width = 0.7, colour = "black") +
        scale_fill_manual(values = cbbPalette)
      }
    } else { 
      p1 <- p1 + 
        annotate('rect', xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = 'white')
      p1 <- p1 + 
        geom_bar(stat = "identity", position = "dodge", width = 0.7, colour = "black") +
        scale_fill_manual(values = cbbPalette)
      }
    # 右侧散点囿
    diff_mean$var <- factor(diff_mean$var, levels = levels(abun_bar$variable))
    diff_mean$p.value <- signif(diff_mean$p.value, 3)
    diff_mean$p.value <- as.character(diff_mean$p.value)
    p2 <- ggplot(diff_mean,aes(var,estimate,fill = group)) +
      scale_x_discrete(limits = levels(diff_mean$var)) +
      labs(x = "", y = "Difference in mean proportions (%)", title = "95% confidence intervals") +
      coord_flip() +
      theme(panel.background = element_rect(fill = 'transparent'),
            panel.grid = element_blank(),
            axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(colour = 'black', size = 12, face = "bold"),
            axis.text = element_text(colour = 'black', size = 12, face = "bold"),
            axis.text.y = element_blank(),
            legend.position = "none",
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(size = 13, face = "bold", colour = "black", hjust = 0.5))
    if (nrow(diff_mean) > 1) {
    for (i in 1:(nrow(diff_mean) - 1)) {
      p2 <- p2 + 
        annotate('rect', xmin = i + 0.5, xmax = i + 1.5, ymin = -Inf, ymax = Inf, fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
      p2 <- p2 + 
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(0.8), width = 0.3, size = 0.3) +
        geom_point(shape = 21, size = 5) +
        geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black') +
        scale_fill_manual(values = cbbPalette)
      p3 <- ggplot(diff_mean, aes(var, estimate, fill = group)) +
        geom_text(aes(y = 0, x = var), label = diff_mean$p.value, hjust = 0, fontface = "bold", inherit.aes = F, size = 4) +
        geom_text(aes(x = nrow(diff_mean)/2 + 0.5, y = 0.3), label = "P-value (corrected)", srt = 90, fontface = "bold", size = 5) +
        coord_flip() +
        ylim(c(0,0.5)) +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              plot.margin = unit(c(0,7,0,0),"mm"))
      }
    } else { 
      p2 <- p2 + 
        annotate('rect', xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = 'white')
      p2 <- p2 + 
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(0.8), width = 0.3, size = 0.3) +
        geom_point(shape = 21, size = 5) +
        geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black') +
        scale_fill_manual(values = cbbPalette)
      p3 <- ggplot(diff_mean, aes(var, estimate, fill = group)) +
        geom_text(aes(y = 0, x = var), label = diff_mean$p.value, hjust = 0, fontface = "bold", inherit.aes = F, size = 4) +
        geom_text(aes(x = 1, y = 0.3), label = "P-value (corrected)", srt = 90, fontface = "bold", size = 4) +
        coord_flip() +
        ylim(c(0,0.5)) +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              plot.margin = unit(c(0,7,0,0),"mm"))
    }
    if (nrow(diff_mean) > 1) {
      p <- p1 + p2 + p3 + plot_layout(widths = c(6, 5, 2))  # 图像拼接
      height <- length(diff_mean$var) * 2 * 0.75
      ggsave(paste0(filename,".pdf"), p, width = 13, height = height)  # 保存图像
      } else { 
      p <- p1 + p2 + p3 + plot_layout(widths = c(6, 5, 2))  # 图像拼接
      ggsave(paste0(filename,".pdf"), p, width = 13, height = 2.5)  # 保存图像
      }
    }
  write.table(data.frame(error_pair), "dt_Error_pair.txt", sep = "\t", row.names = F, quote = F)
  write.table(abun_bar_dt, "dt_abun_bar.txt", sep = "\t", row.names = F, quote = F)
  write.table(diff_mean_dt,"dt_diff_mean.txt", sep = "\t", row.names = F, quote = F)
  }

# All KEGG ----------------------------------------------------------------
load("F:/Database/KEGG/kegg_pathway_database-2021.Rdata")
group <- read.delim("../Data/sample_group.txt", row.names = NULL, sep = '\t', check.names = F)
profile <- read.delim("pred_metagenome_unstrat.tsv", sep = "\t", row.names = 1) %>% select(group$sample)
profile$level2 <- kegg_database$level2[match(rownames(profile), kegg_database$level4)]
profile <- profile[!is.na(profile$level2),]
profile <- profile %>% group_by(level2) %>% summarise_all(sum) %>% column_to_rownames(var = "level2")
profile_relative <- data.frame(apply(profile, 2, function(x) x/sum(x)) * 100)  # 计算相对丰度
group_pair <- data.frame(pair1 = c("HP","YP","HP","HE"), pair2 = c("HE","YE","YP","YE"))  # 设置分组对
error_pair <- c()
abun_bar_dt <- rbind()
diff_mean_dt <- rbind()
for (j in 1:nrow(group_pair)) {  # 生成第j对比较数据集
  pair <- unlist(group_pair[j,])
  filename <- paste0(pair[1], "_vs_", pair[2])
  sample_list_j <- group$sample[group$group%in%pair]
  group_list_j <- group$group[group$group%in%pair]
  
  dt_j <- data.frame(t(profile_relative[,sample_list_j])) %>% filter(apply(., 1, mean) > 1)  # 代谢平均丰度要大于1% 
  dt_j <- cbind(dt_j, group = as.factor(group_list_j))
  
  diff <- dt_j %>% select_if(is.numeric) %>% map_df(~ broom::tidy(t.test(. ~ group, data = dt_j)), .id = 'var')   # t-test
  # diff$p.value <- p.adjust(diff$p.value,"bonferroni")  # FDR
  diff <- diff %>% filter(p.value < 0.05)  # 过滤
  
  # diff <- dt_j %>% select_if(is.numeric) %>% map_df(~ broom::tidy(wilcox.test(. ~ group, data = dt_j, exact = F)), .id = 'var')   # wilcox-test
  # diff$p.value <- p.adjust(diff$p.value,"bonferroni")  # FDR
  # diff <- diff %>% filter(p.value < 0.05)  # 过滤
  
  if (nrow(diff) == 0) {  # 判断是否有差异结果，没有结果则保存item
    error_pair <- c(error_pair, filename)
    next
  } else {  # 生成绘图数据
    dt_level2 <- unique(data.frame(kegg_database$level2,kegg_database$level2des))  # 提取level2编号和描述信息
    abun_bar <- dt_j[,c(diff$var,"group")] %>% gather(key = variable, value = value, -group) %>% group_by(variable, group) %>% summarise_at(.vars = "value", list(Mean = mean))  # 左侧条形囿数据
    abun_bar$variable <- dt_level2$kegg_database.level2des[match(abun_bar$variable,dt_level2$kegg_database.level2)]  # 左侧条形囿数据
    abun_bar_j <- data.frame(abun_bar)
    abun_bar_j$class <- rep(filename, each = nrow(abun_bar_j))
    abun_bar_dt <- rbind(abun_bar_dt, abun_bar_j)
    diff_mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]  # 右侧散点囿数据
    diff_mean$group <- c(ifelse(diff_mean$estimate > 0, levels(dt_j$group)[1], levels(dt_j$group)[2]))
    diff_mean <- diff_mean[order(diff_mean$estimate, decreasing = T),]
    diff_mean$var <- dt_level2$kegg_database.level2des[match(diff_mean$var,dt_level2$kegg_database.level2)]
    diff_mean_j <- data.frame(diff_mean)   # 右侧图保存数据
    diff_mean_j$class <- rep(filename, each = nrow(diff_mean_j))
    diff_mean_dt <- rbind(diff_mean_dt, diff_mean_j)
    
    # 左侧条形囿
    cbbPalette <- c("#E69F00", "#56B4E9")
    abun_bar$variable <- factor(abun_bar$variable, levels = rev(diff_mean$var))
    p1 <- ggplot(abun_bar, aes(variable, Mean, fill = group)) +
      scale_x_discrete(limits = levels(diff_mean$var), labels = function(x) str_wrap(x, width = 20)) +
      labs(x = "", y = "Mean proportion (%)") +
      coord_flip() + 
      theme(panel.background = element_rect(fill = 'transparent'),
            panel.grid = element_blank(),
            axis.ticks.length = unit(0.4, "lines"),
            axis.ticks = element_line(color = 'black'),
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(colour = 'black', size = 12, face = "bold"),
            axis.text = element_text(colour = 'black',size = 12, face = "bold"),
            legend.title = element_blank(),
            legend.text = element_text(size = 10, face = "bold",colour = "black", margin = margin(r = 20)),
            legend.position = c(-0.12, -0.07),
            legend.background = element_rect(fill = "transparent"),
            legend.direction = "horizontal",
            legend.key.width = unit(0.8, "cm"),
            legend.key.height = unit(0.5, "cm"))
    
    if (nrow(diff_mean) > 1) {
      for (i in 1:(nrow(diff_mean)-1)) {
        p1 <- p1 + annotate('rect', xmin = i + 0.5, xmax = i + 1.5, ymin = -Inf, ymax = Inf, fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
        p1 <- p1 + geom_bar(stat = "identity", position = "dodge", width = 0.7, colour = "black") +
          scale_fill_manual(values = cbbPalette)
      }
    } else { 
      p1 <- p1 + 
        annotate('rect', xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = 'white') 
      p1 <- p1 + 
        geom_bar(stat = "identity", position = "dodge", width = 0.7, colour = "black") +
        scale_fill_manual(values = cbbPalette)
    }
    
    # 右侧散点囿
    diff_mean$var <- factor(diff_mean$var, levels = levels(abun_bar$variable))
    diff_mean$p.value <- signif(diff_mean$p.value, 3)
    diff_mean$p.value <- as.character(diff_mean$p.value)
    p2 <- ggplot(diff_mean,aes(var,estimate,fill = group)) +
      scale_x_discrete(limits = levels(diff_mean$var)) +
      labs(x = "", y = "Difference in mean proportions (%)", title = "95% confidence intervals") +
      coord_flip() +
      theme(panel.background = element_rect(fill = 'transparent'),
            panel.grid = element_blank(),
            axis.ticks.length = unit(0.4,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(colour = 'black', size = 12, face = "bold"),
            axis.text = element_text(colour = 'black', size = 12, face = "bold"),
            axis.text.y = element_blank(),
            legend.position = "none",
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(size = 13, face = "bold", colour = "black", hjust = 0.5))
    
    if (nrow(diff_mean) > 1) {
      for (i in 1:(nrow(diff_mean) - 1)) {
        p2 <- p2 + 
          annotate('rect', xmin = i + 0.5, xmax = i + 1.5, ymin = -Inf, ymax = Inf, fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
        p2 <- p2 + 
          geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(0.8), width = 0.4, size = 0.4) +
          geom_point(shape = 21, size = 5) +
          geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black') +
          scale_fill_manual(values = cbbPalette)
        
        p3 <- ggplot(diff_mean, aes(var, estimate, fill = group)) +
          geom_text(aes(y = 0, x = var), label = diff_mean$p.value, hjust = 0, fontface = "bold", inherit.aes = F, size = 4) +
          geom_text(aes(x = nrow(diff_mean)/2 + 0.5, y = 0.3), label = "P-value (corrected)", srt = 90, fontface = "bold", size = 5) +
          coord_flip() +
          ylim(c(0,0.5)) +
          theme(panel.background = element_blank(),
                panel.grid = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                plot.margin = unit(c(0,7,0,0),"mm"))
      }
    } else { 
      p2 <- p2 + 
        annotate('rect', xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = 'white')
      p2 <- p2 + 
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(0.8), width = 0.4, size = 0.4) +
        geom_point(shape = 21, size = 5) +
        geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black') +
        scale_fill_manual(values = cbbPalette)
        p3 <- ggplot(diff_mean, aes(var, estimate, fill = group)) +
        geom_text(aes(y = 0, x = var), label = diff_mean$p.value, hjust = 0, fontface = "bold", inherit.aes = F, size = 4) +
        geom_text(aes(x = 1, y = 0.3), label = "P-value (corrected)", srt = 90, fontface = "bold", size = 4) +
        coord_flip() +
        ylim(c(0,0.5)) +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              plot.margin = unit(c(0,7,0,0),"mm"))
    }
    if (nrow(diff_mean) > 1) {
      p <- p1 + p2 + p3 + plot_layout(widths = c(6, 5, 2))  # 图像拼接
      height <- length(diff_mean$var) * 2 * 0.75
      ggsave(paste0(filename,".pdf"), p, width = 13, height = height)  # 保存图像
    } else { 
      p <- p1 + p2 + p3 + plot_layout(widths = c(6, 5, 2))  # 图像拼接
      ggsave(paste0(filename,".pdf"), p, width = 13, height = 2.5)  # 保存图像
    }
  }
  write.table(data.frame(error_pair), "dt_Error_pair.txt", sep = "\t", row.names = F, quote = F)
  write.table(abun_bar_dt, "dt_abun_bar.txt", sep = "\t", row.names = F, quote = F)
  write.table(diff_mean_dt,"dt_diff_mean.txt", sep = "\t", row.names = F, quote = F)
  
}

