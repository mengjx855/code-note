# Encoding: utf-8
# Creator：Jinxin Meng
# Created date：2022-5-22
# Modified date：2023-12-28

pacman::p_load(ggpubr, ggExtra)

# 相关性分析散点图
data("mtcars")
head(mtcars) %>% select(wt, mpg)
# wt  mpg
# 2.620 21.0
# 2.875 21.0
# 2.320 22.8
ggscatter(data = mtcars, x = "wt", y = "mpg")
ggscatter(data = mtcars, x = "wt", y = "mpg", add = "reg.line", conf.int = T, conf.int.level = .95,
          add.params = list(color = "blue", fill = "lightgray", size = 1))

ggscatter(data = mtcars, x = "wt", y = "mpg", add = "reg.line", conf.int = T, conf.int.level = .95,
          add.params = list(color = "blue", fill = "lightgray", size = 1)) +
  stat_cor(label.sep = "\n", color = "black", label.x = 4, label.y = 30)

p <- ggscatter(data = mtcars, x = "wt", y = "mpg", add = "reg.line", conf.int = T, conf.int.level = .95,
          add.params = list(color = "blue", fill = "lightgray", size = 1)) +
  stat_cor(label.sep = "\n", color = "black", label.x = 4, label.y = 30)
ggMarginal(p, type = "histogram", xparams = list(fill = "#377Eb8"), yparams = list(fill = "#E41A1C"))

# 相关性分析散点图，添加拟合值
library(ggpmisc)
ggscatter(data = mtcars, x = "wt", y = "mpg", add = "reg.line",  conf.int = T, conf.int.level = .95) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(adj.rr.label), after_stat(p.value.label), sep = "~~")), formula = y ~ x, parse = T, size = 2.5)

# 相关性分析psych
# profile行为样本，列为观测值
corr <- psych::corr.test(x = profile_x, y = profile_y, method = "spearman", adjust = "BH")
corr_res <- corr_process(corr, rho = 0.5, padj = 0.05, out_df = T)

# 简单的相关性网络图
edges <- select(plotdat, name_x, name_y, n) %>%
  rename(from = name_x, to = name_y)
nodes <-rbind(data.frame(name = unique(edges$from)) %>% add_column(type = "ARG") %>% 
                mutate(class = metadata_args$drug3[match(name, metadata_args$ARO_short_name)],
                       class = forcats::fct_lump_n(class, 9, ties.method = "last")),
              data.frame(name = unique(edges$to)) %>% add_column(type = "MGE") %>% 
                mutate(class = metadata_mges$type[match(name, metadata_mges$MGE_name)],
                       class = forcats::fct_lump_n(class, 9, ties.method = "last"))) %>% 
  mutate(type = factor(type, labels = c("ARG", "MGE")))
graph <- tbl_graph(nodes = nodes, edges = edges)
nodes$degree <- degree(graph)                  #边的数量给点赋一个值
graph <- tbl_graph(nodes = nodes, edges = edges)
class_color <- c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f", "#fb8072","#80b1d3","#bc80bd","#e5c494","#b3b3b3", "#1b9e77","#d95f02","#7570b3","#66a61e","#e6ab02")
setseed(2023)
ggraph(graph, layout = 'fr') + 
  geom_edge_arc(aes(edge_width = n), color = "grey70", strength = .15, show.legend = T) + 
  scale_edge_width(range = c(.2, 1)) +
  scale_edge_alpha_manual(values = c(1, 2)) +
  geom_node_point(aes(fill = class, size = degree), shape = 21, stroke = .3) +
  scale_size(range = c(2, 4)) +
  scale_fill_manual(values = class_color) +
  geom_node_text(aes(label = name), size = .6, fontface = "italic") +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        legend.key = element_rect(fill = NA),
        aspect.ratio = 1)