####info####
# Encoding: utf-8
# Modified date: 2023-1-5

####option####
pacman::p_load(dplyr,tibble,tidyr,stringr,ggplot2,ggpubr)
options(stringsAsFactors = F, digits = 10, scipen = 20)
setwd("F:/Proj/")

####plot####
dt <- read.delim("viral_attr_plot_dat.txt", sep = "\t", header = T, row.names = NULL)

# 构造饼图数据
# 饼图实际上就是堆叠柱状图转变为极坐标之后出现的，默认的顺序是顺时针，可以调整，所有计算位置的时候要注意，
# perc和ypos的顺序是一致的，填充要是位置不对，那么就在标签位置的位置用1-ypos得到准确的位置
# 可以在堆叠图画完之后先看看，并且把标记表上在转变为饼图
# A tibble: 5 × 7
# Groups:   direction [1]
#  direction taxa         count   perc label                   ypos angle
#  <chr>     <chr>        <int>  <dbl> <chr>                  <dbl> <dbl>
# 1 Case      Microviridae     1 0.0244 Microviridae (n = 1)  0.0122  274.
# 2 Case      Myoviridae       1 0.0244 Myoviridae (n = 1)    0.0366  283.
# 3 Case      Retroviridae     5 0.122  Retroviridae (n = 5)  0.110   310.
# 4 Case      Siphoviridae    10 0.244  Siphoviridae (n = 10) 0.293   375.
# 5 Case      Unclassified    24 0.585  Unclassified (n = 24) 0.707   345.
dt <- dt %>% 
  select(checkv_quality) %>% 
  table(.) %>% 
  data.frame() %>%
  rename(., Tp = .) %>% 
  mutate(Perc = Freq/sum(Freq), # 计算比例
         label = paste0(Tp, " ", round((Perc)*100, digits = 1), "%"), # 标记在饼图中的标签
         ypos = cumsum(Perc) - 0.5 * Perc, # 计算标记的位置，原理上是：标签是在每个块的中间，累加值减去最后一次值的一般
         angle = 360 * ypos + 90, # 计算标签的角度，使用ypos乘以360得到每个位置的角度，然后+90形成发散的状态
         angle = ifelse(ypos < .5, angle + 180, angle)) # 左半部分的标签要朝下，这样符合审美

# 表中如果没有计算angle，可在text部分的aes中设置也是可以的，啥原理还没弄明白
angle = atan(y/x)*360/(2*pi)

# 可视化基础代码
ggplot(dt, aes(x = "", y = Perc , fill = Tp)) +
  geom_bar(stat = "identity", width = 1 , color = "white", show.legend = F) +
  geom_text(data = dt, aes(x = 1.4, y = 1 - ypos, label = label, angle = angle), size = 3) + 
  scale_fill_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3")) +
  coord_polar("y", start = 0) +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.key = element_rect(fill = 'transparent'), 
        legend.position = "right", 
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

# 饼图分面画法
plot_dat <- tib0 %>% 
  group_by(enrich) %>% 
  mutate(Perc = N/sum(N),
         label = paste0(taxo, ", ", N),
         ypos = cumsum(Perc) - 0.5 * Perc,
         angle = 360 * ypos + 90,
         angle = ifelse(ypos < .5, angle + 180, angle))
# 可视化
ggplot(plot_dat, aes(x = "", y = Perc , fill = taxo)) +
  geom_bar(stat = "identity", width = 1 , color = "white", show.legend = F) +
  geom_text_repel(data = plot_dat, aes(x = 1.5, y = 1 - ypos, label = label), size = 3) + 
  facet_grid(cols = vars(enrich)) +
  coord_polar("y", start = 0) +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.key = element_rect(fill = 'transparent'), 
        legend.position = "right", 
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())