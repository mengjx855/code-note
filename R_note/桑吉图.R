# Encoding: utf-8
# Creator：Jinxin Meng
# Created Date:2023-6-19
# Modified Date:2023-12-15

devtools::install_github("davidsjoberg/ggsankey")
pacman::p_load(tidyverse, ggsankey)

# 测试
df <- mtcars %>% make_long(cyl, vs, am, gear, carb)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node,fill = factor(node))) +
  geom_sankey()

# 准备一个文件，列出所有对应关系，不需要统计。例如：
# phylum            species         drug            mechanism
# p__Pseudomonadota s__Escherichia  Fluoroquinolone Efflux
# p__Pseudomonadota s__Enterobacter Fluoroquinolone Efflux
# p__Pseudomonadota s__Hafnia       Fluoroquinolone Efflux
# p__Pseudomonadota s__Escherichia  Fluoroquinolone Efflux

# 数据转换为ggsankey的格式
plotdat <- make_long(select(dat, -1), phylum, species, drug, mechanism)
# 数据是四列，第一列列出了x轴，表明共有几组数据进行流动 （这里是三组数据，分别是phylum, species, drug, mechanism）
# 第二列是node，每组数据的每种元素均会分配一个唯一的node（这里的node就是phylum, species, drug, mechanism数据集合中所有元素种类）
# 第三列是x轴的标签
# 第四列是node的标签

# ggsankey中的图层和元素的划分：
# node指的是每组数据每种元素在途中在图中呈现的柱子；
# flow指的是数据和数据之间的流动图形；
# stage指的是每组数据

# geom_sankey属性:
# geom_sankey()中 flow.xxx (flow.fill, flow.alpha ...) 设置的是流动部分的属性；
# geom_sankey()中 node.xxx (xxx.fill, node.alpha ...) 设置的是每组数据柱子的属性；
# geom_sankey()中 width设置的是柱子的宽度，space设置的每组数据不同属性柱子之间的距离；
# geom_sankey()中 smooth设置的是流动部分的平滑程度；
# geom_sankey()中 fill、color和lwd等参数设置的是柱子和流动图形的属性，前提是没用flow.xxx参数指定flow图形的属性，那么fill、color和lwd只对柱子起作用；
# geom_sankey()中 space设置不同柱子之间的距离
# 每组数据下显示的名字的内容可以在scale_x_discrete()中用lable参数设置，格式可以在theme()中用axis.text.x设置

# theme_sankey() 是sankey图专用的主题属性函数
# geom_sankey_text()和geom_sankey_label() 是sankey图专用的设置标签属性的函数
# geom_sankey_text()中 space设置与柱子同步
ggplot(plotdat, aes(x = x, next_x = next_x, node = node, next_node = next_node, 
                    fill = node, label = node)) + # 这是固定的（x = x, next_x = next_x, node = node, next_node = next_node）
  geom_sankey(node.color = "black", flow.fill = "grey60", flow.alpha = .6, lwd = .3,
              width = .1, space = 300, smooth = 8, show.legend = F) +
  geom_sankey_text(size = 3, color = "black", space = 300) +
  scale_x_discrete(label = c("Cut off", "Drug class", "antibiotic mechanism")) +
  labs(x = "") +
  theme_sankey(base_size = 10) +
  theme(axis.text.x = element_text(color = "black", size = 8))
ggsave("ARG_sankey.pdf", width = 6, height = 4.5)

# 这种画出来的图是按照名字分配顺序的，如何按照柱子的高度重新排序呢？需要对plotdat中的node重新编顺序，用因子类型，例如，对每列进行排序，每列的排序结果合起来，对总体node进行排序
phylum_order <- group_by(dat, phylum) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(phylum) %>% unlist() %>% as.character()
species_order <- group_by(dat, species) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(species) %>% unlist() %>% as.character()
drug_order <- group_by(dat, drug) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(drug) %>% unlist %>% as.character()
mechanism_order <- group_by(dat, mechanism) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(mechanism) %>% unlist %>% as.character()
# 然后对node重排序，这样得到的图就是有顺序的了。
plotdat <- make_long(select(dat, -1), phylum, species, drug, mechanism) %>% 
  mutate(node = as.factor(node), node = forcats::fct_relevel(node, rev(c(phylum_order, species_order, drug_order, mechanism_order))))


