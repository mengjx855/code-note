# Encoding: utf-8
# Creator：Jinxin Meng
# Created Date: 2023-01-15
# Modified Date: 2023-01-15

####ggraph####
####ggraph图形创建####
# https://zhuanlan.zhihu.com/p/375171266
# 主要使用两个包：ggraph和tidygraph
# tidygraph将图对象设计为两个 tidy 表，一个用于表示 node 数据，另一个用于表示 edge 数据。
# tidygraph提供了一些额外的动词用于操作这两个表，同时也提供了很多的图形算法，让图形对象的处理看起来更像是在处理数据框。
# ggraph包含3个核心概念：
# layout：定义图的布局，包含所有的igraph布局以及额外的一些布局。
# nodes：定义节点图形属性，使用geom_node_*()函数来控制
# edges：定义边的图形属性，使用geom_edge_*()函数来控制
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

####创建图####
# 我们使用tbl_graph()函数来创建tbl_graph对象
# as_tbl_graph()函数还可以将很多其他类型的数据转换为tbl_graph对象，例如igraph对象，hclust对象、邻接矩阵、dendrogram对象等等
edges <- data.frame(from = c("v01", "v02", "v02", "v02", "v03", "v04", "v04", "v04", "v05"),
                    to = c("s__Anaerobutyricum_hallii","s__Clostridiales_bacterium","s__Anaerobutyricum_hallii","s__Dorea_sp_AF36_15AT","s__Anaerobutyricum_hallii","s__Clostridiales_bacterium","s__Anaerobutyricum_hallii","s__Dorea_sp_AF36_15AT","s__Clostridiales_bacterium"))
nodes <- data.frame(name = unique(union(edges$from, edges$to)))
graph <- tbl_graph(nodes = nodes, edges = edges)

# 从tbl_graph提取nodes信息
graph %>% activate(nodes) %>% as_tibble()

# 查看类型，tbl_graph对象本质上是一个igraph对象
class(graph)

# tidygraph也提供了很多简便函数，create_*()函数用于创建一些常见的图结构
par(mfrow= c(1, 3))
plot(create_star(10))
plot(create_ring(10))
plot(create_tree(n = 20, children = 3), edge.arrow.size = .4)

# play_*() 函数用于创建模拟图
par(mfrow= c(1, 3))
plot(play_geometry(6, 1, torus = FALSE))
plot(play_islands(4, 10, 0.7, 3))
plot(play_forestfire(20, 0.5), edge.arrow.size = .2)

####布局####
# 使用 ggraph 来绘制图形
# ggraph() 函数相当于 ggplot2::ggplot()，根据传入的图对象以及布局来创建绘图对象
# ggraph 默认会根据图结构自动推断布局，也可以使用指定 layout 参数的值
# 如果布局算法可以接受额外的参数，也可以在ggraph()函数中一并指定
ggraph(graph) +
  geom_edge_link() +
  geom_node_point() +
  theme_graph()

ggraph(graph, layout = 'kk') +
  geom_edge_link() +
  geom_node_point() +
  theme_graph()
ggraph(graph, layout = 'kk', maxiter = 100) +
  geom_edge_link() +
  geom_node_point() +
  theme_graph()

# 或者用create_layout()函数来创建布局，它接受的参数与ggraph()一样，但是返回的是layout_ggraph对象，可以在后续的图结构中使用
# 返回的对象是包含节点位置及属性信息的数据框
layout <- create_layout(graph, layout = 'eigen')
ggraph(layout) +
  geom_edge_link() +
  geom_node_point() +
  theme_graph()
head(layout)

# 查看对象的所有属性
attributes(layout)


####ggraph图形布局####
# ggraph 包提供了许多不同的布局，包括 igraph 所有的布局，同时也提供了一些它自己设计的布局，有超过 20 种布局可供选择。
# 通常来说，不同的布局算法对图的展示差异很大，因此，为自己的数据找到一个合适的布局很重要。
# 也可以自己设计一个布局函数，接受一个 tbl_graph 对象，并返回一个位置数据框。或者直接提供一个位置矩阵或数据框作为布局

# 一些布局既可以显示在笛卡尔坐标系中，也可以在极坐标中有效地进行展示。对于 ggplot2 来说，可以使用 coord_polar() 转换为极坐标轴。但这并不适用于 ggraph，我们通过 circular 参数将布局转换为径向表示。
# 我们先为边添加一列信息，用于表示相关性，-1 表示负相关，1 为正相关，0 为不相关
edges <- edges %>%
  mutate(corr = sample(-1:1, size = n(), replace = TRUE))
graph <- tbl_graph(nodes = nodes, edges = edges)

# 1. linear布局
# 对于线性布局，即将所有节点放置在一条直线上，然后使用 geom_edge_arc 将边绘制成弧形
ggraph(graph, layout = 'linear') + 
  geom_edge_arc(aes(colour = factor(corr))) +
  geom_node_point() +
  theme_graph()

# 设置 circular = TRUE，转换为径向表示
ggraph(graph, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(colour = factor(corr))) +
  geom_node_point() +
  coord_fixed() +
  theme_graph()

# 2. partition 布局
# 分区布局是一种显示分层结构的方式，每一层都会对前一层的切片进行分割。
# 例如，我们将节点的形状设置为条块 geom_node_tile
graph <- tbl_graph(flare$vertices, flare$edges)
ggraph(graph, 'partition') + 
  geom_node_tile(aes(fill = depth), size = 0.25) +
  theme_graph()

# 分区布局的圆形表示，注意节点的形状变成了圆弧条形 geom_node_arc_bar
# 并不是所有的布局都支持圆形表示，下面的布局将会忽略 circular 参数
ggraph(graph, 'partition', circular = TRUE) + 
  geom_node_arc_bar(aes(fill = depth), size = 0.25) +
  coord_fixed() +
  theme_graph()

# 3. node — edge 布局
# 我们可以直接使用 igraph 中定义的布局算法，例如
ggraph(graph, layout = "graphopt") + 
  geom_edge_link(aes(colour = factor(corr)), show.legend = FALSE) +
  geom_node_point() + 
  theme_graph()

# lay <- c('stress', 'fr', 'lgl', 'graphopt')
# plot_fun <- function(g, layout) {
#   p <- ggraph(g, layout = layout) + 
#     geom_edge_link(aes(colour = factor(corr)), show.legend = FALSE) +
#     geom_node_point() + 
#     labs(caption = paste0('Layout: ', layout)) +
#     theme_graph()
#   return(p)
# }
# 
# library(grid)
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(2, 2)))
# 
# for (i in seq_along(lay)) {
#   x = (i - 1) %/% 2 + 1
#   y = (i - 1) %% 2 + 1
#   p = plot_fun(graph, lay[i])
#   print(p, vp = viewport(layout.pos.row = x, layout.pos.col = y))
# }

# 4. 蜂巢图
# 蜂巢图也是一种 node-edge 图，它使用的是节点的信息，将节点进行分类
graph <- graph %>%
  mutate(friends = ifelse(
    centrality_degree(mode = 'all') < 3, "few",
    ifelse(centrality_degree(mode = 'all') > 3, "many", "medium")
  ))

ggraph(graph, 'hive', axis = friends) + 
  geom_edge_hive(aes(colour = factor(corr))) + 
  geom_axis_hive(aes(colour = friends), size = 2, label = FALSE) + 
  coord_fixed() +
  theme_graph()

# 5. 焦点布局
# 将焦点聚集在一个或一组节点上，其他节点则相对于该位置放置
ggraph(graph, 'focus', focus = node_is_center()) + 
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r), data.frame(r = 1:5), colour = 'grey') +
  geom_edge_link(aes(colour = factor(corr)), show.legend = FALSE) + 
  geom_node_point() + 
  coord_fixed() +
  theme_graph()

# 6. 层次布局
# 圆堆积图：以包含的方式来展示层次结构
graph <- tbl_graph(flare$vertices, flare$edges)
ggraph(graph, 'circlepack', weight = size) + 
  geom_node_circle(aes(fill = factor(depth)), size = 0.25, n = 50) + 
  coord_fixed() +
  theme_graph()

# 堆积树状图
ggraph(graph, 'circlepack', weight = size) + 
  geom_edge_link() + 
  geom_node_point(aes(colour = factor(depth))) +
  coord_fixed() +
  theme_graph()

# 树状图 
# 矩形层次关系
ggraph(graph, 'treemap', weight = size) + 
  geom_node_tile(aes(fill = factor(depth)), size = 0.25) +
  theme_graph()

ggraph(graph, 'treemap', weight = size) + 
  geom_edge_link() + 
  geom_node_point(aes(colour = factor(depth))) +
  theme_graph()

# 根据树的不同深度进行绘制
ggraph(graph, 'tree') + 
  geom_edge_diagonal() +
  theme_graph()

# 层次聚类树状图
dendrogram <- hclust(dist(iris[,1:4]))
ggraph(dendrogram, 'dendrogram', height = height) + 
  geom_edge_elbow() +
  theme_graph()

# 圆形树状图
ggraph(dendrogram, 'dendrogram', circular = TRUE) + 
  geom_edge_elbow() + 
  coord_fixed() +
  theme_graph()

# 系统发育树，不存在根节点，使用无根布局
tree <- create_tree(100, 2, directed = FALSE) %>% 
  activate(edges) %>% 
  mutate(length = runif(n()))

ggraph(tree, 'unrooted', length = length) + 
  geom_edge_link() +
  theme_graph()

# 7. 矩阵布局
# 矩阵布局是将节点放置在对角线，如果对应位置的两个节点之间有交叠，那矩阵中对应的行列将会绘制一个点或矩形。
# 不同的节点顺序，会影响矩阵布局的形状
graph <- create_notable('zachary')
ggraph(graph, 'matrix', sort.by = node_rank_leafsort()) + 
  geom_edge_point(mirror = TRUE) + 
  coord_fixed() +
  theme_graph()

ggraph(graph, 'matrix', sort.by = node_rank_spectral()) + 
  geom_edge_point(mirror = TRUE) + 
  coord_fixed() +
  theme_graph()

# 7. Fabric 布局
# Fabric 布局是一种可扩展的特殊的 BioFabric 布局。其特殊的地方在于，它将节点表示为水平线，边表示为连接两个水平线的竖直线
ggraph(graph, 'fabric', sort.by = node_rank_fabric()) + 
  geom_node_range(colour = 'grey') + 
  geom_edge_span(end_shape = 'square') + 
  coord_fixed() +
  theme_graph()

# 添加一份重复的阴影边
ggraph(graph, 'fabric', sort.by = node_rank_fabric(), shadow.edges =TRUE) + 
  geom_node_range(colour = 'grey') + 
  geom_edge_span(aes(filter = shadow_edge), colour ='lightblue' , end_shape = 'square') + 
  geom_edge_span(aes(filter = !shadow_edge), end_shape = 'square') + 
  coord_fixed() +
  theme_graph()

####应用####
# https://zhuanlan.zhihu.com/p/375755699