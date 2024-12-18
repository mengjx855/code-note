# Encoding: utf-8
# Creator：Jinxin Meng
# Created Date: 2022-10-21
# Modified Date: 2022-10-21

library(igraph)

#### 使用数据框创建igraph对象 ####
# 可以使用graph_from_data_frame() (或from_data_frame()) 函数从一个或多个包含节点或边的数据框中创建igraph对象。
# graph_from_data_frame(d, directed = TRUE, vertices = NULL)
# from_data_frame(...)
# 根据 vertices 参数是否为 null，会有两种不同的操作模式;
# 如果 vertices = NULL，则数据框的前两列表示边列表，剩下的列表示边的属性;
# 如果 vertices 不为空，则必须是一个包含节点信息的数据框。第一列将作为节点的名称，显示在图中，其他列作为节点的属性，且只保留边列表 d 中包含这些节点的边

# 定义节点
actors <- data.frame(name = c("Alice", "Bob", "Cecil", "David", "Esmeralda"),
                     age = c(48, 33, 45, 34, 21),
                     gender = c("F", "M", "F", "M", "F"))

# 定义边
relations <- data.frame(from = c("Bob", "Cecil", "Cecil", "David", "David", "Esmeralda"),
                        to = c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
                        friendship = c(4, 5, 5, 2, 1, 1),
                        advice = c(4, 5, 5, 4, 2, 3))

# 将数据转换为igraph对象
igraph <- graph_from_data_frame(relations, directed = TRUE, vertices = actors)
plot(igraph)

igraph <- graph_from_data_frame(relations, directed = TRUE, vertices = NULL)
plot(igraph)

# as_data_frame用于将igraph对象转换为一个或多个数据框。
# 根据 what 参数来进行转换，可选的值为 c("edges", "vertices", "both")
as_data_frame(igraph, what = "edges")
as_data_frame(igraph, what = "vertices")



####使用边列表创建igraph对象####
# graph_from_edgelist用于从边列表来创建igraph对象，接受一个包含两列的矩阵，每行定义了一条边。
mat <- matrix( c("foo", "bar", "bar", "foobar"), nc = 2, byrow = TRUE)
igraph <- graph_from_edgelist(mat)
plot(igraph)

####邻接矩阵创建igraph对象####
# graph_from_adjacency_matrix是一个灵活的函数，能够使用邻接矩阵来创建igraph对象
# graph_from_adjacency_matrix(adjmatrix, 
#                             mode = c("directed", "undirected", "max", "min", "upper", "lower", "plus"), 
#                             weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
# 节点的顺序将依照矩阵的顺序排序，mode 参数各个值的含义为：
# directed：有向图
# undirected：无向图，与 max 一样
# max：将 max(A(i,j), A(j,i)) 作为无向图边的值
# upper：使用上三角来构建无向图
# lower：使用下三角来构建无向图
# min：将 min(A(i,j), A(j,i)) 作为无向图边的值
# plus：将 A(i,j)+A(j,i) 作为无向图边的值
# 如果 weighted 参数为 TRUE，会创建加权图
adj <- matrix(sample(0:1, 100, replace = TRUE, prob = c(0.9, 0.1)), nc = 10)
igraph <- graph_from_adjacency_matrix(adj, weighted = T)
plot(igraph)

# 获取igraph对象节点的权重
# E()取边
E(igraph)$weight
# as_data_frame(igraph, what = "edges")



####图层的方式创建网络####
# 我们还可以使用类似ggplot2图层的方式，使用+来添加节点和边
# 我们首先使用 make_empty_graph 创建一个空的图，默认为有向图，设置 directed = FALSE 来绘制无向图
# 然后，使用 + 来添加 vertices() 创建的结点和 edges() 创建的边
igraph <- make_empty_graph(directed = FALSE) +
  vertices(LETTERS[1:10]) +
  edges('A','B', 'B','D', 'C','D', 'D','E', 'E','G', 'F','G', 'G','H', 'H','I', 'I','J')
plot(igraph)

# 还可以使用 - 来删除某些节点或者边
# 删除边的使用，传递给 edges() 的值略有不同，要删除的边为使用 | 连接两个节点名称的字符串形式
igraph <- make_empty_graph(directed = FALSE) + 
  vertices(LETTERS[1:10]) +
  edges('A','B', 'B','D', 'C','D', 'D','E', 'E','G', 'F','G', 'G','H', 'H','I', 'I','J') -
  vertices(LETTERS[5]) -
  edges('A|B', 'I|J')
plot(igraph)



####查看igraph对象####
# 在第一行中，16进制数后面的字符分别代表 (DN-- 10 9 --)
# 第一个-：D或U表示有向和无向
# 第二个-：N表示命名图，即节点有name属性
# 第三个-：W表示加权图，即边存在weight属性
# 第四个-：B表示bipartite (two-mode) graph，即节点存在type属性
# 后面的两个数值 (10 9) 表示节点和边的数量
# 最后是网络的名称
# 第二行：表示节点和边的属性，其中括号中的字符分别表示：
# g: graph; v: vertex; e: edge; n: numeric; c: character; l: logical; x: complex
# 例子中的name为节点的字符型属性
igraph <- graph(edges = c('A', 'B', 'B','D', 'C','D', 'D','E', 'E','G', 'F','G', 'G','H', 'H','I', 'I','J'))
plot(igraph)
igraph



####获取图的大小####
# 获取图的顺序，即节点数
gorder(igraph)
# 获取图的大小，即边的数量
gsize(igraph)



####提取子图####
# 以传递节点 ID 的方式，获取指定节点及这些节点之间的边，构建子图
induced_subgraph(igraph, 1:5) %>% plot(.)
# 以指定边 ID 的方式，获取包含指定边及边所涉及到的节点，构建子图
subgraph.edges(igraph, 1:5) %>% plot(.)



####邻接矩阵查询节点和边的属性####
# 最直接的表示图的方法是邻接矩阵，行列与节点对应，如果两个节点之间存在连接，则矩阵中对应位置的值为 1，否则为 0。对于加权图来说，值代表的是边的权重。
# 对于无向图来说，邻接矩阵总是对称的，而对于有向图来说并不一定是对称的。
# 使用[]可以返回对象的邻接矩阵，也可以使用get.adjacency(g)来获取
igraph[]

# 获取一条对应的边
# 0表示不存在A到F的边，我们可以将其设置为 1，意味着会添加这条边
# 如果要删除一条边，可以将邻接矩阵对应位置的值设置为0
igraph['A', 'F']
igraph['A', 'F'] <- 1
igraph['D', 'E']
igraph['D', 'E'] <- 0

# 查询多条边会返回一个矩阵
igraph[c('A', 'B', 'C'), c('A', 'B', 'C', 'D')]

# 获取一个节点与其他所有节点之间的连接
# 可以用name索引，也能以数字索引的方式
igraph['A', ]
igraph[ ,'G']
igraph[2, ]
igraph[2, -1]

# 条件查询边
degree(igraph)
igraph[degree(igraph) > 2, degree(igraph) < 3]

# 配对节点
# 将会查询的边为 A-C、B-D、C-E，需要保证 from 和 to 的长度一致
igraph[from = c('A', 'B', 'C'), to = c('C', 'D', 'E')]

# 获取邻接节点
neighbors(igraph, 'G')



####邻接表查询节点和边的属性####
# 邻接表是以表的形式来表示图的，只存储每个节点中与其相连的节点，用[[]]查看邻接表
igraph[[]]

# 获取节点的邻接节点
igraph[[c('A', 'B')]]

# 不同方向的邻接节点
igraph[['G', ]]
igraph[[, 'G']]

# 显示边
igraph[['A', edges = TRUE]]

# 获取两个节点集合之间的边
igraph[[c('A', 'B', 'C'),c('C', 'D', 'E', 'F'), edges = TRUE]]



####节点和边的属性####
# 使用 E() 和 V() 函数，可以获取 igraph 对象的边和节点
# V() 的返回值是根据结点的 ID 进行了排序。由于我们的节点是字符型的，会根据节点添加的顺序自动为节点设置 ID，第一个节点的 ID 为 1
E(igraph)
V(igraph)

# 可以使用两个节点名称之间添加一个竖线组成的字符串来引用连接两个节点的边，如果是无向图，则两个节点的顺序不影响边的引用
# 使用 ends() 可以获取边的矩阵表示
igraph %>% ends('A|B')
igraph %>% ends(E(igraph))

# 对于有向图，head_of 和 tail_of 可以获取边的头尾两个端点，有箭头的一端为头
igraph %>% tail_of('A|B')
igraph %>% head_of('A|B')

# 使用 neighbors 来获取节点的邻接节点
igraph %>% neighbors('D', mode = "in")
igraph %>% neighbors('D', mode = "out")

# 对于无向图，两种模式的邻接节点是一样的
# 计算边和节点的数量
ecount(igraph)
vcount(igraph)



####节点和边的序列操作####
# 节点序列操作
V(igraph)[1:4]
V(igraph)[1:3, 5:7]
V(igraph)[c('A', 'C', 'D', 'G')]

# 根据邻接节点索引
V(igraph)[.nei('D')]
V(igraph)[.innei('D')]
V(igraph)[.outnei('D')]

# 查询条件
V(igraph)[degree(igraph) > 2]

# 反转、去重和集合操作
rev(V(igraph))
unique(V(igraph)['A', 'A'])
union(V(igraph)[1:3], V(igraph)[6:8])
intersection(V(igraph)[1:5], V(igraph)[4:8])
difference(V(igraph), V(igraph)[1:5])

# 获取边涉及的结点
V(igraph)[.inc('A|B'), .inc('D|E')]

# 边序列操作与节点类似
E(igraph)[1:4]
E(igraph)[c('A|B', 'G|H', 'H|I')]

# 获取涉及到节点的边
E(igraph)[.inc('D')]
E(igraph)[.from('D')]
E(igraph)[.to('D')]

# 条件索引
E(igraph)[seq_len(gsize(igraph)) %% 2]

# 获取两个节点集合之间存在的边
# 不管方向
E(igraph)[V(igraph)['D'] %--% V(igraph)['B', 'C', 'E']]
# 左到右
E(igraph)[V(igraph)['D'] %->% V(igraph)['B', 'C', 'E']]
# 右到左
E(igraph)[V(igraph)['D'] %<-% V(igraph)['B', 'C', 'E']]

# 集合操作、去重和反转和节点操作一样
# 注意：我们在中括号中使用的简写函数只能在对应条件下使用，无法单独使用



###网络元数据属性####
# 我们可以将节点和边的元数据（附加信息）作为属性值的方式添加到节点和边的属性中
# 使用$可以获取和设置属性，例如，获取节点的 name 属性
V(igraph)$name

# 为节点新建一个名为 class 的属性
V(igraph)$class <- rep(c("I", "II"), each = 5)
neighbors(igraph, 'G', mode = "all")$class

# 为边设置属性值
# 也可以使用 set_edge_attr 和 set_vertex_attr 函数来设置，效果是一样的
# 可以使用 delete_edge_attr 和 delete_vertex_attr 删除属性
E(igraph)$type <- rep(c('activation', 'repression', 'inhibition'), each = 3)
E(igraph)$weight <- sample(3:6, 9, replace = TRUE)

# 获取所有节点和边的属性
edge_attr(igraph)
vertex_attr(igraph)

# 获取图属性
igraph$name
graph_attr(igraph)
graph_attr_names(igraph)
graph_attr(igraph, 'name')

# 设置图的属性
igraph$name <- "demo"
igraph <- set_graph_attr(igraph, 'name', 'demo')

# 删除图属性
igraph <- set_graph_attr(igraph, 'some', 'something')
graph_attr(igraph)
igraph <- delete_graph_attr(igraph, 'some')
graph_attr(igraph)

# 上面的属性设置都是在已经构建完图之后，我们也可以在构建图的时候添加属性。例如
# 在这里，我们使用add_edges和add_vertices函数来添加边和节点，对应的，可以使用delete_edges和delete_vertices来删除边和节点
igraph <- make_empty_graph(n = 5) %>%
  add_edges(c(1,2, 2,3, 3,4, 4,5)) %>%
  set_edge_attr("color", value = "red") %>%
  add_edges(c(5,1), color = "green") %>%
  set_vertex_attr("color", value = "orange") %>%
  add_vertices(3, color = "red") %>%
  add_vertices(2, color = "green")
plot(igraph)

# 有几个特殊的属性：
# 属性    类型    描述
# layout  graph   图的布局，可以是矩阵或者函数
# color   vertex  节点的颜色
# name    vertex  节点的名称
# shape   vertex  节点的形状
# type    vertex  节点的分组
# color   edge    边的而颜色
# weight  edge    边的权重



####igraph布局和绘图
# 图是一种抽象的数学结构，不同对象之间通过线条连接起来，而对象在图中并没有固定的位置表示，不同的放置位置显示出的效果通常是不一样的。
# 选择一种优秀的布局方式，可以让图形呈现出更好的效果，而 igraph 的工作方式是通过一类 node-edge 的算法来进行布局的。
# 算法会将节点作为二维或三维空间上的点，使用直线或曲线来连接两个相邻的节点。对于有向图来说，带箭头的线表示连接方向。在边两端的节点可以由不同的几何图形来表示，而一些重要的节点或边的属性可以用来设置图形参数值。
# 图的可视化通常由三个步骤组成：
# 1. 找到节点在二维或三维空间上合适的排列方式，这一步可能是最重要的。好的布局往往能够解释一些有趣的现象，如对称性、密集连接区域。
# 2.将节点和边的重要属性映射到图像上
# 3.安排好节点与边的绘制顺序

# 布局
# 布局算法用于寻找合适的排列方式，但是找到优秀的排列确实不容易，因此，大部分的算法都是通过间接测量来评估布局的好坏。而且是一种启发式的方法来求解，所以并不是每次都能找到最优解。
# 大部分的布局算法只适用于较小的图，较大的图需要一些特殊技术来处理。
# igraph 中所有布局算法都是layout_*()形式，每个布局算法都会返回一个layout实例，类似列表型的对象，包含了每个节点在图中的 x、y 坐标。
# 布局函数layout_ 功能
# as_brpartitle   节点有type属性才能用
# as_star         星形布局
# as_tree         树状    
# in_circle       圆形
# nicely          自动选择合适的算法
# on_grid         网格布局
# on_sphere       球形
# randomly        随机
# with_dh         Davidson-Harel算法
# with_fr         Fruchterman-Reingold算法
# with_gem        GEM算法
# with_graphopt   graphopt算法
# with_kk         Kamada-Kawai算法
# with_lgl        大图布局
# with_mds        多维缩放布局
# with_sugiyama   分层定向无循环图的Sugiyama算法
# 布局是通过 plot() 函数的 layout 参数来控制的
igraph <- sample_gnm(n = 15, m = 25)
plot(igraph)
plot(igraph, layout = layout_randomly)
plot(igraph, layout = layout_in_circle)
plot(igraph, layout = layout_with_fr) # 力导向布局，最常用的是 Fruchterman-Reingold 算法
plot(igraph, layout = layout_with_kk) # 另一种比较常用的力导向算法 Kamada Kawai

par(mfrow=c(1, 3))
tree <- make_tree(20, 3)
plot(tree, layout=layout_as_tree)
plot(tree, layout=layout_as_tree(tree, flip.y=FALSE)) # 自底向上
plot(tree, layout=layout_as_tree(tree, circular=TRUE)) # 圆形排列

# 设置根节点
par(mfrow = c(1, 2))
tree2 <- make_tree(10, 3) + make_tree(10, 2)
plot(tree2, layout = layout_as_tree)
plot(tree2, layout = layout_as_tree(tree2, root = c(1, 11), rootlevel = c(2, 1)))

# 传递位置矩阵
igraph <- sample_gnm(n = 15, m = 25)
l <- cbind(1:vcount(igraph), c(1, vcount(igraph):2))
plot(igraph, layout = l)



####igraph函数的应用####
# https://zhuanlan.zhihu.com/p/374836002

