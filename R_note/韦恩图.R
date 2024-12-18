# Encoding: utf-8
# Creator：Jinxin Meng
# Created Date: 2021-10-19
# Modified Date: 2023-12-15

# venn ----
set.seed(2021)
genes <- paste("gene", 1:1000, sep = "")
x <- list(A = sample(genes, 300), B = sample(genes, 525), C = sample(genes, 440), D = sample(genes, 350))

# ggvenn package
library(ggvenn)
ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4)

# ggVennDiagram
library(ggVennDiagram)
ggVennDiagram(x, label_alpha = 0)

genes <- paste("gene",1:1000,sep="")
set.seed(2022)
x <- list(A=sample(genes,300),
          B=sample(genes,525),
          C=sample(genes,440))

library(ggplot2)
ggVennDiagram(x,category.names = c("A","B","C"),
              label = "count", 
              label_color = "black",
              label_alpha = 0,
              edge_lty = "dashed", 
              edge_size = 1) +
  scale_fill_gradient(low="white",high = "#b9292b",name = "gene count")

# VennDiagram
library(VennDiagram)
venn.diagram(x, filename = "venn-4-dimensions.png")
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(x)
display_venn(x, category.names = c("Set 1" , "Set 2 " , "Set 3", "Set 4"), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))


# Note: most examples are listed as dontrun to meet CRAN requirements, 
# but all should work as-is!

# compact and minimal notation
## Not run: 
venn.plot <- venn.diagram(
  list(A = 1:150, B = 121:170), 
  filename = tempfile(
    pattern = 'Venn_2set_simple',
    fileext = '.tiff'
  )
)
venn.plot <- venn.diagram(
  list(A = 1:150, B = 121:170, C = 101:200), 
  filename = tempfile(
    pattern = 'Venn_3set_simple',
    fileext = '.tiff'
  )
)
## End(Not run)
# a more elaborate two-set Venn diagram with title and subtitle
venn.plot <- venn.diagram(
  x = list(
    "A" = 1:100,
    "B" = 96:140
  ),
  filename = tempfile(
    pattern = 'Venn_2set_complex',
    fileext = '.tiff'
  ),
  scaled = TRUE,
  ext.text = TRUE,
  ext.line.lwd = 2,
  ext.dist = -0.15,
  ext.length = 0.9,
  ext.pos = -4,
  inverted = TRUE,
  cex = 2.5,
  cat.cex = 2.5,
  rotation.degree = 45,
  main = "Complex Venn Diagram",
  sub = "Featuring: rotation and external lines",
  main.cex = 2,
  sub.cex = 1
);
## Not run: 
# sample three-set Euler diagram
venn.plot <- venn.diagram(
  x = list(
    "Num A" = paste("Num", 1:100),
    "Num B" = c(paste("Num", 61:70), paste("Num", 71:100)),
    "Num C" = c(paste("Num", 41:60), paste("Num", 61:70))),
  euler.d = TRUE,
  filename = tempfile(
    pattern = 'Euler_3set_simple',
    fileext = '.tiff'
  ),
  cat.pos = c(-20, 0, 20),
  cat.dist = c(0.05, 0.05, 0.02),
  cex = 2.5,
  cat.cex = 2.5,
  reverse = TRUE
);
# sample three-set Euler diagram
venn.plot <- venn.diagram(
  x = list(
    A = c(1:10),
    B = c(11:90),
    C = c(81:90)
  ),
  euler.d = TRUE,
  filename = tempfile(
    pattern = 'Euler_3set_scaled',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.pos = 0
);
## End(Not run)
# sample four-set Venn Diagram
A <- sample(1:1000, 400, replace = FALSE);
B <- sample(1:1000, 600, replace = FALSE);
C <- sample(1:1000, 350, replace = FALSE);
D <- sample(1:1000, 550, replace = FALSE);
E <- sample(1:1000, 375, replace = FALSE);
venn.plot <- venn.diagram(
  x = list(
    A = A,
    D = D,
    B = B,
    C = C
  ),
  filename = tempfile(
    pattern = 'Venn_4set_pretty', 
    fileext = '.tiff'
  ),
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", 
                "white", "white", "white", "white", "darkblue", "white", 
                "white", "white", "white", "darkgreen", "white"),
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.07,
  cat.fontfamily = "serif",
  rotation.degree = 270,
  margin = 0.2
);
# sample five-set Venn Diagram
venn.plot <- venn.diagram(
  x = list(
    A = A,
    B = B,
    C = C,
    D = D,
    E = E
  ),
  filename = tempfile(
    pattern = 'Venn_5set_pretty',
    fileext = '.tiff'
  ),
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  alpha = 0.50,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.05
);
# Complex three-way Venn with labels & sub-/super-scripts	
venn.plot <- venn.diagram(
  x = list(
    I = c(1:60, 61:105, 106:140, 141:160, 166:175, 176:180, 181:205, 
          206:220),
    II = c(531:605, 476:530, 336:375, 376:405, 181:205, 206:220, 166:175, 
           176:180),
    III = c(61:105, 106:140, 181:205, 206:220, 221:285, 286:335, 336:375, 
            376:405)
  ),
  category.names = c(
    expression( bold('A'['1: subscript']) ),
    expression( bold('B'^'2: going up') ),
    expression( paste(bold('C'^'3'), bold('X'['i' <= 'r'^'2']^'2') ) )
  ),
  filename = tempfile(
    pattern = 'Fig3-1_triple_labels_sub_and_superscripts',
    fileext = '.tiff'
  ),
  output = TRUE,
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = 'lzw',
  units = 'px',
  lwd = 6,
  lty = 'blank',
  fill = c('yellow', 'purple', 'green'),
  cex = 3.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
);
# Complex 3-way Venn using expressions
venn.plot <- venn.diagram(
  x = list(
    "Num A" = paste("Num", 1:100),
    "Num B" = c(paste("Num", 61:70), paste("Num", 71:100)),
    "Num C" = c(paste("Num", 41:60), paste("Num", 61:70))),
  category.names = c(
    expression( bold('A'['1']) ),
    expression( bold('A'['2']) ),
    expression( bold('A'['3']) )
  ),
  euler.d = TRUE,
  filename = tempfile(
    pattern = 'Fig3-2_Euler_3set_simple_with_subscripts',
    fileext = '.tiff'
  ),
  cat.pos = c(-20, 0, 20),
  cat.dist = c(0.05, 0.05, 0.02),
  cex = 2.5,
  cat.cex = 2.5,
  reverse = TRUE
);
## Not run: 
# Example to print to screen
venn.plot <- venn.diagram(
  x = list(
    sample1 = c(1:40),
    sample2 = c(30:60)
  ),
  filename = NULL,
  disable.logging = TRUE
);
# Save picture to non-TIFF file type
# currently working on adding this functionality directly into venn.diagram
venn.plot <- venn.diagram(
  x = list (
    A = 1:10,
    B = 6:25
  ),
  filename = NULL,
  disable.logging = TRUE
);
jpeg(tempfile(pattern = 'venn_jpeg', fileext = '.jpg'));
grid.draw(venn.plot);
dev.off();
## End(Not run)
#dontrun-starts-here
### NB: All figures from the paper can be run, but are turned off from
###     automatic execution to reduce burden on CRAN computing resources.
## Not run: 
# Figure 1A
venn.plot <- venn.diagram(
  x = list(
    Label = 1:100
  ),
  filename = tempfile(
    pattern = '1A-single_Venn',
    fileext = '.tiff'
  ),
  col = "black",
  lwd = 9,
  fontface = "bold",
  fill = "grey",
  alpha = 0.75,
  cex = 4,
  cat.cex = 3,
  cat.fontface = "bold",
);
# Figure 1B
venn.plot <- venn.diagram(
  x = list(
    X = 1:150,
    Y = 121:180
  ),
  filename = tempfile(
    pattern = '1B-double_Venn',
    fileext = '.tiff'
  ),
  lwd = 4,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.75,
  label.col = "white",
  cex = 4,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("cornflowerblue", "darkorchid1"),
  cat.cex = 3,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 14)
);
# Figure 1C
venn.plot <- venn.diagram(
  x = list(
    R = c(1:70, 71:110, 111:120, 121:140),
    B = c(141:200, 71:110, 111:120, 201:230),
    G = c(231:280, 111:120, 121:140, 201:230)
  ),
  filename = tempfile(
    pattern = '1C-triple_Venn',
    fileext = '.tiff'
  ),
  col = "transparent",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkblue", "white",
                "white", "white", "darkgreen"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkblue", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06, 0.03),
  cat.pos = 0
);
# Figure 1D
venn.plot <- venn.diagram(
  x = list(
    I = c(1:60, 61:105, 106:140, 141:160, 166:175, 176:180, 181:205, 
          206:220),
    IV = c(531:605, 476:530, 336:375, 376:405, 181:205, 206:220, 166:175, 
           176:180),
    II = c(61:105, 106:140, 181:205, 206:220, 221:285, 286:335, 336:375, 
           376:405),
    III = c(406:475, 286:335, 106:140, 141:160, 166:175, 181:205, 336:375, 
            476:530)
  ),
  filename = tempfile(
    pattern = '1D-quadruple_Venn',
    fileext = '.tiff'
  ),
  col = "black",
  lty = "dotted",
  lwd = 4,
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
                "white", "white", "darkblue", "white",
                "white", "white", "white", "darkgreen", "white"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 2.5,
  cat.fontfamily = "serif"
);
# Figure 2-1
venn.plot <- venn.diagram(
  x = list(
    A = 1:105,
    B = 101:115
  ),
  filename = tempfile(
    pattern = '2-1_special_case_ext-text',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.pos = c(-20, 20),
  ext.line.lty = "dotted",
  ext.line.lwd = 2,
  ext.pos = 12,
  ext.dist = -0.12,
  ext.length = 0.85
);
# Figure 2-2
venn.plot <- venn.diagram(
  x = list(
    A = 1:100,
    B = 1:10
  ),
  filename = tempfile(
    pattern = '2-2_special_case_pairwise-inclusion',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.pos = 0
);
# Figure 2-3
venn.plot <- venn.diagram(
  x = list(
    A = 1:150,
    B = 151:250
  ),
  filename = tempfile(
    pattern = '2-3_special_case_pairwise-exclusion',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.pos = c(0, 0),
  cat.dist = 0.05
);
# Figure 2-4
venn.plot <- venn.diagram(
  x = list(
    A = c(1:50, 101:140, 141:160, 161:170),
    B = c(171:230, 101:140, 161:170, 291:320),
    C = c(141:160, 161:170, 291:320)
  ),
  filename = tempfile(
    pattern = '2-4_triple_special_case-001',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.dist = c(0.05, 0.05, -0.1)
);
# Figure 2-5
venn.plot <- venn.diagram(
  x = list(
    A = c(1:100),
    B = c(61:70, 71:100),
    C = c(41:60, 61:70)
  ),
  filename = tempfile(
    pattern = '2-5_triple_special_case-012AA',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.pos = c(-25, 0, 30),
  cat.dist = c(0.05, 0.05, 0.02)
);
# Figure 2-6
venn.plot <- venn.diagram(
  x = list(
    A = c(1:90),
    B = c(1:25),
    C = c(1:5)
  ),
  filename = tempfile(
    pattern = '2-6_triple_special_case-022AAAO',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.pos = 0,
  cat.dist = c(0.03, 0.03, 0.01)
);
# Figure 2-7
venn.plot <- venn.diagram(
  x = list(
    A = c(1:20),
    B = c(21:80),
    C = c(81:210)
  ),
  filename = tempfile(
    pattern = '2-7_triple_special_case-100',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.dist = 0.05
);
# Figure 2-8
venn.plot <- venn.diagram(
  x = list(
    A = c(1:80),
    B = c(41:150),
    C = c(71:100)
  ),
  filename = tempfile(
    pattern = '2-8_triple_special_case-011A',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.dist = c(0.07, 0.07, 0.02),
  cat.pos = c(-20, 20, 20)
);
# Figure 2-9
venn.plot <- venn.diagram(
  x = list(
    A = c(1:10),
    B = c(11:90),
    C = c(81:90)
  ),
  filename = tempfile(
    pattern = '2-9_triple_special_case-121AO',
    fileext = '.tiff'
  ),
  cex = 2.5,
  cat.cex = 2.5,
  cat.pos = 0,
  cat.dist = c(0.04, 0.04, 0.02),
  reverse = TRUE
);
#dontrun-ends-here
## End(Not run)

# upset ----
upset(plot_dat, # 输入数据，可以为list(使用函数fromList()), 可以用邻接矩阵那种格式的
      nsets = 6, # 集合数
      # intersections = group_order, # 指定某些交集

      # 右上，条形图参数
      number.angles = 0, # 上方柱子文字的倾斜程度
      mainbar.y.label = "ARGs number",  # 上方柱状图的注释
      queries = list(# 条形图的颜色
        list(query = elements, params = c("direction","Control"), color = "#74add1", active = T)),
      
      # 左下，条形图的参数
      sets.x.label = "ARGs number",  # 条形图x轴标签
      sets.bar.color = group_color,  # 左下条形图颜色 有名字的颜色向量
      sets = group_order, # 指定特殊的集合, 有排序的作用，但实现排序，需要使用keep.order参数 
      keep.order = T, # 排序，按照sets参数设置的顺序
      
      # 右下，交集散点图的参数
      nintersects = 50, # 右下 交集的数量
      point.size = 1.5,  # 右下点大小
      line.size = 0.4, # 右下线的粗细
      order.by = c("freq", "degree"), # 对交集进行排序，默认为升序，freq指定的是按照交集数量的多少排序，degree指的是按照交集集合数量的进行排序。这两个参数是先后顺序
      decreasing = c(T, F), # 调整交集排序的方向，和order.by一一对应，决定是否排序。
      set.metadata = list( # 交集的背景颜色
      data = upset_metadata, 
      plots = list(list(type = "matrix_rows", column = "set", alpha = .4,
                        colors = c(HuangR_2020_AS = "#e41a1c", ZhouC_2020_AS = "#377eb8", ZhuQ_2021_GD = "#984ea3",
                                   ChuY_2021_Gout = "#FFD320", LiuP_2021_MG = "#ff7f00", ChenB_2020_SLE = "#f781bf")) ) ),
      # 总体参数
      text.scale = c(1.2, 1, 1.1, 1, 1.3, 0.75), # text.scale 参数值的顺序为:柱状图的轴标签和刻度,条形图的轴标签和刻度,集合名称,柱子上方表示交集大小的数值
      mb.ratio = c(.5, .5),  # 控制上下部分图形所占的比例
      )

# 参数 usage
# data 数据
# nsets 从数据中选择 n 个最大的集合
# nintersects 设置需要绘制的交集数目，如果设置为 NA，则绘制所有的交集
# sets 设置集合的名称，如:c("Name1"，"Name2")
# keep.order 是否保持 sets 参数中集合的顺序，默认为 FALSE，根据集合大小排序
# set.metadata 绘制 metadata
# intersections 选择需要显示的交集。如:list(list("Set name1"，"Set name2"), list("Set name1"，"Set name3"))
# matrix.color 矩阵中交集点的颜色
# main.bar.color 柱状图柱子的颜色
# mainbar.y.label 柱状图的 y 轴标签
# mainbar.y.max 柱状图 y 轴最大值
# sets.bar.color 条形图条形的样色
# sets.x.label 条形图的 x 轴标签
# point.size 矩阵中圆圈的大小
# line.size 矩阵中连接圆圈的线条宽度
# mb.ratio 柱状图与矩阵点图之间的比例大小 (如 c(0.55，0.45))
# expression 用于选取交集或元素的子集的表达式("ColName >3")
# att.pos 属性图的位置。如果为 NULL 或"bottom"在下方绘制.如果是"top"将会绘制在上方
# att.color 为查询到的数据属性直方图或散点图颜色，默认为上方柱状图的颜色
# order.by 矩阵点图中交集的排序方式，排序方式有 freq、degree 或两者
# decreasing order.by变量的排序方式."freq"为降序, "degree"为升序
# show.numbers 在柱状图上显示交集的大小
# number.angles 设置柱状图上的柱形上方数字标签的角度
# group.by 数据的分组方式 ("degree"或"sets")
# cutoff 设置每个 group.by 分组中交集的数目
# queries 设置查询条件
# query.legend 查询图例的位置
# shade.color 矩阵点图阴影颜色
# shade.alpha 矩阵点图阴影透明度
# matrix.dot.alpha 空交集点的透明度
# empty.intersections 如果要绘制空交集，将参数值设置为"on"
# color.pal 设置属性图的颜色画板
# boxplot.summary 以交集组合作为 x 轴，绘制指定属性的箱线图，一次只能绘制两个
# attribute.plots 添加自定义 ggplot2 图片
# scale.intersections 交集大小的转换函数，"identity"，"log10"，"log2"
# scale.sets 集合大小的转换函数，"identity"，"log10"，"log2"
# text.scale 设置文本大小，包舍6个值 c(intersection size title,intersection size tick labels, set size title, set size ticklabels，set names,numbers above bars)
# set size.angles 设置条形图 x 轴文本的角度
# set size.show 是否在条形图的条形上显示集合大小
# set size.numbers size set_size.show = TRUE时，调整显示的数字的大小
# set size.scale max 设置条形图 x 轴最大值

# 排序问题 order.by()
# order.by参数可以传入两个参数，一个是freq，代表着每个交集内元素数量的大小，也就是上方柱子的高低，另一个参数是degree，这个参数可以指定交集的数量，也就是下方几个点连在一起
# order.by指定的参数和decreasing对应，默认是升序

# 如果想给下方集合点的背后加颜色
# 需要构建一个数据集 (upset_metadata), id为set的名字，set可以这是分组啥的
upset_metadata <- data.frame(id = colnames(plot_dat)[1:6],set = colnames(plot_dat)[1:6])
#               id            set
# ChenB_2020_SLE ChenB_2020_SLE
# ChuY_2021_Gout ChuY_2021_Gout
# HuangR_2020_AS HuangR_2020_AS
#   LiuP_2021_MG   LiuP_2021_MG
#  ZhouC_2020_AS  ZhouC_2020_AS
#   ZhuQ_2021_GD   ZhuQ_2021_GD
set.metadata = list(data = upset_metadata, plots = list(
                    list(type = "matrix_rows", column = "set", alpha = .4,
                         colors = c(HuangR_2020_AS = "#e41a1c", ZhouC_2020_AS = "#377eb8", ZhuQ_2021_GD = "#984ea3",
                                   ChuY_2021_Gout = "#FFD320", LiuP_2021_MG = "#ff7f00", ChenB_2020_SLE = "#f781bf"))))

# 给柱子加颜色或者给集合加颜色
# 在邻接矩阵的最后一列标记每个元素的属性，可以在这里根据属性加颜色，例如
#                                    ChenB_2020_SLE ChuY_2021_Gout HuangR_2020_AS LiuP_2021_MG ZhouC_2020_AS ZhuQ_2021_GD direction
# ChenB_2020|SRR8898460_k141_63782_1              0              0              0            0             1            1   Control
# ChenB_2020|SRR8898461_k141_194380               0              0              0            1             1            0      Case
# ChenB_2020|SRR8898463_k141_35192_1              1              1              0            0             0            0   Control
# ChenB_2020|SRR8898465_k141_202808               0              0              0            0             1            1   Control
# ChenB_2020|SRR8898471_k141_135945               0              1              0            0             1            0   Control
# ChenB_2020|SRR8898472_k141_46439_1              0              0              0            0             1            1   Control
# 指定query为elements, params传入一个向量，第一个参数是属性的列名，第二个参数是属性的名称，然后指定颜色，active指定箭头（试一下就知道了） 
queries = list(list(query = elements, params = c("direction","Control"), color = "#74add1", active = T))

# 指定query为intersects，传入的参数要求是列表格式的，输入两个或两个以上的集合名称，然后指定颜色就可以了
queries = list(list(query = intersects, params = list("ChenB_2020_SLE","ChuY_2021_Gout"), color = "white", active = T))

