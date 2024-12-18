# Encoding: utf-8
# Creator：Jinxin Meng
# Created Date:2022-7-19
# Modified Date:2023-12-10

library(ComplexHeatmap)

# 所有函数
pheatmap(as.matrix(otu),  # 热图，保存8-4
        cluster_rows = F, cluster_cols = F, # 是否聚类
        cluster_row_slices = F, cluster_column_slices = F, # 分块内部进行聚类，最后分块聚类的顺序是指定的
        treeheight_col = 30, treeheight_row = 30, # 聚类树高
        show_row_dend = F, show_column_dend = F, # 聚类但是不显示列树状图，Heatmap包中的参数，完全不显示聚类树
        breaks = c(-4, -2, 0, 2, 4),  # 设置区间,这个和颜色呼应
        color = c("#f46d43","#fee08b","#ffffff","#d9ef8b","#66bd63"),  # 设置区间颜色
        legend_break = c(-4, -2, 0, 2, 4),
        legend_label = c(-10, -5, 0, 5, 10),
        heatmap_legend_param = list(title = "Relative Abundance (%)"), # legend 标题
        border_color = "white",
        border_gp = gpar(col = "black"), # 分块聚类后边框
        show_rownames = T, show_colnames = T, # 是否显示名字
        fontsize_col = 15, fontsize_row = 15, # 显示名称字体大小
        cellwidth = 20, cellheight = 20, # 格子的大小
        split = row_split, column_split = col_split,  # 分块聚类
        annotation_row = row_annotation, # 注释 
        annotation_colors = colors, 
        annotation_names_col = T,
        row_title_rot = 0, # 分块的标题方向，分块标题的方向的参数是row/columm_title_*
        row_title_side = "right", # 指定行列标签的位置
        row_title_gp = gpar(fontsize = 13.2), # 标题的大小
        column_title_rot = 0,
        column_title_side = c("top", "bottom"),
        column_title_gp = gpar(fontsize = 13.2),
        display_numbers = t(plot_mat_p),
        row_gap = unit(2, "mm"), # 设置行分块之间的间距
        column_gap = unit(1, "mm"),
        use_raster = F, # 行超过2000，保存成了位图了，用此参数禁止格式转变。
        heatmap_legend_param = list(border = "black",  # 修改图例
                                    title = "Distance", 
                                    title_gp = grid::gpar(fontface = "italic", fontsize = 10),
                                    title_position = "topcenter",
                                    legend_direction = "vertical",
                                    legend_width = unit(4, "cm"),
                                    labels_gp = grid::gpar(fontsize = 8)),
)

# 按照行来分块聚类
# 指定和行名顺序一样的一个向量，然后转为因子类型。
split_row <- factor(row_annotation$disease)
split_column <- factor(col_annotation$study)

# 行和列的注释表
row_annotation <- data.frame(disease = group$disease, group = group$group, row.names = group$sample)
#     disease group   列为分类变量，行名为样本名，分类变量指定样本的类型
# s1  gout    control
# s2  gout    control
# s3  gout    disease

# legend_breaks指定图例的颜色区间，与legend_lable一起使用。
# breaks指定颜色的区间，与color对应，每个区间设置颜色。

# color指定的多种方法, 反正就是输入颜色序列，要是指定了breaks，序列的长度要和区间的数量一样
color = colorRampPalette(c("#F3DA7D","#A13043"))(100) # 渐变色 黄-褐
color = colorRampPalette(c("#3288bd", "#ffffff", "#d53e4f"))(100) # 渐变色 蓝-白-红
color = colorRampPalette(c("#a6611a","#dfc27d","#f5f5f5","#80cdc1","#018571"))(100) # 渐变色 褐-白-青
color = colorRampPalette(c("#e66101","#fdb863","#f7f7f7","#b2abd2","#5e3c99"))(100) # 渐变色 黄-白-紫
color = colorRampPalette(c("#ca0020","#f4a582","#f7f7f7","#92c5de","#0571b0"))(100) # 渐变色 红-白-蓝
color = colorRampPalette(c("#d7191c","#fdae61","#ffffbf","#abd9e9","#2c7bb6"))(100) # 渐变色 红-黄-蓝
color = colorRampPalette(c("#f46d43","#fee08b","#ffffff","#d9ef8b","#66bd63"))(100) # 渐变色 红-白-绿
color = paletteer::paletteer_c("grDevices::Temps", 30)
color = circlize::colorRamp2(c(-4, 0, 4), c("#3288bd","white","#d53e4f")) # 这个颜色不会被弱化，可以同时设置break
color = c("#d53e4f","#fc8d59","#fee08b","#ffffbf","#e6f598","#99d594","#3288bd") # 渐变色 红-黄-蓝
color = c("#c51b7d","#e9a3c9","#fde0ef","#f7f7f7","#e6f5d0","#a1d76a","#4d9221") # 渐变色 紫-白-绿
color = c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac") # 渐变色 红-白-蓝
color = c("#8c510a","#d8b365","#f6e8c3","#f5f5f5","#c7eae5","#5ab4ac","#01665e") # 渐变色 褐-白-青

# 指定分组的颜色
colors <- list(group = c(d1 = "#8dd3c7",d3 = "#ffffb3", d5 = "#bebada",d7 = "#fb8072",d14 = "#80b1d3"))

# 设置行列注释并配色
annotation_row <- data.frame(row.names = core) %>% mutate(class = dat$class[match(rownames(.), dat$gene)])
annotation_col <- data.frame(row.names = colnames(plotdat)) %>% mutate(group = group$group[match(rownames(.), group$sample)])
class_color <- structure(c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494"), names = levels(factor(annotation_row$class)))
group_color <- structure(c("#729ECE","#FF9E4A","#67BF5C","#ED665D","#AD8BC9", "#A8786E","#ED97CA","#CDCC5D","#A2A2A2"), names = levels(factor(annotation_col$group)))
color_list <- list(class = class_color, group = group_color)

# display_number 可以直接输入一个与元数据一致的矩阵
# 也可直接对p值矩阵做判断，代码如下：
matrix(ifelse(cor_p_spread <= 0.01, "**", ifelse(cor_p_spread < 0.05 & cor_p_spread > 0.01, "*", "")),nrow(cor_p_spread))

# 给边框上色 
border = "black"  # 这个指定的是每个小方框（Cell）的边框颜色
border_gp = gpar(col = "black")   # 可以指定热图整体边框的颜色

# 按照行或者列分开聚类后，调整间隔距离
row_gap = unit(0,'mm'), column_gap = unit(0,'mm'), 

# 复杂热图添加柱状图
# 热图的原始表格如下
# name  s1  s2  s3
# f1    10  20  30
# f2    23  21  18
# f3    15  13  19
# 准备一个这样的数据表名称为dat_col_annotation，例如，我想给每一个feature添加计算好的LDA值
# name  LDA   color 
# f1    2.10  #729ece
# f2    2.45  #ed665d
# f3    2.83  #67bf5c
# 制作附图的对象
row_bar <- HeatmapAnnotation(LDA = anno_barplot(dat_col_annotation$LDA, gp = gpar(fill = dat_col_annotation$color, col = NA)))
# 绘制主图，把注释的图的对象加入进来，加入注释图的参数包括right_annotation/top_annotation/bottom_annotation/left_annotation
pheatmap(plotdat, right_annotation = col_bar)

# Heatmap()
Heatmap(profile, 
        width, # 热图主体的宽度
        height, # 热图主体的高度
        heatmap_width, # 整个热图的宽度（包括注释）
        heatmap_height, # 整个热图的宽度（包括注释）
        show_row_names = F, # 是否展示行名
        show_column_names = F, # 是否展示列名
        row_split = split_row, # 按照行分块
        row_title_side = "right", # 行分块标题在左还是右
        row_title_gp = gpar(fontsize = 8), # 行标题的格式
        row_title_rot = 0, # 行标题旋转的角度，0就是水平
        row_gap = unit(2, "mm"), # 设置行分块之间的间距
        show_row_dend = F, # 聚类但是不显示行树状图 
        show_column_dend = F, # 聚类但是不显示列树状图 
        )
# 添加简单注释柱状图
bar <- HeatmapAnnotation(bar = dat$status[match(colnames(profile), dat$name)], # 对应列的分组向量
                         col = list(bar = c("depleted" = "#327cc0", "enriched" = "#e43589")), # 分组颜色
                         name = "bar", # 简单注释的名称
                         show_legend = F, # 不展示图例 
                         gp = gpar(fontsize = 8), # 设置标签的格式
                         annotation_legend_param = list(), # 设置图例的格式
                         show_annotation_name = F, # 不展示简单柱状的标记
                         border = T, # 是否展示简单注释边框
                         gap = unit(), # 简单注释和主图的距离
                         annotation_name_gp = gpar(), # 简单注释名称字体格式
                         annotation_name_offset = unit(), # 简单注释与柱子的距离
                         annotation_name_side = "right", # 简单注释的左右位置
                         annotation_name_rot = 0, # 简单注释是否旋转
                         annotation_name_align = T, # 简单注释是否和行列注释对齐
                         annotation_height = 10, # 简单注释的高度
          )

######################################
# 我的复杂热图的思路如下
# 原始表格 dat
# name  s1  s2  s3  s4
# f1    10  20  30  20
# f2    23  21  18  22
# f3    15  13  19  21
# f4    23  24  19  24
# feature的信息表格 feature_info
# name  class
# f1    A
# f2    A
# f3    B
# f4    B
# map表格 sample_group
# sample  group
# s1      g1
# s2      g1
# s3      g2
# s4      g2
# 注释数据，就拿LDA分析作为代表，我计算了每个feature在两个组之间的LDA得分 lefse_res
# name  LDA   pval  enriched
# f1    2.10  0.02  g1
# f2    2.45  0.001 g1
# f3    2.83  0.03  g2
# f4    2.33  0.01  g2
# 接下来进行数据处理
# 颜色设定
group_color <- c("#729ece","#ff9e4a")
names(group_color) <- c("g1","g2")
class_color <- c("#6ee2ff","#f7c530")
names(class_color) <- c("A","B")
# 注释数据
row_annotation_dat <- data.frame(name = rownames(dat)) %>% 
  mutate(class = feature_info$class[match(name, feature_info$name)],
         class_color = class_color[match(class, names(class_color))]) %>% 
  merge(., lefse_res, by = "name", all.x = T) %>% 
  mutate(enriched_color = group_color[match(enriched, names(group_color))])
# 数据长这样
# name  LDA   pval  enriched  class class_color enriched_color
# f1    2.10  0.02  g1        A     #6ee2ff     #729ece
# f2    2.45  0.001 g1        A     #6ee2ff     #729ece
# f3    2.83  0.03  g2        B     #f7c530     #ff9e4a
# f4    2.33  0.01  g2        B     #f7c530     #ff9e4a
col_annotation_dat <- data.frame(name = colnames(dat)) %>% 
  mutate(group = map$group[match(name, map$sample)],
         group_color = group_color[match(group, names(group_color))])
# 数据长这样
# sample  group group_color
# s1      g1    #729ece
# s2      g1    #729ece
# s3      g2    #ff9e4a
# s4      g2    #ff9e4a
# 行列注释，然后给颜色
annotation_row <- select(row_annotation_dat, name, class) %>% column_to_rownames(var = "name")
annotation_col <- select(col_annotation_dat, name, group) %>% column_to_rownames(var = "name")
annotation_colors <- list(class = select(row_annotation_dat, class, class_color) %>% 
                            unique() %>% pull(name = class), 
                          group = select(col_annotation_dat, group, group_color) %>% 
                            unique() %>% pull(name = group))
# 分块
col_split <- factor(col_annotation_dat$group)
row_split <- factor(row_annotation_dat$class)
# 注释图对象的生成
col_bar <- HeatmapAnnotation(LDA = anno_barplot(row_annotation_dat$LDA, gp = gpar(fill = row_annotation_dat$enriched_color, col = NA)))
# 堆叠柱状图
plot_bar <- columnAnnotation(bar = anno_barplot(plot_bar, gp = gpar(fill = 2:4)))
# 最终热图
pheatmap(dat, scale = "none", cluster_rows = T, cluster_cols = T, treeheight_row = 15, treeheight_col = 15,
        show_rownames = F, show_colnames = F, cellwidth = 4, cellheight = 4, color = color, border_color = "#ffffff",
        row_split = row_split, column_split = col_split, column_gap = unit(0.5, "mm"), row_gap = unit(0.5, "mm"),
        annotation_row = annotation_row, annotation_colors = annotation_colors, top_annotation = col_bar)
