# Encoding: utf-8
# Creator：Jinxin Meng
# Created Date:2022-1-19
# Modified Date:2022-1-19

# 添加百分比轨道
circos.track(track.index = 2, # 轨道顺序由外向内，最外圈是1，中间为2，最内为3
             panel.fun = function(x, y) {
               xlim = get.cell.meta.data("xlim")
               ylim = get.cell.meta.data("ylim")
               sector.name = get.cell.meta.data("sector.index")
               xplot = get.cell.meta.data("xplot")
               circos.lines(x = xlim, 
                            y = c(mean(ylim), mean(ylim)), 
                            lty = 3) # 添加点线
               by = 0.5
               # 绘制百分比刻度标签
               for (p in seq(by, 1, by = by)) {
                 circos.text(p * (xlim[2] - xlim[1]) + xlim[1],
                             mean(ylim) + 0.1,
                             paste0(p * 100, "%"),
                             cex = 0.4,
                             adj = c(0.5, 0),
                             niceFacing = TRUE)
               }
             }, bg.border = NA
)

# df : 至少有两列的数据帧。前两列指定连接，第三列（可选）包含映射到链接宽度的数值以及颜色（如果将col指定为颜色映射）功能地块中的扇区将为并集（df$$1￥，df$$2￥）。
# grid.col : 对应于扇区的网格颜色。向量的长度应为1或部门。It他更喜欢这样网格.col名称对应的命名向量部门。如果它不是一个命名向量，而是网格.col对应于扇区的顺序。
# grid.border : 网格的边框。如果为空，则边框颜色与网格颜色相同
# transparency : 链接颜色的透明度，0表示没有透明度，1表示完全透明透明度。如果透明度已在列或中设置行.列或者列.col，此参数将已忽略。NAalso忽略此参数。
# col : 链接的颜色。它可以是与df中的连接相对应的向量，也可以是根据df中的值（第三列）生成颜色的函数，也可以是表示所有链接的颜色相同的单个值。可以使用colorRamp2生成一个函数，将值映射到颜色。
# order : 扇区顺序。默认顺序为union（df$$1￥，df$$2￥）。
# directional : 链接是否有方向。1表示方向从df中的第一列到第二列，-1表示相反方向，0表示无方向，2表示双向。该值可以是与df中的行数长度相同的向量。
# xmax : x轴上的最大值，该值应为命名向量。
# direction.type : 表示方向的类型。可以是“diffHeight”和“arrows”中的一个或两个值。如果该值包含“diffHeight”，则使用链接的不同高度来表示起始根具有较长高度的方向，以使人们感觉到有东西在发生。如果该值包含“arrows”，则用户可以使用以下命令自定义箭头争论。争论值可以是与df中的行数具有相同长度的向量。注意：如果要为某些链接同时设置diffheight和arrows，则需要将这两个选项嵌入一个字符串中，例如“高度+箭头".
# diffHeight : 如果directive设置为TRUE，则两个“根”之间的高度差。如果该值设置为正值，则起始根比结束根短；如果该值设置为负值，则起始根比结束根长。该值可以是与df中的行数长度相同的向量。
# link.target.prop : 如果弦图是方向图，对于每个源扇区，是否绘制显示目标扇区比例的条形图。
# target.prop.height : 钢筋的高度link.target.prop链接已打开。
# reduce : 如果某个网格的宽度与整个圆的宽度之比小于此值，则会在上删除该网格绘图集如果你想保持所有的小网格，它的值应该小于零。
# self.link : 如果一个扇区中存在自链接，1表示链接将退化为“山”，宽度与此连接的值相对应。2表示起始根和结束根的宽度都与连接的值相对应。
# preAllocateTracks : 在绘制弦图之前预先分配空音轨。它可以是一个数字，表示需要创建多少空磁道，也可以是一个包含空磁道设置的列表。有关详细信息，请参阅渐晕图。
# annotationTrack : 应该绘制哪个注释轨迹？默认情况下，将创建包含扇区名称的轨迹和包含栅格的轨迹。
# annotationTrackHeight : 轨迹高度与annotationTrack中的值相对应。
# link.border : 链接的边框，单个标量或与df或数据帧的nrows长度相同的向量
# link.lwd : 链接边界的宽度，单个标量或与df或数据帧的nrows长度相同的向量
# link.lty : 链接边框样式，单个标量或与df或数据帧的nrows长度相同的向量
# link.auto : 忽略。
# link.sort : 是否根据每个扇区上链接的宽度对其进行排序。如果设置为“总体”，则所有链接都将按顺序排序，无论它们来自第一列还是第二列。
# link.decreasing : 为了链接.排序
# link.arr.length : 传给循环链接. 此参数的格式与链接.lwd.
# link.arr.width : 传给箭头。此参数的格式与链接.lwd.
# link.arr.type : 传给循环链接，设置与相同链接.lwd. 默认值为三角形。
# link.arr.col : 颜色或放在皮带中心的单线连接。此参数的格式与链接.lwd.
# link.arr.lwd : 位于皮带中心的单线连接的线宽。此参数的格式与链接.lwd.
# link.arr.lty : 位于皮带中心的单线链环的线型。此参数的格式与链接.lwd.
# link.largest.ontop : 控制添加链接的顺序，是否基于绝对值？
# link.rank : 此参数已删除。
# link.visible : 是否绘制链接。该值是逻辑值，如果设置为FALSE，则不会绘制相应的链接，但仍会占用空间。此参数的格式与链接.lwd
# link.zindex : 为了添加到圆的链接，一个大的值意味着以后添加它。
# link.overlap : 如果是方向弦图，那么在同一扇区中出现或结束的链路是否重叠？
# scale : 将每个扇区缩放到相同的宽度
# group : 它包含组标签，扇区名称用作向量中的名称。
# big.gap : 数据框第一列扇区和第二列扇区之间的间距。
# small.gap : 行业之间的差距很小。