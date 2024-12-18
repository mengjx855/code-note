# info ----
# Creator：Jinxin Meng
# Created Date: 2022-9-28
# Modified Date: 2022-9-28

# ggpar ----
ggpar(
  p,
  palette = NULL,
  gradient.cols = NULL,
  main = NULL, # 图形标题内容
  submain = NULL, # 图形子标题内容
  caption = NULL, # 图形说明文字
  xlab = NULL, # x轴标题内容
  ylab = NULL, # y轴标题内容
  title = NULL, # 图形标题内容
  subtitle = NULL, # 图形副标题内容
  font.main = NULL, # 标题字体，长度为3的向量，分别表示的字体大小、样式（plain|bold|italic|bold.italic）和颜色。
  font.submain = NULL, # 副标题字体
  font.x = NULL, # xlab字体格式
  font.y = NULL, # ylab字体格式
  font.caption = NULL,
  font.title = NULL,
  font.subtitle = NULL,
  font.family = "",
  font.tickslab = NULL, # 坐标轴字体格式 (size, face, color), e.g.: c(14, "bold", "red"). 
  font.xtickslab = font.tickslab, # 字体 ("plain", "italic", "bold", "bold.italic")
  font.ytickslab = font.tickslab, # 字体 ("plain", "italic", "bold", "bold.italic")
  xlim = NULL,
  ylim = NULL,
  xscale = c("none", "log2", "log10", "sqrt"),
  yscale = c("none", "log2", "log10", "sqrt"),
  format.scale = FALSE,
  legend = NULL, # c("top", "bottom", "left", "right", "none")
  legend.title = NULL, # 图例的标题。c(fill = "", color = "group")
  font.legend = NULL, # 图例的字体
  ticks = TRUE, # 是否加tick
  tickslab = TRUE, # 是否加tick标签
  x.text.angle = NULL, # 坐标轴字的角度
  y.text.angle = NULL, # 坐标轴字的角度
  xtickslab.rt = x.text.angle, # 坐标轴字的角度
  ytickslab.rt = y.text.angle, # 坐标轴字的角度
  xticks.by = NULL,
  yticks.by = NULL,
  rotate = FALSE,
  orientation = c("vertical", "horizontal", "reverse"),
  ggtheme = NULL,
  ...
)

# ggline ----
ggline(
  data,
  x,
  y,
  group = 1,
  numeric.x.axis = FALSE,
  combine = FALSE,
  merge = FALSE,
  color = "black",
  palette = NULL,
  linetype = "solid",
  plot_type = c("b", "l", "p"),
  size = 0.5,
  shape = 19,
  stroke = NULL,
  point.size = size,
  point.color = color,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  facet.by = NULL,
  panel.labs = NULL,
  short.panel.labs = TRUE,
  select = NULL,
  remove = NULL,
  order = NULL,
  add = "none",
  add.params = list(),
  error.plot = "errorbar",
  label = NULL,
  font.label = list(size = 11, color = "black"),
  label.select = NULL,
  repel = FALSE,
  label.rectangle = FALSE,
  show.line.label = FALSE,
  position = "identity",
  ggtheme = theme_pubr(),
  ...
)

# ggscatter ----
ggscatter(
  data,
  x,
  y,
  combine = FALSE,
  merge = FALSE,
  color = "black",
  fill = "lightgray",
  palette = NULL,
  shape = 19,
  size = 2,
  point = TRUE,
  rug = FALSE,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  facet.by = NULL,
  panel.labs = NULL,
  short.panel.labs = TRUE,
  add = c("none", "reg.line", "loess"),
  add.params = list(),
  conf.int = FALSE,
  conf.int.level = 0.95,
  fullrange = FALSE,
  ellipse = FALSE, # 置信区间
  ellipse.level = 0.95,
  ellipse.type = "norm",
  ellipse.alpha = 0.1,
  ellipse.border.remove = FALSE,
  mean.point = FALSE, # 每组样本的均值点
  mean.point.size = ifelse(is.numeric(size), 2 * size, size),
  star.plot = FALSE, # 样本和中点的连线
  star.plot.lty = 1,
  star.plot.lwd = NULL,
  label = NULL,
  font.label = c(12, "plain"),
  font.family = "",
  label.select = NULL,
  repel = FALSE,
  label.rectangle = FALSE,
  parse = FALSE,
  cor.coef = FALSE,
  cor.coeff.args = list(),
  cor.method = "pearson",
  cor.coef.coord = c(NULL, NULL),
  cor.coef.size = 4,
  ggp = NULL,
  show.legend.text = NA,
  ggtheme = theme_pubr(),
  ... # passed to geom_point() and ggpar()
)

# ggscatterhist ----
ggscatterhist(
  data,
  x,
  y,
  group = NULL,
  color = "black",
  fill = NA,
  palette = NULL,
  shape = 19,
  size = 2,
  linetype = "solid",
  bins = 30,
  margin.plot = c("density", "histogram", "boxplot"), # margin区域的图型
  margin.params = list(),
  margin.ggtheme = theme_void(),
  margin.space = FALSE,
  main.plot.size = 2,
  margin.plot.size = 1,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  legend = "top",
  ggtheme = theme_pubr(),
  print = TRUE,
  ... # passed to ggscatter()
)


# ggbarplot ----------------------------------------------------------------------------
ggbarplot(
  data,
  x,
  y,
  combine = FALSE,
  merge = FALSE,
  color = "black",
  fill = "white",
  palette = NULL,
  size = NULL,
  width = NULL,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  facet.by = NULL, # 分面
  panel.labs = NULL,
  short.panel.labs = TRUE,
  select = NULL,
  remove = NULL,
  order = NULL,
  add = "none",
  add.params = list(),
  error.plot = "errorbar",
  label = FALSE, # 是否添加标签
  lab.col = "black", # 标签颜色
  lab.size = 4, # 标签大小
  lab.pos = c("out", "in"), # 标签位置
  lab.vjust = NULL, # 标签垂直移动
  lab.hjust = NULL, # 标签水平移动
  lab.nb.digits = NULL,
  sort.val = c("none", "desc", "asc"),
  sort.by.groups = TRUE,
  top = Inf,
  position = position_stack(),
  ggtheme = theme_pubr(),
  ...
)

# ggerrorplot --------------------------------------------------------------------------------------
ggerrorplot(
  data,
  x,
  y,
  desc_stat = "mean_se",
  numeric.x.axis = FALSE,
  combine = FALSE,
  merge = FALSE,
  color = "black",
  fill = "white",
  palette = NULL,
  size = NULL,
  width = NULL,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  facet.by = NULL,
  panel.labs = NULL,
  short.panel.labs = TRUE,
  select = NULL,
  remove = NULL,
  order = NULL,
  add = "none",
  add.params = list(),
  error.plot = "pointrange",
  ci = 0.95,
  position = position_dodge(),
  ggtheme = theme_pubr(),
  ...
)

# ggboxplot ----
ggboxplot(
  data,
  x,
  y,
  combine = FALSE,
  merge = FALSE,
  color = "black",
  fill = "white",
  palette = NULL,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  bxp.errorbar = FALSE, # 误差线
  bxp.errorbar.width = 0.4, # 误差线宽度
  facet.by = NULL,
  panel.labs = NULL,
  short.panel.labs = TRUE,
  linetype = "solid", # 线型
  size = NULL, # 线宽
  width = 0.7,
  notch = FALSE,
  outlier.shape = NA, # 不要离散点
  select = NULL,
  remove = NULL,
  order = NULL,
  add = "none",
  add.params = list(),
  error.plot = "pointrange",
  label = NULL,
  font.label = list(size = 11, color = "black"),
  label.select = NULL,
  repel = FALSE,
  label.rectangle = FALSE,
  ggtheme = theme_pubr(),
  ...
)


# stat_compare_means ----
stat_compare_means(
  mapping = NULL,
  data = NULL,
  method = NULL,
  paired = FALSE,
  method.args = list(), # 
  ref.group = NULL,
  comparisons = NULL,
  hide.ns = FALSE,
  label.sep = ", ",
  label = NULL,
  label.x.npc = "left",
  label.y.npc = "top",
  label.x = NULL,
  label.y = NULL,
  vjust = 0,
  tip.length = 0.03,
  bracket.size = 0.3,
  step.increase = 0,
  symnum.args = list(), # 	a list of arguments to pass to the function symnum for symbolic number coding of p-values. For example, symnum.args <- list(cutpoints = c(0, 0.0001, 0001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")).
  geom = "text",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
)

# 棒棒糖图
ggdotchart(
  data,
  x,
  y,
  group = NULL,
  combine = FALSE,
  color = "black",
  palette = NULL,
  shape = 19,
  size = NULL,
  dot.size = size,
  sorting = c("ascending", "descending", "none"),
  add = c("none", "segment"), # 是否加线
  add.params = list(), # 线条属性
  x.text.col = TRUE,
  rotate = FALSE,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  facet.by = NULL,
  panel.labs = NULL,
  short.panel.labs = TRUE,
  select = NULL,
  remove = NULL,
  order = NULL,
  label = NULL,
  font.label = list(size = 11, color = "black"),
  label.select = NULL,
  repel = FALSE,
  label.rectangle = FALSE,
  position = "identity",
  ggtheme = theme_pubr(),
  ...
)


# 饼图
# 数据准备，最起码得有一列值
ggpie(
  data,
  x,
  label = x,
  lab.pos = c("out", "in"),
  lab.adjust = 0,
  lab.font = c(4, "plain", "black"),
  font.family = "",
  color = "black",
  fill = "white",
  palette = NULL,
  size = NULL,
  ggtheme = theme_pubr(),
  ...
)

"#58BBD0","#EB9256","#9B99CB","#AD9F2A","#87CCBA","#ACDEF3","#F9C6B3","#F5EAF2","#F9C851","#D3EDE4"

