#### Info ####
# Jinxin Meng, 20231228, 20240126
# Version: 1.0
# mantel分析|ggcor|linkET|两个矩阵的相关性分析
# linkET除了做mantel,还可以做简单的相关性分析

library(linkET)
library(vegan)
data("varespec")
data("varechem")

# 数据介绍
# 物种丰度表
# Callvulg Empenigr Rhodtome Vaccmyrt Vaccviti  ...
# 18     0.55    11.13     0.00     0.00    17.80
# 15     0.67     0.17     0.00     0.35    12.13
# 24     0.10     1.55     0.00     0.00    13.47
# 27     0.00    15.13     2.42     5.92    15.97
# 23     0.00    12.68     0.00     0.00    23.73
# ...
# 环境因子
# N    P     K    Ca    Mg  ...
# 18 19.8 42.1 139.9 519.4  90.0
# 15 13.4 39.1 167.3 356.7  70.7
# 24 20.2 67.7 207.1 973.3 209.1
# 27 20.6 60.8 233.7 834.0 127.2
# 23 23.8 54.5 180.6 777.0 125.8
# ...


# 相关性热图，两个矩阵
correlate(varespec[1:14], varechem) %>% 
  qcorrplot() +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))

# 相关性热图，单个矩阵，半张图
qcorrplot(varespec[1:14], type = "lower") +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))

# 设置热图的过渡色板
set_corrplot_style()
qcorrplot(mtcars) + geom_square()

# mantel test
# 物种数据，环境数据，spec_select和env_select取子矩阵，是子矩阵在做分析
# spec_select参数设置为 x:y列是某个分类单元。例如，一个otu表中，前两个分类单元是OTU1,OTU2是乳酸菌属，可以写成Lac=1:2
mantel_test <- mantel_test(varespec, varechem, spec_select = list(Spec01 = 1:7, Spec02 = 8:18, Spec03 = 19:37, Spec04 = 38:44))
# 检验后，在按照区间划分，后续需要设置属性
mantel_df <- mantel_test %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

qcorrplot(correlate(varechem), type = "upper", diag = F) +
  geom_square() +
  geom_couple(data = mantel_df, aes(colour = pd, size = rd), curvature = nice_curvature()) +
  scale_fill_viridis_c() +
  scale_size_manual(values = c(0.6, 0.9, 1.2)) +
  scale_colour_manual(values = c("#b1bb52", "#9794a4", "#aaaba6")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
ggsave("mantel_test_linkET.pdf", width = 10, height = 8)

