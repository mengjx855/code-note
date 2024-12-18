# Jinxin Meng, 20220906, 20220911 ------------------
setwd("/proj/proj_2022/20220906_virome_ACVD_Wuhan/PCoA/")
pacman::p_load(dplyr, tidyr, tibble)
library(ggpubr)
library(vegan)

# betadisper 评估数据离散程度 判断PREMANOVA的结果 ----------------------
# ref: https://blog.csdn.net/qazplm12_3/article/details/120559125
# ref: https://cloud.tencent.com/developer/article/1761255
# ref: https://www.jianshu.com/p/2e74b7df023c
# 使用函数betadisper()测试组同质性
# 差异分析判断不同分组样品检测指标的离散度(方差)差异显著性。
# 如果差异不显著，那么表明分组指标是影响分组数据的主要原因，可以解释adonis的结果。
# betadisper检验Pr(>F)<0.05表明不同组的数据在空间分布的离散度显著不同。这是导致adonis结果显著的主要原因。
# 不同分组之间物种的构成的显著不同不是体现在物种空间中心点的变化，而是物种空间离散度的变化。
distance <- readRDS("distance_vOTUs.rds")
group <- read.delim("../data/sample_group.txt")
group <- group[match(rownames(as.matrix(distance)), group$sample),] %>% 
  as.data.frame(row.names = NULL)

# 分散度评估
dispersion <- betadisper(distance, group$group)
boxplot(dispersion, main = "vOTU level\nANOVA p<0.001\npermutest p<0.001\nTukeyHSD p<0.001")

# 组间分散度差异分析，
anova(dispersion)
permutest(dispersion)
TukeyHSD(dispersion)

# Age
dispersion <- betadisper(distance, group$Gender)
boxplot(dispersion)
anova(dispersion)
permutest(dispersion)
TukeyHSD(dispersion)

# family
otu <- readRDS("../data/vOTUs_rarefied.family.tpm.rds")
distance <- vegdist(t(otu), method = "bray")
group <- read.delim("../data/sample_group.txt")
group <- group[match(rownames(as.matrix(distance)), group$sample),] %>% 
  as.data.frame(row.names = NULL)

# 分散度评估
dispersion <- betadisper(distance, group$group)

boxplot(dispersion, main = "family level\nANOVA p=0.506\npermutest p=0.514\nTukeyHSD p=0.506")

# 组间分散度差异分析，

anova(dispersion)

permutest(dispersion)

TukeyHSD(dispersion)

dispersion$distances %>% 
  data.frame(value = .) %>% 
  rownames_to_column("sample") %>%
  left_join(group, by = "sample") %>% 
  ggboxplot("group", "value", fill = "group", xlab = "", ylab = "Distance to centroid", 
            legend = "none")
