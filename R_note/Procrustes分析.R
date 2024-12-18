

# Procrustes分析（Procrustes Analysis，普鲁克分析）是一种通过分析形状分布，比较两组数据一致性的方法。
# 数学上来讲，就是不断迭代，寻找标准形状（canonical shape），并利用最小二乘法寻找每个对象形状到这个标准形状的仿射变化方式（Gower, 1975）。
# 该过程也称为最小二乘正交映射（least-squares orthogonal mapping）。

# Procrustes分析的过程
# 简单地说，Procrustes分析基于匹配两个数据集中的对应点（坐标），通过平移、旋转和缩放其中一个数据集中点的坐标以匹配另一数据集中对应点的坐标，并最小化点坐标之间的偏差平方和（表示为M2）。
# 对应点坐标之间的偏差称为矢量残差（vector residuals），小的矢量残差代表了两数据集具有更高的一致性。

# Procrustes分析M2统计量的显著性检验（PROTEST）
# Procrustes分析对两个数据集中点的坐标进行描述性总结和图形化比较，尽管M2统计量提供了对一致性的度量（M2越小表示两数据集关联度越高），但是并未评估M2是否比预期的要好（M2是否显著，而不是由偶然所致）。
# 可通过置换检验的原理实现对M2显著性的检验，称为PROTEST或PROcrustean randomization test。
# 通过在其中一个数据集中随机地对观测值进行置换，同时保持数据集的共变结构，之后使用Procrustes分析计算随机置换后数据集与另一数据集的M2值（称为Mp2）。
# 如此随机置换共N次，并记录原始观测值的M2 < Mp2的次数，由此获得p值：p=(n+1)/(N+1)
# 如果p达到显著性水平（如p<0.05），则可以认为原始观测值的M2并非由偶然因素所致，M2显著小于随机数据集的Mp2，两个数据集的期望值比在随机情况下表现出更大的一致性。

# 例如在组学分析中，经常需要分析来自相同样本不同组学数据集之间是否存在相似性关系，
# 就微生物组而言，例如开篇时提到的那个示例，物种组成丰度和ARGs丰度是否存在潜在一致性。
# 当然分析潜在关系的方法有很多，Procrustes分析就是一种很好的选择。
# 再如群落分析中，常见通过Procrustes分析物种组成与环境的关系，以及物种形态、遗传组成、空间结构、行为特征等更多类型的数据集类似的案例。

# # 提到群落分析，Procrustes分析与Mantel test都是用于分析物种组成和环境属性关系的常见方法，当然二者的具体关注点还是有区别的，方法各有自身的优点。
# 但如果只聚焦在评估两数据集一致性上，似乎Procrustes分析更直观一些。
# M2统计量及其显著性检验p值提供了两个数据集之间一致性的总体度量，同时数据集的图形匹配和相关残差提供了比Mantel test更丰富的信息源。
# 在对应点的坐标匹配度较好时，两个数据集表现出良好的一致性。
# 坐标匹配度越差表明这些点与整体趋势不匹配，这类似于回归分析中残差较大的点，这些点不符合样本的总体趋势。
# 此外，PROTEST的统计功效也被证明优于Mantel test的统计功效（Peres-Neto and Jackson, 2001）。
# 因此，如果两组数据之间存在潜在关系，则Procrustes分析更有能力检测到，并且鉴于结果的图形性质，它还提供了出色的解释性准则。

# chatgpt
# Procrustes 分析的目的是评估两个数据集中对象形状的相似性，通过旋转、平移和缩放来最小化两个配准形状之间的差异。在进行 Procrustes 分析时，可以选择使用对称模式或非对称模式。
# 对称模式 (Symmetric Procrustes Analysis): 这种模式不区分两个数据集中的形状，即它将两个形状视为等同重要，因此在变换过程中两个形状都可以被平移、旋转和缩放。对称 Procrustes 分析常用于形状对称性的评价，或者当两个形状同等重要时，例如在形态学对比研究中。
# 非对称模式 (Asymmetric Procrustes Analysis or Generalized Procrustes Analysis): 这种模式中，一组形状（通常是目标或参照形状）在分析过程中保持不变，而另一组形状会被变换以最小化到参照形状的差异。非对称 Procrustes 分析适用于那些有明确的基准或参考数据集的情况，如将个体形状与某个标准或平均形状进行对比。
# 选择哪种模式取决于研究问题和数据类型。对于一些分析，可能需要尝试两种模式，并根据结果和研究目标来决定哪种方法更适合。

# 作者回应 https://stats.stackexchange.com/questions/563911/what-is-the-difference-between-symmetric-and-non-symmetric-in-procrustes-protest
# 将矩阵B拟合到目标矩阵A：我想获得与A类似的旋转/缩放/平移矩阵B。如果这是您的目标，您应该使用非对称(B到A)旋转。 
# 如果您不想获得结果（旋转），但您只是对两个矩阵之间的相似性以及该相似性的统计量感兴趣，则应该使用对称旋转。 
# 非对称分析类似于回归分析 (B on A），对称 Procrustes 类似于相关分析（A和B之间的相似性，不考虑方向）。

# https://blog.csdn.net/woodcorpse/article/details/106554527

#### info ####
# 生信小白鱼，20231228
setwd("F:/R_proj/Rstat/Procrustes_analysis/")
library(vegan)

##样方-环境属性矩阵
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#环境变量的 PCA 需要标准化，详情 ?rda
env_pca <- rda(env, scale = TRUE)

##样方-物种丰度矩阵
otu <- read.delim('spe_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#物种数据 Hellinger 预转化，详情 ?decostand  
# 在PCA和RDA及K-均一值划分分析前需要对数据进行Hellinger转化。
otu_hel <- decostand(otu, method = 'hellinger')
 
#对转化后的物种数据执行 PCA，无需标准化，详情 ?rda
otu_pca <- rda(otu_hel, scale = FALSE)

##排序图比较，以 PCA 的 I 型标尺为例
par(mfrow = c(1, 2))
biplot(env_pca, choices = c(1, 2), scaling = 1, 
    main = '环境组成的PCA', col = c('red', 'blue'))
biplot(otu_pca, choices = c(1, 2), scaling = 1, 
    main = '物种组成的PCA', col = c('red', 'blue'))

##Procrustes 分析
#提取两个 PCA 中的样方排序坐标，均以 I 型标尺为例
site_env <- summary(env_pca, scaling = 1)$site
site_otu <- summary(otu_pca, scaling = 1)$site

#执行 Procrustes 分析，详情 ?procrustes
#以对称分析为例（symmetric = TRUE）
proc <- procrustes(X = env_pca, Y = otu_pca, symmetric = TRUE)
summary(proc)

#旋转图
plot(proc, kind = 1, type = 'text')

#一些重要的结果提取
names(proc)

head(proc$Yrot)  #Procrustes 分析后 Y 的坐标
head(proc$X)  #Procrustes 分析后 X 的坐标
proc$ss  #偏差平方和 M2 统计量
proc$rotation  #通过该值可获得旋转轴的坐标位置

#残差图
plot(proc, kind = 2)
residuals(proc)  #残差值

#PROTEST 检验，详情 ?protest
#以 999 次置换为例
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(123)
prot <- protest(X = env_pca, Y = otu_pca, permutations = how(nperm = 999))
prot

#重要统计量的提取
names(prot)
prot$signif  #p 值
prot$ss  #偏差平方和 M2 统计量

##ggplot2 作图
library(ggplot2)

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#添加分组信息
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')

#ggplot2 作图
p <- ggplot(Y) +
geom_point(aes(X1, X2, color = groups), size = 1.5, shape = 16) +
geom_point(aes(PC1, PC2, color = groups), size = 1.5, shape = 1) +
scale_color_manual(values = c('red2', 'purple2', 'green3'), limits = c('A', 'B', 'C')) +
geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.1, 'cm')), 
    color = 'blue', size = 0.3) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
    legend.key = element_rect(fill = 'transparent')) +
labs(x = 'Dimension 1', y = 'Dimension 2', color = '') + 
geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) + 
annotate('text', label = sprintf('M^2 == 0.2178'), 
    x = -0.21, y = 0.42, size = 3, parse = TRUE) +
annotate('text', label = 'P < 0.001', 
    x = -0.21, y = 0.38, size = 3, parse = TRUE) 

p

#输出图片
ggsave('procrustes.pdf', p, width = 6, height = 5)
ggsave('procrustes.png', p, width = 6, height = 5)

# 原文链接：https://blog.csdn.net/qq_42830713/article/details/128474690
# 总体的过程是：
# 1、加载R包与数据集
# 2、降维分析
# 在进行普鲁克分析前，需要先对两个数据进行降为分析，这里使用的是NMDS
# 也可以使用PCA或者PCoA，然后进行普鲁克分析，计算显著性
# 普氏分析
# 4、选择模式
# 通过对该函数中各参数的了解，可知X为目标矩阵也就是降维后的环境（功能基因等）坐，Y为降维后的物种数据的坐标，因为后续普氏分析中旋转和缩放操作是针对Y，将Y匹配给X。
# 另外，当symmetric=FALSE时，处于"非对称"模式，X和Y的分配值调换后，普氏分析的偏差平方和（M^2^）也会随之改变。
# 而当symmetric=TRUE时，从而给出更合适比例的对称统计。“对称”模式下，X和Y的分配值调换后，普氏分析的偏差平方和（M2）不会发生改变，但注意旋转仍将是非对称的。
# 5、评价一致性
plot(pro.s.e, kind = 2)
residuals(pro.s.e)
# 6、事后检验999次
# 事后检验，对偏差平方和M2统计量进行置换999次的普氏检验。由于置换999次检验会存在细微的误差，我们设定了种子数，避免重复存在差异。
# 7、提取结果
# Procrustes Anaylsis ----------
# 评估物种群落结构与环境因子间是否具存在显著的相关性（两个数据集）
# 加载R包
library(vegan)
# 加载数据
data(varespec)
head(varespec)
data(varechem)
head(varechem)
# 需要先对数据计算其样本间的距离
# 一般物种群落使用bray-curtis距离，而环境因子使用欧式距离
spe.dist <- vegdist(varespec) # 默认Bray-Curtis
env.dist <- vegdist(scale(varechem), "euclid")
?procrustes
# 在进行普鲁克分析前，需要先对两个数据进行降为分析，这里使用的是NMDS
# 也可以使用PCA或者PCoA，然后进行普鲁克分析，计算显著性
mds.s <- monoMDS(spe.dist)
mds.e <- monoMDS(env.dist)
# 以对称模式为例进行普氏分析（symmetric = TRUE）
pro.s.e <- procrustes(mds.s,mds.e, symmetric = TRUE)
summary(pro.s.e)
# 一致性
plot(pro.s.e, kind = 2)
residuals(pro.s.e)
# 普氏分析中M2统计量的显著性检验
set.seed(1)
pro.s.e_t <- protest(mds.s,mds.e, permutations = 999)
pro.s.e_t
 偏差平方和（M2统计量）
pro.s.e_t$ss
# 对应p值结果
pro.s.e_t$signif
library(ggplot2)
# 获得x和y轴的坐标及旋转过的坐标
Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X))
Pro_X <- data.frame(pro.s.e$rotation)

# 绘图
ggplot(data=Pro_Y) +
  geom_segment(aes(x = X1, y = X2,
               xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2),
               # geom_segment 绘制两点间的直线
               arrow = arrow(length = unit(0, 'cm')),
               color = "#9BBB59", size = 1) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2,
               xend = MDS1, yend = MDS2),
               arrow = arrow(length = unit(0.2, 'cm')),
               color = "#957DB1", size = 1) +
  geom_point(aes(X1, X2), color = "#9BBB59", size = 3, shape = 16) +
  geom_point(aes(MDS1, MDS2), color = "#957DB1", size = 3, shape = 16) +
  theme(panel.grid = element_blank(), # 绘制背景
        panel.background = element_rect(color = 'black',
        fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Correlation between community and environment") +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[1,2]/Pro_X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[2,2]/Pro_X[2,1], size = 0.3) +
  annotate('text', label = 'Procrustes analysis:\n
                    M2 = 0.6297, p-value = 0.001',
           x = -0.3, y = 0.3, size = 4,hjust = 0) +
  theme(plot.title = element_text(size=14,colour = "black",
                      hjust = 0.5,face = "bold"))

library(export)
graph2ppt(file="Procrustes.ppt", append=T, height=5, width=5)


# 自己测试
#### 普鲁克分析 ####
library(vegan)
profile_args <- read.delim("../ARGs/ARGs_tpm", row.names = 1)
profile_mges <- read.delim("../taxa/species_tpm", row.names = 1)
profile_taxa <- read.delim("../taxa/species_tpm", row.names = 1)

# 距离
dist_1 <- vegdist(t(profile_args))
dist_2 <- vegdist(t(profile_mges))

# 降维
PCoA_1 <- cmdscale(dist_1)
PCoA_2 <- cmdscale(dist_2)

# 以对称模式进行普氏分析（symmetric = TRUE）
proc <- procrustes(PCoA_1, PCoA_2, symmetric = T)
summary(proc)

# 评价一致性
# 如果样本中物种与环境一致性（相似性）越近，则对应的残差越小，
# 反之物种与环境的相似性越远，则残差越大（三条辅助线对应的位置分别为残差25%、50%和75%）
plot(proc, kind = 2)
residuals(proc)

# 事后检验999次
# 普氏分析中M2统计量的显著性检验
# 在 Procrustes 分析中，M2 常常指代的是 Procrustes 统计量，它是度量两个形状间差异程度的一个指标。
# 更具体地说，M2 是源数据集的点经过旋转、缩放和/或平移后与目标数据集中对应点的平方距离和。
# 它代表了变换后源数据集的点与目标数据集点之间的不匹配程度。M2 的数值越小，表明两组数据集的形状越相似；反之，则表明它们之间的差异越大。
set.seed(2024)
proc_test <- protest(PCoA_1, PCoA_2, permutations = 999)
proc_test
# Call:
#   protest(X = PCoA_args, Y = PCoA_taxa, permutations = 999) 
# Procrustes Sum of Squares (m12 squared):        0.5268  # M2统计值
# Correlation in a symmetric Procrustes rotation: 0.6879 
# Significance:  0.001  # 结果显示: p = 0.001具有很强显著性
# Permutation: free
# Number of permutations: 999

# 提取结果
#偏差平方和（M2统计量）
proc_test$ss
# 对应p值结果
proc_test$signif

# 基于 ggplot2 包ggplot函数将其结果绘制成图，并利用 export 包导出 ppt 。
# 首先提取降维后的数据轴1和2的坐标，并且提取转换的坐标；然后进行绘制。
library(ggplot2)
# 获得x和y轴的坐标及旋转过的坐标
# Pro_Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
# Pro_X <- data.frame(proc$rotation)
proc_point <- cbind(
  data.frame(proc$Yrot) %>% rename(X1_rotated = X1, X2_rotated = X2),# Y-矩阵旋转后的坐标，就是物种矩阵为了逼近X做了调整，调整后的坐标。
  data.frame(proc$X) %>% rename(X1_target = Dim1, X2_target = Dim2) # X-目标矩阵，就是proc分析输入的第一个矩阵，作为目标矩阵没变化。
  )
proc_coord <- data.frame(proc$rotation) # Y旋转后的坐标轴

# 绘图
ggplot(proc_point) + # 旋转坐标到目的坐标的一半上一个颜色，另一半上不同的颜色
  geom_segment(aes(x = X1_rotated, y = X2_rotated, xend = (X1_rotated + X1_target)/2, yend = (X2_rotated + X2_target)/2), 
               arrow = arrow(length = unit(0, 'cm')), color = "#9BBB59", size = .4) +
  geom_segment(aes(x = (X1_rotated + X1_target)/2, y = (X2_rotated + X2_target)/2, xend = X1_target, yend = X2_target),
               arrow = arrow(length = unit(0.2, 'cm')), color = "#957DB1", size = .4) +
  geom_point(aes(X1_rotated, X2_rotated), color = "#9BBB59", size = 1.6, shape = 16) + # 旋转坐标点
  geom_point(aes(X1_target, X2_target), color = "#957DB1", size = 1.6, shape = 16) + # 目的坐标点
  labs(x = 'Dim 1', y = 'Dim 2',
       subtitle = paste0("coefficients: M2 = ", round(proc_test$ss, 4), ", p = ", proc_test$signif)) +
  labs(title = "Correlation analysis by Procrustes analysis") +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.4) +
  geom_abline(intercept = 0, slope = proc_coord[1,2]/proc_coord[1,1], size = 0.4) +
  geom_abline(intercept = 0, slope = proc_coord[2,2]/proc_coord[2,1], size = 0.4) +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.line = element_blank(),
        plot.title = element_text(size = 10, color = "black"),
        plot.subtitle = element_text(size = 10, color = "black"),
        panel.border = element_rect(linewidth = .4, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 3/4)
