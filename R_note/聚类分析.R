# Encoding: utf-8
# Creator：Jinxin Meng
# Created Date:2021-08-10
# Modified Date:2023-11-15

#### kmeans ####
# 利用stats包中的kmeans函数就能实现kmeans聚类方法
# cluster:返回的聚类判别结果， centers：最终的聚类中心，totss：总的方差，withinss：各组组内方差
# tot.withinss：总的组内方差，betweenss：各组的组间方差，size：各组的样本数
# 聚类结果的好坏可以通过组内方差与组间方差、总方差的关系来评价。其目的是尽可能的保证聚类后，各组内的数据进来同质、类与类间的数据尽量异质。
# 因此，聚类后的数据，组内数据方差要尽可能的小，组间的方差要尽可能的大。也就是tot.withinss越小越好、betweenss(totss-tot.withinss)越大越好。
# 进一步的，可以写作F＝betweenss/tot.withinss，当满足F值越大，则可以认为聚类效果越好。
# nstart参数的理解：
# 当我们使用kmeans()函数进行k均值聚类时，算法在开始时需要选择一些初始点作为聚类中心。这些初始点的选择可能会影响最终的聚类结果。为了得到更好的聚类结果，可以多次尝试不同的初始点。
# nstart参数就是用来指定要尝试的不同初始点的个数。默认情况下，nstart的值是1，也就是只尝试一组初始点。如果你希望增加尝试的次数，可以将nstart设置为一个大于1的整数。
# 为什么要尝试多组初始点呢？因为k均值聚类算法是基于随机选择初始点的，而随机性可能导致结果不够稳定。通过多次尝试不同的初始点，我们可以增加找到全局最优解（最好的聚类结果）的机会，并减少陷入局部最优解（较差的聚类结果）的可能性。
set.seed(2023)
x <- rbind(matrix(rnorm(30, sd = 0.3), ncol = 2), matrix(rnorm(30, mean = 1, sd = 0.3), ncol = 2)) # 随机生成值
colnames(x) <- c("x", "y")
km <- kmeans(x, centers = 2, iter.max = 10, nstart = 10) # 聚类分析 centers设置聚类的类群个数，nstart设置起始点的个数， 输入的数据行时样本，列时特征
plot(x, col = km$cluster) # 聚类结果

# 如何选择最佳聚类的类中心个数
# 利用factoextra包的fviz_nbclust() 函数寻找最佳的聚类数
# factoextra包用于提取并可视化多变量数据分析的结果
library(factoextra)
dat <- iris[-5]
fviz_nbclust(dat, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2) # WSS这里的wss(within-cluster sum of squares)是组内平方和
# 在聚类中心为3时，组内平方和降低最多，接下来组内平方和变化不是很大了，这个时候的类中心个数被选择聚类分析的聚类中心数
km <- kmeans(dat, centers = 3, nstart = 20, iter.max = 20) # 执行聚类
# 可视化聚类结果
fviz_cluster(km, data = dat)
# 计算聚类中心
aggregate(dat, by = list(cluster = km$cluster), mean)


#### 层次聚类 ####
# 层次聚类，也称为层次聚类分析，是一种无监督学习方法，它能够根据数据点之间的相似度或距离来将数据分组成一个层次结构的簇。
# 一般来说，对于小样本数据(n<100)，可以考虑使用层次聚类方法，多选择欧氏距离或类平均法。
# dist 函数是R语言中用于计算两个或多个向量之间的距离的函数，method 参数用于指定计算距离所使用的距离度量方法。
# dist() 计算距离的方法：
# "euclidean" ：欧氏距离。它是最常见的距离度量方法，计算数据点之间的直线距离。适用于连续型数据。（主要）
# "maximum"：切比雪夫距离（Chebyshev Distance）。它计算两个数据点之间各维度上的最大差值。适用于连续型数据。
# "manhattan" ：曼哈顿距离。它计算两个数据点之间沿着坐标轴的距离总和。适用于连续型数据。
# "canberra" ：坎贝拉距离（Canberra Distance）。它是一种适用于连续型数据的距离度量方法，特别适用于具有不同尺度的特征。它计算两个数据点之间各维度上的加权距离，其中权重由数据点的值决定。
# "binary" ：二进制距离。它用于处理二进制数据（0和1）。对于两个二进制数据点，二进制距离等于它们不同的位数的个数。
# "minkowski" ：闵可夫斯基距离。它是一种通用的距离度量方法，可以根据参数来调整度量方式。参数 p 决定了闵可夫斯基距离的类型，当 p = 1 时，变成曼哈顿距离，当 p = 2 时，变成欧氏距离。
dat <- iris[-5]
dist <- dist(dat, method = "euclidean")
# hclust 是R语言中用于层次聚类的函数，它是“层次聚类”（Hierarchical Clustering）的缩写。层次聚类是一种无监督学习方法，用于将数据点分组成层次结构的集群或簇。它将数据点逐渐合并成越来越大的集群，直到最终形成一个包含所有数据点的集群。
# hclust 函数的工作方式如下：
# 计算数据点之间的相似度（距离）：首先，通过选择一种距离度量方法（例如欧几里德距离、曼哈顿距离等），计算数据集中每对数据点之间的相似度或距离。这些距离值通常被组织成一个称为距离矩阵的矩阵。
# 创建初始聚类：将每个数据点视为一个单独的聚类，初始时有与数据点数量相同的聚类。
# 合并最近的聚类：从初始聚类开始，选择距离最近的两个聚类，并将它们合并为一个新的聚类。这个过程不断迭代，每次选择距离最近的两个聚类合并，直到所有数据点都被合并到一个聚类为止。
# 构建层次树：将合并的过程记录下来，形成一个树状结构，称为层次聚类树或谱系树。这个树状结构显示了聚类的层次结构，可以用于不同层次的集群分析。
# 剪枝：通过指定一个高度阈值或聚类数目，可以从树的不同层次截取出最终的聚类结果。这称为剪枝操作，它决定了将数据分成多少个最终的簇。
# hclust 函数中的 method 参数用于指定层次聚类中的凝聚方法。不同的凝聚方法会影响层次聚类的结果，因为它们决定了如何计算簇之间的距离和合并方式。以下是 method 参数的各个选项及其含义：
# "ward.D" ：使用Ward's最小方差法。这种方法尝试最小化簇内的方差增加量。它在合并簇时倾向于选择方差增加最小的簇，因此可以产生相对平衡的簇。
# "ward.D2" ：与 "ward.D" 类似，但是使用了方差的平方来计算距离。
# "single" ：单链接法（single linkage）。合并两个簇时，计算两个簇中最近的两个数据点之间的距离。这种方法可能会导致链状效应，形成长而不平衡的簇。
# "complete" ：完全链接法（complete linkage）。合并两个簇时，计算两个簇中最远的两个数据点之间的距离。这种方法倾向于产生更紧凑且球形的簇。
# "average" ：平均链接法（average linkage），也称为UPGMA（Unweighted Pair Group Method with Arithmetic Mean）。合并两个簇时，计算两个簇中所有数据点之间的平均距离。它在某些情况下可以产生相对平衡的簇。
# "mcquitty" ：WPGMA（Weighted Pair Group Method with Arithmetic Mean）。与 "average" 类似，但使用了加权平均距离。
# "median" ：WPGMC（Weighted Pair Group Method with Centroid Mean）。与 "average" 类似，但使用了加权平均距离。
# "centroid" ：UPGMC（Unweighted Pair Group Method with Centroid Mean）。合并两个簇时，计算两个簇的质心之间的距离。这种方法可以产生平衡的簇。
hclust <- hclust(dist, method = "average")
# 展示聚类树状图
par(cex = 0.5)  # 设置文字大小
plot(hclust, main = "hclust", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, hang = -1)
abline(h = 3, col = "red") # 在图上绘制红线
clust = cutree(hclust, h = 3) # 确定红线下的聚类，裁剪树结果
table(clust) # 显示剪枝结果
# ggplot风格的树状图,使用ggdendro包, dendrology 树木学
library(ggdendro)
plotdat <- dendro_data(hclust, type = "rectangle")
plotdat$labels <- mutate(plotdat$labels, species = iris$Species[match(label, rownames(iris))]) 
ggplot() +
  geom_segment(data = plotdat$segments, aes(x, y, xend = xend, yend = yend), lwd = .4) +
  geom_text(data = plotdat$labels, aes(x, y, label = label, color = species), size = 1.5,
            nudge_y = -.1) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

#### PAM 聚类分析 ####
# partitioning around medoids (PAM) clustering method
# PAM（Partitioning Around Medoids）是一种用于聚类分析的算法，类似于K-means，但具有更高的稳健性。PAM通过选择数据集中实际存在的点作为中心（称为medoids），来代表每个簇，而不是计算均值。这使得PAM在处理噪声和离群点方面比K-means更为有效。
# PAM的基本步骤：初始化：选择k个初始medoids。一般是随机选择或通过启发式方法选择。分配对象到最近的medoid：将每个数据点分配到离它最近的medoid所在的簇中。更新medoids：迭代地尝试替换medoids，以找到能使聚类代价（通常是所有数据点到其对应medoid的距离之和）最小化的medoid。迭代：重复步骤2和3，直到medoids不再变化，或聚类代价不再显著下降。
# PAM的优点：稳健性：由于使用实际数据点作为medoids，PAM对离群点不太敏感。解释性：medoids是实际数据点，因此更容易解释每个簇的代表性。
# PAM的缺点：计算复杂度：由于需要计算所有数据点与所有medoids之间的距离，PAM在处理大规模数据集时可能较慢，计算复杂度为 O(k(n−k)^2)，其中n是数据点总数。
# 应用：PAM适用于需要稳健性的聚类任务，尤其是那些可能包含噪声和离群点的数据集。它在生物信息学、市场研究和其他领域都得到了应用。为了提升PAM在大数据集上的效率，可以使用改进算法如CLARA（Clustering Large Applications）和CLARANS（Clustering Large Applications based on Randomized Search）。这些方法通过抽样和随机化技术来减少计算量。
# https://blog.csdn.net/liujh845633242/article/details/103665123 python版本的计算方法

# 安装并加载cluster包
install.packages("cluster")
library(cluster)
# 生成示例数据
set.seed(123)
data <- matrix(rnorm(100), ncol=2)
# 进行PAM聚类，指定簇的个数为3
pam_result <- pam(data, k=3)
# 查看聚类结果
print(pam_result$clustering)
plot(data, col=pam_result$clustering)
points(pam_result$medoids, col=1:3, pch=8, cex=2)

# 在Python中，PAM聚类可以使用scikit-learn-extra库中的KMedoids类来实现。以下是一个Python示例
# 安装scikit-learn-extra库
!pip install scikit-learn-extra
from sklearn_extra.cluster import KMedoids
import numpy as np
import matplotlib.pyplot as plt
# 生成示例数据
np.random.seed(123)
data = np.random.randn(50, 2)
# 进行PAM聚类，指定簇的个数为3
kmedoids = KMedoids(n_clusters=3, random_state=0).fit(data)
# 查看聚类结果
labels = kmedoids.labels_
medoids = kmedoids.cluster_centers_
# 绘制结果
plt.scatter(data[:, 0], data[:, 1], c=labels, cmap='viridis')
plt.scatter(medoids[:, 0], medoids[:, 1], c='red', marker='x')
plt.show()