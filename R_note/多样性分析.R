# 一、β-多样性分析介绍
# β（Beta）Diversity：是对不同样品/不同组间样品的微生物群落构成进行比较分析。
# β多样性分析前的数据“来源”：OTUs的丰度信息表；OTUs之间的系统发生关系，
# 计算UnweightedUnifrac及Weighted Unifrac距离。通过多变量统计学方法主成分分析(PCA， Principal Component Analysis)，主坐标分析(PCoA，Principal Co-ordinates Analysis)，非加权组平均聚类分析(UPGMA，Unweighted Pair-group Method withArithmetic Means)等分析方法，从中发现不同样品（组）间的差异。

# UniFrac距离
# 由于微生物极其多样，不同微生物彼此之间的系统发育关系往往千差万别，仅仅将群落中不同微生物成员视为相互独立的变量显然并不合理。
# 因此，在比较不同群落样品之间的差异时，需要考虑两个群落成员之间的系统发育关系是否相似。
# 基于这个思想，计算微生物群落样品间距离的UniFrac距离应运而生，通过比较两个群落各自独有的微生物成员之间系统发育关系的远近，更为客观地反映两个群落样品之间的相似程度。

# UniFrac距离有：
# 1）非加权（Unweighted）：仅仅考虑微生物成员在群落中存在与否，而不考虑其丰度高低。
# 2）加权（Weighted）：兼顾群落成员之间的系统发育关系以及它们在各自群落中的丰度高低。
# 两种距离算法侧重于不同的群落结构特征：究竟是由于群落成员的截然不同导致样品的差异，还是由于同一组成员在不同样品中丰度梯度的改变导致样品的差异。由于主坐标分析是以“无监督”的方式降维分解样品距离矩阵，因此，合理运用非加权和加权两种UniFrac距离，可以较全面地揭示微生物群落数据背后隐含的生态学意义（即UniFrac PCoA分析）。

# 如何计算UniFrac距离：
# https://www.jianshu.com/p/c8f3a4467dae

# 四种方法
picante::unifrac()
PhyloMeasures::unifrac.query()
GUniFrac::GUniFrac()
phyloseq::UniFrac()

# comm: 每一行是一个样品，每一列是一个OTU
# tree: 有根的系统发育树
# picante package
picante::unifrac(comm, tree)
# PhyloMeasures package, 下面只是最简单的用法，还有其它参数可以调节
PhyloMeasures::unifrac.query(tree, comm)
# GUniFrac package，下面只是最简单的用法，还有其它参数可以调节
GUniFrac::GUniFrac(comm, tree)
# phyloseq package，对数据类型有特殊要求，必须是phyloseq类
# 虽然phyloseq的安装稍微有些麻烦，并且在计算unfirac之前还需要先转换一下数据类型，但其计算unfirac的效率最高。
# 需要特别注意的是，使用phyloseq计算unfirac时，tree必须是二分叉的树，
# 意思就是每一个内部的node必须有且仅有两个子node。
# 判断一个树是不是二分叉的
ape::is.binary(tree)
# 如果不是，则使用下面的方法进行转变
tree <- ape::multi2di(tree)
# 然后基于新树计算unifrac
df <- phyloseq(otu_table(comm, taxa_are_rows = F), phy_tre(tree))
phyloseq::UniFrac(df, weighted=T)

# alpha diversity - Calculate Faith’s Phylogenetic Diversity
# 需要进化树，导入phylo object。
library(phyloseq)
tr <- ape::read.tree("../08.tree/faa.tre")
library(picante)
pd <- pd(otu_table, tr, include.root = F)
