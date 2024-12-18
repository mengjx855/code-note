#### info ####
# Jinxin Meng, 20231225, 20231228
pacman::p_load(tidyr, dplyr, tibble, purrr)
setwd("F:/Microbiome/proj_2023/02.Broiler_fecal_resistome_20230302/SEM/")

#### code ####
library(lavaan)
source("F:/Code/R_func/plot_PCoA.R")
source("F:/Code/R_func/plot_PCA.R")
source("F:/Code/R_func/diversity.R")
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
profile_args <- read.delim("../ARGs/ARGs_tpm", row.names = 1)
profile_mges <- read.delim("../MGEs/MGEs_tpm", row.names = 1)
profile_ko <- read.delim("../kegg/KO_tpm", row.names = 1)
profile_taxa <- read.delim("../taxa/species_tpm", row.names = 1)

dat <- list(
  args = calu_PCA(profile_args, cumulative_eig = 80, prefix = "args", add_eig = F),
  mges = calu_PCA(profile_mges, cumulative_eig = 80, prefix = "mges", add_eig = F),
  ko = calu_PCA(profile_ko, cumulative_eig = 80, prefix = "ko", add_eig = F),
  taxa = calu_PCA(profile_taxa, cumulative_eig = 80, prefix = "taxa", add_eig = F)
)
dat <- purrr::reduce(dat, \(x, y) merge(x, y, by = "sample"))

dat <- list(
  args = calu_alpha(profile_args, method = "shannon", out_colnames = "args"),
  mges = calu_alpha(profile_taxa, method = "shannon", out_colnames = "taxa"),
  ko = calu_alpha(profile_mges, method = "shannon", out_colnames = "mges"),
  taxa = calu_alpha(profile_ko, method = "shannon", out_colnames = "ko")
)
dat <- purrr::reduce(dat, \(x, y) merge(x, y, by = "sample"))


# 基本语法简介
# 语法一：f3~f1+f2（路径模型）
# 结构方程模型的路径部分可以看作是一个回归方程。而在R中，回归方程可以表示为y~ax1+bx2+c，“~”的左边的是因变量，右边是自变量，“+”把多个自变量组合在一起。那么把y看作是内生潜变量，把x看作是外生潜变量，略去截距，就构成了结构方程模型的语法一。
# 语法二：f1 =~ item1 + item2 + item3（测量模型）
# "=~"的左边是潜变量，右边是观测变量，整句理解为潜变量f1由观测变量item1、item2和item3表现。
# 语法三：item1 ~~ item1 , item1 ~~ item2
# "~~"的两边相同，表示该变量的方差，不同的话表示两者的协方差
# 语法四：f1 ~ 1
# 表示截距
# 此外还有其它高阶的语法，详见lavaan的help文档，一般的结构方程建模分析用不到，便不再列出。
# "=~" 利用被测变量(右)定义潜在变量(左)：测量模型。意思就是这几个指标构成了它们
# "~~"的两边相同，表示该变量的方差，不同的话表示两者的协方差，构建回归方程：路径模型。各测量模型之间的关系，本文中它们是两两相互影响的
model = '
  microbiome =~ taxa_1
  functions =~ ko_1
  resistome =~ args_1
  mobilome =~ mges_1
  
'
# functions ~ microbiome
# resistome ~ microbiome
# mobilome ~ microbiome
# resistome ~ functions
# functions ~ mobilome
# resistome ~ mobilome
# 适配模型
fit <- sem(model, data = dat)

# 查看拟合结果
# 潜在变量的路径系数大小: latent variables
# 协方差的路径系数大小：covariances
# 方差的路径系数大小: variances
summary(fit)

# 计算拟合系数
# 卡方值 | CHISQ | 矩阵整体相似程度 | P>0.05
# 拟合指数 | GFl | 说明模型解释力 |>0.90
# 相对拟合指数 | CFl | 说明模型较虚无模型的改善程度 | >0.95
# 未标准化残差 | RMR | 未标准化假设模型整体残差 | 越小越好
# 标准化残差 | SRMR | 标准化模型整体残差 | <0.08
# 近似均方根误差 |RMSEA | 理论模型与饱和模型的差异
fitMeasures(fit, c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea"))


#调用SEM绘图包
# library(semPlot)
semPaths(fit, what = "std", layout = "tree2", fade=F, nCharNodes = 0, intercepts = F, residuals = F, thresholds = F, )
# semPaths(fit, layout = "tree2", intercepts = F, residuals = F, thresholds = F)
# semPaths(fit, what = "std", layout = "tree2", fade=F, nCharNodes = 0)


# 令很多同学头疼的问题之一就是一个SEM中各种名称，这儿统一给大家总结一下：
# 自变量（Independent），预测因子（predictor），外生变量（exogenous (external)）都是一个东西，在模型中都是去影响别的变量的。
# 因变量（Dependent），标准变量（criterion），内生变量（endogenous (internal)）都是一个东西，表示别的变量的效果，在模型中受别的变量影响。
# 潜变量（Latent variable），因子（factor），构象（construct）都是一个东西，指你模型中的变量。
# 模型model指的就是你研究的变量之间的关系的统计表达。
# 如果你想把你的模型直观地画出来，画出来的这个东西就叫路径图path diagram
# 模型设定Specification就是指你对模型参数和整个模型的规定，这儿需要注意：
# 所有的模型都是错的，没有百分百吻合数据的模型，我们做SEM就是想要找到一个尽可能符合我们数据同时还符合理论解释的通的模型。
# 我们对于模型参数的设定要么是固定的，要么是设定为自由的
# 固定参数Fixed parameters就是指这个参数不从我们的数据中来估计，而是将它固定为0或者1。
# 自由参数Free parameters是指需要模型从数据中估计的参数。
# 拟合指数Fit indices反映模型拟合的如何，指的是我们设定固定和自由参数是不是和原始数据的方差协方差一致，chi-square, CFI, NNFI, RMSEA都是常见的拟合指数。
# 一个SEM可以划分为两个部分，一个部分叫做测量模型measurement model，指的是对潜变量和显变量关系的设定，另一个叫做结构模型structural model，指的是潜变量间或者潜变量与其余变量关系的设定。
# 识别Identification模型识别的意思是对模型中的自由参数能否获得一个特定的解，模型有解是需要满足一定条件的：就是就是自由参数的个数q必须要小于或等于你的样本协方差矩阵中非冗余元素的个数p*，这个P*=p(p + 1)/2，其中p为测量变量的个数。
# 识别的类型
# 不识别underidentified：就是从我们的数据中给一个或者多个参数找不出来一个特定的解。比如x + y = 5，这个xy可以来回变我们是找不到这个方程的特定解的，这个适合就叫做不识别，就是自由参数的个数q大于了独立方程的个数。
# 恰好识别just identified：就是对模型中的自由参数恰好可以从我们的数据中求得一个独立的解，比如我们有方程组x + y = 5 and 2x + y = 8，这个时候xy刚刚好有一个解，但是恰好识别的情形下，模型是无法被检验的。
