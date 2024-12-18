# Encoding: utf-8
# Creator：Jinxin Meng
# Created Date: 2023-11-15
# Modified Date: 2023-11-15

# https://www.jianshu.com/p/22fdbbc5b5c7
#
#
#
#


library(ggtree)
library(tidytree)

# 怎么查找父子节点呢？
# 先得到treedata，在treedata中找每个tip的节点编号
tr_dat <- tr %>% as.treedata() %>% as_tibble()
parent(); child(); offspring(); ancestor(); sibling(); MRCA()
# 这些函数参数需要输入一个树数据，一个node编号
# 查到node编号，可以进行对一个tip或者一个clade注释，例如
geom_cladelab(node = 125, label = "label", offset = .02, barsize = 1, barcolour = "black", hjust = -.1)
# 查询label的节点编号
nodeid(tr_dat, label)

# 对tip进行修饰和注释
geom_tiplab() # 注释节点的标签，必须是映射，不能单独指定
geom_tippoint() # 注释节点的形状，必须是映射，不能单独指定
# 想要对某个tip单独进行表示，可以用这个命令
geom_point2(aes(subset=(node==20))) # node指定节点的编号

# 对bratch背景进行修饰
geom_hilight() # 对分支背景上色
df = data.frame(node = c(1, 2), type = c("A1", "A2"))
ggtree(tr) + geom_hilight(data = df, aes(node = node, type = type), align = "both") # align:[both|right|left] 控制背景填充的位置


# re-root树
tr <- root(tr, outgroup = 109)
