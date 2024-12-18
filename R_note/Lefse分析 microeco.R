# lefse --------
# microeco是基于R6class开发的
library(ggtree)
library(microeco)
library(magrittr)

profile <- read.delim("table_rarefied.tsv", row.names = 1, check.names = F) %>%
  profile_transRA()
#       s1  s2  s3
#OTU1   1   2   3
#OTU2   4   9   0   
#OTU3   3   3   3
group <- read.delim("sample_group", row.names = 1) %>% 
  filter(rownames(.) %in% colnames(profile))
#   group
#s1 g1
#s2 g2
#s3 g2
tax <- read.delim("taxonomy.tsv", row.names = 1)
#       domain  phylum  class   order   family  genus   species
#OTU1   d1      p1      c1      o1      f1      g1      s1
#OTU2   d2      p2      c2      o2      f2      g2      s2 
#OTU3   d1      p1      c3      o3      f3      g3      s3
tr <- read.tree("tree.nwk")

set.seed(2024)
#清除物种表中哪些unknown, unclassified, uncultured .... 
# %<>%的用法就等于 tax = tidy_taxonomy(tax)
tax %<>% tidy_taxonomy

#构建microeco dataset
dataset <- microtable$new(sample_table = group, otu_table = profile, tax_table = tax, phylo_tree = tr)
#lefse分析
lefse <- trans_diff$new(dataset = dataset, method = "lefse", group = "group", alpha = 0.05, lefse_subgroup = NULL)
#lefse柱状图
lefse$plot_diff_bar(use_number = 1:30, width = 0.6, group_order = c("Control", "IBD"), color_values = c("#69a7bc","#E69F00"))
ggsave("lefse.bar.pdf", width = 8, height = 10)
#lefse进化树
#此数据集中的分类单元太多。例如，我们只选择树中丰度最高的200个分类群和50个差异特征。
lefse$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, clade_label_level = 5, 
                          group_order = c("Control","IBD"), color = c("#69a7bc","#E69F00"))
ggsave("lefse.cladogram.pdf", width = 20, height = 8)
#可以指定标记的类群名称
use_labels <- c("p__Proteobacteria", "p__Bacteroidetes", "p__Acidobacteria", "p__Firmicutes")
lefse$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, select_show_labels = use_labels)
#可以对检测到的差异菌进行丰度可视化
# 然后，可以很容易地将LEfSe检测到的生物标志物的丰度可视化。
lefse$plot_diff_abund(use_number = 1:20, group_order = c("Control", "IBD"), color_values = c("#69a7bc","#E69F00"))
#lefse结果表格，并保存
lefse$res_diff %>% head

