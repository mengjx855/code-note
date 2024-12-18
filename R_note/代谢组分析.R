# Jinxin Meng, 20240801, 20240811 ---------------------

setwd("F:/proj/proj_2024/20240731_BC_IBD_Zhangyn/analysis_metabolites/")

pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)
library(ropls)
library(clusterProfiler)

source("/code/R_func/difference_analysis.R")

# difference and enrichment analysis ----------

metadata <- read.delim("metadata.txt")
profile <- read.delim("profile.txt", row.names = 1)
group <- read.delim("group.txt")

gps <- list(
  c("bc_0d", "pbs_0d"),
  c("bc_3d", "pbs_3d"))

gp <- gps[[2]]

samples <- group %>% 
  filter(group %in% gp) %>% 
  pull(sample)
profile_x <- profile %>% 
  select(all_of(samples))
group_x <- group %>% 
  filter(sample %in% samples)

# 监督分组以疾病-健康为例，orthoI = NA 时执行 OPLS-DA, 自动计算正交分量的数量；默认设置为0，执行的是PLS分析
oplsda <- opls(x = data.frame(t(profile_x), check.names = F), pull(group_x, group), orthoI = NA, predI = 1)

diff <- difference_analysis(profile_x, group_x, gp = gp, diff_method = "wilcox")

out <- oplsda@vipVn %>%
  data.frame(vip = .) %>% 
  rownames_to_column(var = "cpd_id") %>% 
  left_join(diff, ., by = c("name" = "cpd_id")) %>% 
  mutate(enriched = ifelse(vip > 1 & log2FC > 1, gp[1], 
                           ifelse(vip > 1 & log2FC < -1, gp[2], "none")))
write.table(out, paste0("diff.", gp[1], "_vs_", gp[2], ".tsv"), sep = "\t", quote = F, row.names = F)
table(out$enriched)

# 富集分析
bg_data <- read.delim("/database/KEGG/v20230401/cpd2path_enrichment.tsv")
cpds <- out %>% 
  filter(enriched == gp[1]) %>% 
  select(cpd_id = name) %>% 
  left_join(metadata, by = "cpd_id") %>% 
  filter(cpd_KEGG != "") %>% 
  pull(cpd_KEGG)

eKEGG <- enricher(gene = cpds, TERM2GENE = bg_data, minGSSize = 1, pvalueCutoff = 1, qvalueCutoff = 1)
eKEGG_out <- data.frame(eKEGG@result)
write.table(eKEGG_out, paste0("eKEGG.", gp[1], "_vs_", gp[2], ".up.tsv"), sep = "\t", quote = F, row.names = F)

plot_data <- out %>% 
  select(name, log2FC, vip, enriched) %>% 
  mutate(enriched = factor(enriched, c("bc_3d", "none", "pbs_3d")),
         log2FC = ifelse(log2FC > 3, 3, log2FC),
         log2FC = ifelse(log2FC < -3, -3, log2FC)) %>% 
  left_join(metadata %>% select(1:2, cpd_class), by = c("name" = "cpd_id")) %>% 
  mutate(label = ifelse(enriched != "none", cpd_name, ""))
  # mutate(label = ifelse((grepl("Fatty acid", cpd_class) & enriched != "none") | cpd_name == "3-Hydroxybutyric acid", cpd_name, ""))

# 火山图
ggscatter(plot_data, "log2FC", "vip", color = "enriched", legend = "right",
          palette = c("#eeacec", "grey", "#21c1dc"), size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text(aes(label = label), size = 1, fontface = "italic")
# ggsave("diff.bc_3d_vs_pbs_3d.vol.v2.pdf", width = 5, height = 3.5)