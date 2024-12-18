####info####
# Encoding: utf-8
# modifier：Jinxin Meng
# Modified date：2022-9-25

####option####
options(stringsAsFactors = F, scipen = 5)
pacman::p_load(dplyr, tidyr, tibble, stringr, ggplot2, gggenes)
setwd("F:/Proj/25.Virome_IBS_stat_2022_8_29/04./")

####data####
KO <- read.delim("KO_list2", sep = "\t", header = F) %>% 
  rename(gene = V1, KO = V2)

gene_type <- read.delim("hmmsearch_gene_info", sep = "\t", header = F) %>% 
  rename(gene = V1, anno = V2, type = V3) %>% 
  select(gene, type) %>% 
  group_by(gene, type) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  slice_max(count, n = 1, with_ties = F)

KO_type <- read.delim("KO_type.txt", sep = "\t", header = T) 

genome <- read.delim("genome_gene_local.txt", sep = "\t", header = T) %>% 
  mutate(label = gene_type$type[match(gene, gene_type$gene)],
         label = ifelse(is.na(label), "Unknown", label),
         label2 = KO$KO[match(gene, KO$gene)],
         label3 =  KO_type$type[match(label2, KO_type$KO)],
         label3 = ifelse(is.na(label3), label, label3),
         label2 = ifelse(is.na(label2), "", label2))

genome2 <- genome %>% 
  subset(!genome%in%"ZhuQ_2021|SRR11992801_k141_169070")

var <- c("viral", "microbial","Unknown","Butanoate_metabolism","CAZy_CBM50","CAZy_GH24",
         "Fatty_acid_biosynthesis","Folate_biosynthesis","Lipopolysaccharide_biosynthesis","Sulfur_metabolism")

color <- c("#d9d9d9", "transparent", "transparent", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3","#fdb462", "#b3de69")

ggplot(genome2, aes(xmin = start, xmax = end, y = genome, forward = direction, fill = factor(label3, levels = var))) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), color = "grey40") +
  geom_text(aes(x = start + (end - start)/2, label = label2, color = label2), angle = 90, size = 2) +
  scale_fill_manual(values = color) +
  scale_x_continuous(expand = c(.01, .01)) +
  labs(x = "", y = "", fill = "Gene type") +
  theme_genes()


genome3 <- genome %>% 
  subset(genome%in%"ZhuQ_2021|SRR11992801_k141_169070")

ggplot(genome3, aes(xmin = start, xmax = end, y = genome, forward = direction, fill = factor(label3, levels = var))) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), color = "grey40") +
  geom_text(aes(x = start + (end - start)/2, label = label2, color = label2), angle = 90, size = 2) +
  scale_fill_manual(values = color) +
  scale_x_continuous(expand = c(.01, .01)) +
  labs(x = "", y = "", fill = "Gene type") +
  theme_genes()
