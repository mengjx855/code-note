####info####
# created date: 2022-9-22
# modified date: 2022-10-21

####option####
options(stringsAsFactors = F)
pacman::p_load(dplyr, stringr, tidyr, tibble, ggplot2, ggrepel)
setwd("/share/data1/mjx/proj/12.MDD_Virome/07.network/")

####code####
bac_info <- read.delim("../06.metaphlan/bacteria_taxa.txt", sep = "\t")
vir_info <- read.delim("../03.comparison/marker_vOTUs_taxa_host.txt", sep = "\t")
vir_name <- read.delim("../03.comparison/marker_trans_names.txt", sep = "\t")

edges <- rbind(corr %>% 
                 select(virus, bacteria, corr) %>% 
                 add_column(cohort = "disease") %>% 
                 rename(from = virus, to = bacteria) %>% 
                 mutate(benriched = bac_info$enriched[match(to, bac_info$bacteria)]),
               corr2 %>% 
                 select(virus, bacteria, corr) %>% 
                 add_column(cohort = "control") %>% 
                 rename(from = virus, to = bacteria) %>% 
               mutate(benriched = bac_info$enriched[match(to, bac_info$bacteria)]))

nodes <- rbind(data.frame(name = unique(edges$from)) %>% 
                 mutate(family = vir_info$family[match(name, vir_info$vOTUs)],
                        enriched = vir_info$enriched[match(name, vir_info$vOTUs)],
                        name = vir_name$marker2[match(name, vir_name$marker)]) %>% 
                 add_column(type = "virus"),
               data.frame(name = unique(edges$to)) %>% 
                 mutate(family = bac_info$family[match(name, bac_info$bacteria)],
                        enriched = bac_info$enriched[match(name, bac_info$bacteria)]) %>% 
                 add_column(type = "bacteria"))
edges <- mutate(edges, from = vir_name$marker2[match(from, vir_name$marker)])


####ggraph####
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

edges2 <- subset(edges, cohort%in%"disease")

graph <- tbl_graph(nodes = nodes, edges = edges2)

color <- c("#8dd3c7","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf",
           "#6a3d9a","#666666","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5")

ggraph(graph, layout = "linear", circular = T) + 
  geom_edge_arc(aes(color = factor(benriched, levels = c("Disease", "Control"))), width = .4) +
  geom_node_point(aes(shape = type, fill = family), color = "black", size = 1.7, show.legend = F) +
  geom_node_text(aes(x = x*1.04, y = y*1.04, label = name, angle = atan(y/x)*360/(2*pi), hjust = "outward",
                     color = factor(enriched, levels = c("Disease", "Control"))), size = 2) +
  scale_color_manual(values = c("#f77d4d", "#1f78b4")) +
  scale_edge_color_manual(values = c("#f77d4d", "#1f78b4")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = color) +
  labs(color = "enriched", edge_color = "enriched") +
  coord_fixed() +
  theme_graph() +
  theme(aspect.ratio = 1)

ggplot(unique(select(nodes, family)), aes(x = 1, y = family , fill = family)) +
  geom_point(shape = 21, size = 8, show.legend = F) +
  geom_text(aes(x = 1.1, label = family), hjust = 0) +
  scale_fill_manual(values = color) +
  xlim(c(1,4)) +
  theme_void()
ggsave("legend.pdf", width = 4.5, height = 6)