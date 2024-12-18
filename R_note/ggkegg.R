#### Jinxin Meng, 20241109, 20241109 ####

pacman::p_load(dplyr, tidyr, tibble, purrr, stringr, ggplot2)
pacman::p_load(ggkegg, tidygraph)

# https://noriakis.github.io/software/ggkegg/index.html
# ggkegg fetches information from KEGG and parse, analyze and visualize them using ggplot2 and ggraph, 
# combined with the other packages investigating biological functions using KEGG. 
# This package aims to visualize the complex components of KEGG using the grammar of graphics. 
# For Python, please use pykegg using plotnine, which offers the almost same functionality as ggkegg 
# used in conjuction with the package such as gseapy and PyDESeq2 and single-cell transcriptomics analysis library scanpy.

#### code ####
# One of the main aims of ggkegg is manupilating KEGG information in tidy ways using tidygraph, 
# and offers the customized visualization of KEGG information including KEGG PATHWAY, MODULE, and NETWORK.

pathway("hsa04110") %>%  ## Obtain and parse the KEGG pathway
  activate(nodes) %>%  ## node manipulation
  mutate(convert_hsa=convert_id("hsa"),
         convert_map=convert_id("pathway")) %>% ## convert IDs for organism hsa and pathway
  ggraph(x=x, y=y)+ ## ggraph plot
  geom_edge_parallel(arrow=arrow(length=unit(1,"mm")),
                     aes(linetype=subtype_name),
                     end_cap=circle(7.5,"mm"))+ ## Parallel edges
  geom_node_rect(aes(filter=type=="gene",
                     fill=I(bgcolor)),
                 color="black")+ ## rectangular nodes
  geom_node_text(aes(label=convert_hsa),
                 size=2, family="serif")+ ## text
  geom_node_text(aes(label=convert_hsa,
                     filter=!is.na(convert_hsa) & convert_hsa=="TP53"),
                 size=2, color="red", family="serif")+ ## highlight text
  theme_void()

highlight_entities("hsa04110", "CDKN2A")

# 获取数据信息，返回是一个tbl_graph对象，需要使用ggraph可视化
pathway('hsa04110')
module('M00002')
network('N01592') %>% network_graph() %>% plot_kegg_network() # 可以使用对应的函数构建igraph对象

# 从这些信息中获取nodes或者edges表
pathway('hsa04110') %>% activate(nodes) %>% as_tibble()
