# Encoding: utf-8
# Creator：Jinxin Meng
# Created Date: 2023-12-14
# Modified Date: 2023-12-14

pacman::p_load(tidyr, dplyr, tibble, purrr)
pacman::p_load(ggplot2, ggpubr)
pacman::p_load(clusterProfiler, org.Ss.eg.db, enrichplot)
db=org.Ss.eg.db

# AnnotationDb对象和如何去访问其信息
# AnnotationDb对象包括OrgDb和ChipDb，OrgDb就是常见的org.Ss.eg.db（pig），org.Hs.eg.db（human）等
keytypes(db); columns(db) # 查看OrgDb中存储了基因的哪类型信息，不同物种的，不同版本的信息完整度不同。人的是最全的。
AnnotationDbi::select(db, keys = gene_list, keytype = "SYMBOL", columns = c("GENENAME", "ENTREZID")) # 使用gene_symbol去对应gene_name和entrezid


# 读入差异分析的结果
diff <- read.delim("diff_Treat_vs_ETEC.txt")  # gene|Treat_1|Treat_2|Treat_3|ETEC_1|ETEC_2|ETEC_3|baseMean|log2FoldChange|lfcSE|stat|pvalue|padj|enriched|gene_name
# 拿到差异基因向量，后边要做基因名称的转换，富集分析需要gene_id，这里用的是gene_symbol去对应gene_id，也可以用ensemb_id去对应gene_id，或者其他。
gene_list <- filter(diff, enriched != "None") %>% dplyr::select(gene_name) %>% unlist() %>% unique()
# 基因格式转换，使用AnnotationDBi:select函数进行转换
gene_dat <- AnnotationDbi::select(db, keys = gene_list,  keytype = "SYMBOL", columns = c("GENENAME", "ENTREZID"))

# GO 富集分析
# ont (ontology)可选BP，MF，CC，ALL等四个内容。
ego <- enrichGO(gene = gene_dat$ENTREZID, OrgDb = org.Ss.eg.db, ont = "ALL", qvalueCutoff = 0.2, pAdjustMethod = "BH")
# 输出注释结果，默认输出结果只有geneID这一列，我们从基因注释结果中，新加一列基因名，方便后续查找基因。这个表就是富集分析的结果了。
ego_res <- data.frame(ego) %>%  
  mutate(geneName = lapply((strsplit(x = geneID, "/")), 
                           function(x) gene_dat$SYMBOL[match(x, gene_dat$ENTREZID)]) %>% 
           sapply(., function(x) paste(x, collapse = "/")), .after = geneID)

# GO 富集分析 针对top30的BP条目画气泡图
ego_res %>% filter(ONTOLOGY == "BP") %>% head(n = 30) %>% mutate(GeneRatio2 = sapply(strsplit(GeneRatio, split = "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))) %>% 
  ggplot(., aes(x = GeneRatio2, y = reorder(ID, Count))) +geom_point(aes(fill = p.adjust, size = Count), shape = 21, stroke = .4) +
  scale_fill_gradient(high = "#fcbba1", low = "#cb181d", trans = "log10") + scale_size_continuous(range = c(2, 5)) +
  labs(x = "GeneRatio", y = "", title = paste0("Biological Process top 30 ", paste0(gp[[i]][-1], collapse = " vs."))) +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"), axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"), axis.line = element_blank(),
        plot.title = element_text(size = 10, color = "black"), panel.grid = element_line(linewidth = .4, color = "grey90"),
        panel.border = element_rect(linewidth = .4, color = "black", fill = NA), legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"), legend.background = element_blank(),
        legend.key.height = unit(4, "mm"), legend.key.width = unit(4, "mm"), legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(2, "mm"),
        aspect.ratio = 2) +
  guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = .3, ticks.colour = "black", ticks.linewidth = .3, order = 1))
# Y叔的包也可以快速出图，基于整体的画图，过滤不方便了。
dotplot(ego, showCategory = 10, split = "ONTOLOGY")

# GO 富集分析 针对各类top10条目画柱状图
ego_res %>% group_by(ONTOLOGY) %>% top_n(., n = 10) %>% mutate(GeneRatio2 = sapply(strsplit(GeneRatio, split = "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))) %>% 
  ggbarplot(., x = "Description", y = "GeneRatio2", fill = "ONTOLOGY", sort.val = "asc", sort.by.groups = T, orientation = "horizontal",
            ylab = "Gene Ratio", xlab = "", title = paste0(paste0(gp[[i]][-1], collapse = "_vs_"), " GO enrichment analysis")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"), axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"), plot.title = element_text(size = 10, color = "black"),
        panel.grid = element_line(linewidth = .4, color = "grey90"), legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"), legend.position = "right", aspect.ratio = 2)
# Y叔的包也可以快速出图，基于整体的画图
barplot(ego, showCategory = 20)

# GO term与差异基因关系网络图，叫做Gene-Concept Network，就是差异基因和term的网络图，属于某一个term的基因和这个term的节点连起来
cnetplot(ego, showCategory = 3, layout = "kk", circular = F) +  theme(aspect.ratio = 1) # 这个我感觉没啥必要做
# 热图样功能分类
heatplot(ego, showCategory = 5) # 这个我感觉没啥必要做
# 树状图
ego_x <- pairwise_termsim(ego) # 先计算一个相似性结果，默认是jaccard相似度
treeplot(ego_x, showCategory = 10)
# 富集图将富集的term组织成一个网络，其边缘连接重叠的基因集。这样，相互重叠的基因集往往会聚集在一起，从而很容易识别功能模块。
ego_x <- pairwise_termsim(ego)
emapplot(ego_x, cex_category = .5,layout = "kk", cex.params = list(category_label = .4))

# KEGG 富集分析
# 每次都会从KEGG官网下载最新的数据
ekegg <- enrichKEGG(gene = na.omit(gene_dat$ENTREZID), organism = "ssc", qvalueCutoff = 0.2, pAdjustMethod = "BH")

# kegg 通路图，要准备一个这样的文件，行名为基因名，第一列是log2FC
#           log2FoldChange
# 100524987 0.2969393
# 100525064 -0.2808995
library(pathview)
pathview(gene.data = dat, pathway.id = "ssc04010", species = "ssc", out.suffix = "fmt1") 
# kegg.native=F可以做出pdf版本的通路图，sign.pos可以画出通路符号的含义，
pathview(gene.data = dat, 
  pathway.id = "ssc04010", # 通路id
  species = "ssc", # 物种
  kegg.native = F, # 是否采用graphviz视图
  sign.pos= "bottomleft", 
  out.suffix = "fmt2",
  low = list("#74add1", "#74add1"), mid = list("white", "white"), high = list("#f46d43","#f46d43"), # 设置中高低的颜色
  limit = list(1, 1), # 图例范围
  same.layer = T, # 颜色和底层的图是否在一个图层
  pdf.size = c(10, 10))

# emapplot也是一种网络图，不过可以把相似的条目聚集到一起，便于识别不同的功能模块。
# 这个函数同样也是可以直接使用enrichResult、gseaResult、compareClusterResult3种结果。
# 不过在使用前，必须用pairwise_termsim函数添加相似性矩阵才行~
ekegg_pt <- pairwise_termsim(ekegg)
emapplot(ekegg_pt, showCategory = 30 ,
  color = "p.adjust", # 映射给条目颜色 ’pvalue’, ’p.adjust’ or ’qvalue’
  shadowtext = T, # 显示标签阴影
  repel = FALSE, # 解决标签重叠问题
  node_label = "category", #显示谁的标签  ’category’,’group’,’all', ’none’
  #形状控制参数
  layout.params = list(layout = NULL, #和cnetplot的形状参数一样，后面也是，就不多说了
    coords = NULL), #控制位置，需要含2列的data.frame,x是x轴坐标,y是y轴坐标
  #控制连线
  edge.params = list(show = TRUE, # 是否显示连线
    min = 0.2), #判断两个节点是否相似的阈值
  #大小控制参数
  cex.params = list(category_node = 1, # 节点大小
    category_label = 1, #节点标签大小
    line = 1, # 线的粗细
    pie2axis, #饼图大小
    label_group), #分组标签大小
  #控制高亮，和cnetplot同，不再多说
  hilight.params = list(category = NULL,
    alpha_hilight = 1,
    alpha_no_hilight = 0.3),
  #控制聚类
  cluster.params = list(cluster = FALSE, #对条目聚类
    method = stats::kmeans, #聚类方法
    n = NULL, #聚类个数
    legend = FALSE, #显示聚类图例
    label_style = "shadowtext",
    label_words_n = 4, #聚类标签个数
    label_format = 30)#最长的字符数，控制折叠
  )
