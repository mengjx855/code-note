# Jinxin Meng ####

# 针对个别的基因展示其丰度
FeaturePlot(seurat, reduction = "tsne", features = c("Lgr5", "Hmgcs2")) + coord_fixed()
FeaturePlot(seurat, features = "Hopx", split.by = "orig.ident") # 也可以按照分组分开几个图


# 提取每个基因的表达数据
data <- FetchData(seurat, vars = c("Lgr5", "Hmgcs2"))
data <- FetchData(seurat, vars = c("Lgr5", "Hmgcs2")) %>% 
  rownames_to_column("name") %>% 
  left_join(rownames_to_column(data.frame(seurat@reductions$tsne@cell.embeddings), "name"), by = "name")


# 计算并绘制模块得分
seurat <- AddModuleScore(seurat, list(c("Lgr5", "Hmgcs2")), name = "combine")
FeaturePlot(seurat, reduction = "tsne", features = "combine1", pt.size = .6)

# 使用 blend 参数（仅支持两个基因），如果你只需要展示两个基因的共表达，可以使用 blend = TRUE 的方法：
# 在UMAP或t-SNE图中展示两个基因的共表达
FeaturePlot(seuratect, features = c("GeneA", "GeneB"), blend = TRUE)

# 展示不同样本来源的细胞，着不同颜色
# 如果orig.ident中没有对应的样本信息的话，需要自己补充
seurat$source <- c(rep("sample1", time = 1000), rep("sample2", time = 800))
DimPlot(seurat, reduction = "umap", group.by = "orig.ident", cols = c("#3e4faa", "#ea3726")) # 展示在一张图中，着色不同
DimPlot(seurat, reduction = "umap", split.by = "orig.ident") # 分面展示


# 比较基因表达情况
VlnPlot(seurat, features = "Hopx", group.by = "orig.ident")


# 寻找细胞间的差异基因作为分群的marker
# 找到所有差异表达基因
allmarkers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 提取每个簇中具有最显著差异表达的前5个基因
top5_markers <- allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
# 创建基因的Dot Plot
DotPlot(seurat, features = unique(top5_markers$gene)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
  NoLegend()
# 创建基因的Heatmap
DoHeatmap(seurat, features = unique(top5_markers$gene)) +
  theme(axis.text.y = ggplot2::element_text(size = 8)) + 
  NoLegend()
# 比较指定的组
# 找到 cluster 1 相对于其他所有类群的差异基因
differential_genes <- FindMarkers(seurat, ident.1 = 1)
# 如果你想比较两个特定类群，比如 cluster 1 和 cluster 2
differential_genes <- FindMarkers(seurat, ident.1 = 1, ident.2 = 2)
FindMarkers(seurat_mono, ident.1 = x, ident.2 = "control", group.by = "group2")



# 添加 metadata
group <- read.delim("data/sample_group.txt")
seurat <- data.frame(sample = seurat@meta.data$orig.ident) %>% 
  left_join(group, by = "sample") %>%
  dplyr::select(-sample) %>% 
  AddMetaData(seurat, .)


# marker基因在各个类群中表达小提琴图
seurat_mono <- readRDS('mono_cell/seurat.mono.out.rds')
plot_data <- FetchData(seurat_mono, vars = c("CD14","FCGR3A")) %>% 
  add_column(cluster = paste0('clu_', seurat_mono@active.ident)) %>% 
  gather(key = "gene", value = "value", -cluster)
ggplot(plot_data, aes(value, cluster), color = factor(cluster)) +
  geom_violin(aes(fill = cluster), scale ="width") +
  facet_grid(cols = vars(gene), scales = "free") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 1, vjust = NULL, color = "black", size = 10),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(0, "cm"),
        strip.text.x = element_text(angle = 0, size = 10, hjust = .5, face = "italic"),
        strip.background.x = element_blank())



#### Jinxin Meng, 20240711, 20240805 全部分析内容 ####
pacman::p_load(dplyr, tidyr, tibble, purrr, stringr, ggplot2, ggpubr)
pacman::p_load(Seurat, harmony, data.table, cowplot)

# read raw data
path <- list.files("../data/", full.names = T)
name <- list.files('../data/')

seurat_list <- map2(path, name, \(x, y) {
  message(paste0("Processing: ", y)) 
  count <- Read10X(x)
  CreateSeuratObject(counts = count, project = y, min.cells = 5, min.features = 100)
})

# 查看各个样本中检测到基因和细胞的数量
do.call(rbind, lapply(seurat_list, dim))

# 合并多个对象后再合并layers
seurat <- merge(seurat_list[[1]], seurat_list[-1], add.cell.ids = name) %>% JoinLayers()
saveRDS(seurat, "seurat.raw.rds")
# seurat <- readRDS("seurat.rds")




#### Seurat QC质控 ####
# 计算线粒体基因比例，小鼠数据基因名为小写"^mt-"
# 线粒体基因表达比例较高的细胞通常被认为是质量较低或处于应激状态的细胞。
# 这些细胞可能正在经历凋亡或者其他形式的细胞死亡，因此通常会在预处理步骤中进行过滤。
mit_genes <- rownames(seurat)[grep("^mt-", rownames(seurat), ignore.case = T)]
seurat <- PercentageFeatureSet(seurat, features = mit_genes, col.name = "percent_mit") # 获取这些目的基因的所有counts占总count的比例
ggdensity(seurat@meta.data$percent_mit, xlab = "percent_mit")

# 计算核糖体基因比例，
# 核糖体基因（ribosomal genes）在细胞中的表达通常与细胞的生长和增殖状态有关。
# 高水平的核糖体基因表达可能表明细胞处于活跃的蛋白质合成状态，这通常与细胞增殖相关。
# 在单细胞 RNA 测序（scRNA-seq）数据分析中，核糖体基因的高表达有时被用来识别增殖活跃的细胞。
# 秦哲文：一般不去除核糖体的基因；
rib_genes <- rownames(seurat)[grep("^Rp[sl]", rownames(seurat), ignore.case = T)]
seurat <- PercentageFeatureSet(seurat, features = rib_genes, col.name = "percent_rib") # 获取这些目的基因的所有counts占总count的比例
ggdensity(seurat@meta.data$percent_rib, xlab = "percent_rib")

# 计算红血细胞基因比例
hb_genes <- rownames(seurat)[grep("^Hb[^(p)]", rownames(seurat), ignore.case = T)]
seurat <- PercentageFeatureSet(seurat,  features = hb_genes, col.name = "percent_hb") # 获取这些目的基因的所有counts占总count的比例
ggdensity(seurat@meta.data$percent_hb, xlab = "percent_hb")

# 可视化细胞的上述比例情况
p1 <- VlnPlot(seurat, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA"), 
              pt.size = 0, ncol = 2) + NoLegend()
p2 <- VlnPlot(seurat, group.by = "orig.ident", features = c("percent_mit", "percent_rib", "percent_hb"), 
              pt.size = 0, ncol = 3) + scale_y_continuous(breaks = seq(0, 100, 10)) + NoLegend()
cowplot::plot_grid(p1, p2, ncol = 1, align = 'v')
ggsave('qc/qc_before.pdf', width = 6, height = 10)

# 过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
# https://www.sciencedirect.com/science/article/pii/S2352345X19300797#sec3
# 这里，我们筛选的标准是：percent_mito < 25，percent_ribo > 3，percent_hb < 1。
# 不同的数据集，不同的组织，需要根据特定情况来调整筛选阈值。
cells <- reduce(list(WhichCells(seurat, expression = percent_mit < 25),
                     WhichCells(seurat, expression = nFeature_RNA > 100),
                     WhichCells(seurat, expression = nFeature_RNA < 10000),
                     WhichCells(seurat, expression = percent_rib > 3),
                     WhichCells(seurat, expression = percent_hb < 1)),
                \(x, y) base::intersect(x, y))
seurat <- subset(seurat, cells = cells)

# 过滤后查看结果
p1 <- VlnPlot(seurat, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA"), 
              pt.size = 0, ncol = 2) + NoLegend()
p2 <- VlnPlot(seurat, group.by = "orig.ident", features = c("percent_mit", "percent_rib", "percent_hb"), 
              pt.size = 0, ncol = 3) + NoLegend()
cowplot::plot_grid(p1, p2, ncol = 1, align = 'v')
ggsave('qc/qc_after.pdf', width = 6, height = 10)

# 查看质控后的单细胞矩阵信息
dim(seurat)
#     基因数  细胞数量
# [1] 17783   16707




#### 去批次，降维，聚类 ####
# 在拿到下游单细胞矩阵前，样本经历了多个实验环节的处理，时间、处理人员、试剂以及技术平台等变量都会成为混杂因素。
# 以上因素混合到一起，就会导致数据产生批次效应（batch effect）。
# 为了尽可能避免混杂因素，我们可以严格把控测序的技术流程，同时也需要在下游分析中进行事后补救（算法去批次）。
# 目前单细胞测序常用的去批次算法有Harmony，CCA，RPCA,FastMNN,scVI等。在这里，我们采用Harmony进行演示。

# 数据标准化
# 数据标准化的方法是通过对原始表达值进行对数转换"LogNormalize"，使其总体更加符合正态分布。
# 经过对数转换之后，不同基因或细胞之间也更具有可比性，从一定程度上消除不同细胞之间的技术差异。
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = T) 

# 筛选高变基因
# 高变异基因就是highly variable features（HVGs）是指在某些细胞中高度表达，而在其他细胞中低度表达的基因。
# 默认情况下，此步骤回筛选2000个HVGs用于下游分析。
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000) # 筛选高变基因
VariableFeaturePlot(seurat, cols = c("grey","red"))

# 数据归一化
# scale归一化：将每个基因在所有细胞中的均值变为0，方差标为1，
# 赋予每个基因在下游分析中同样的权重，不至于使高表达基因占据主导地位
seurat <- ScaleData(seurat)

# PCA线性降维分析
# 单细胞测序是一种高通量测序技术，它产生的数据集在细胞和基因数量上都具有很高的维度。
# 这立即指向了一个事实，即单细胞测序数据受到了"维度诅咒"的困扰（单细胞最好的教程（四）：降维）。
# "维度诅咒"这个概念最早是由R. Bellman提出的，它描述的是理论上高维数据包含更多的信息，但在实践中并非如此。
# 更高维度的数据往往包含更多的噪声和冗余，因此增加更多的信息并不利于后续的分析步骤。
# 为了在数据的处理和可视化更加便捷的同时，保留数据重要的信息，
#  在这里我们需要应用PCA（principal components analysis）即主成分分析技术，降低数据维度（单细胞PCA降维结果理解）。
seurat <- RunPCA(seurat, features = VariableFeatures(seurat))

# 流石图协助选择PC维度
# Elbow Plot用于选择合适数量的主成分时，通常寻找“肘部”位置，即图中解释方差增益开始显著减少的拐点。
# 这个拐点代表了从该主成分开始，增加更多的主成分对总体解释方差的贡献显著降低，因此可以作为选择主成分数目的参考点。
# 每新增一个PC轴，这些PC对整体的差异解释的越多，方差会越来越小。
ElbowPlot(seurat, ndims = 30)

# Harmony去批次
# Harmony应用主成分将转录组表达谱嵌入到低纬度空间，然后进行迭代过程去除特有数据集的影响。
# 原理可以查看原文，目的就是减少不同样本或不同批次的数据之间的差异。
seurat <- RunHarmony(seurat, "orig.ident")

# 然后使用UMAP/TSNE可视化
seurat <- RunUMAP(seurat, dims = 1:20, reduction = "harmony")
DimPlot(seurat, reduction = "umap", label = F)
seurat <- RunTSNE(seurat, dims = 1:20, reduction = "harmony")
DimPlot(seurat, reduction = "tsne", label = F)




#### 细胞聚类 ####
# Seurat软件使用基于图论的聚类算法对细胞进行聚类和分群，要包括以下步骤：
# 1.构建细胞间的聚类关系：利用PCA空间中的欧几里得距离构建KNN聚类关系图；
# 2.优化细胞间聚类关系距离的权重值：利用Jaccard相似性优化任意两个细胞间的边缘权重；
# 3.聚类和分群：使用Louvain 算法进行细胞群聚类优化。
# 分群标准的确定使用了两个函数FindNeighbors和FindClusters。
# FindNeighbors函数用于计算给定数据集的最近邻图，可以返回包含KNN和SNN的对象列表。
# 还可以通过迭代分组来优化聚类结果，以最大化标准模块度函数。
# FindClusters函数是基于共享最近邻（SNN）模块化优化的聚类算法识别细胞簇。
# 该函数的重要参数为分辨率resolution,resolution最小值为0，分为1类；值越大，分群数越多；
# 在0.4-1.2之间通常会对3k 左右的单细胞数据集产生良好的结果。

# 运行FindNeighbors:
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:12) 

# 然后再运行FindClusters，设置了不同的分辨率，观察分群效果
seurat <- FindClusters(seurat, resolution = c(seq(0, 1, .1))) # 探究0~1分辨率，间隔为0.1

map2(paste0("RNA_snn_res.", seq(0, 1, .1)), paste0("Resolution: ", seq(0, 1, .1)),
     \(x, y) DimPlot(seurat, reduction = "tsne", group.by = x, label = T) +
       ggtitle(y) + 
       theme(aspect.ratio = 1) + 
       guides(color = "none")) %>% 
  plot_grid(plotlist = ., nrow = 3)
ggsave("cluster/tsne.resolution.dimplot.jpg", width = 16, height = 13)




#### 细胞注释 ####
# 细胞基因标记
markers <- list(
  Stem_cells = c("Lgr5","Smoc2","Ascl2"),
  Paneth_cells = c("Ang4", "Lyz1"),
  Goblet_cells = c("Muc2", "Ccl9"),
  Enteroendocrine_cells = c("Chga", "Neurog3", "Reg4"),
  Tuft_cells = c("Dclk1", "Cd24a", "Trpm5"),
  Enterocyte = c("Aldob", "Apoa1", "Gsta1","Fabp1", "Prap1")
)

# 标志基因的featureplot
map2(markers, names(markers), \(x, y) {
  ggtext(data.frame(x = 1, y = 1), "x", "y", label = stringr::str_to_title(gsub("_", " ", y)), 
         size = 14, face = "bold") + 
    theme_void() + 
    FeaturePlot(seurat, reduction = "tsne", features = x, ncol = length(x), coord.fixed = T) +
    patchwork::plot_layout(widths = c(1, length(x)))
  ggsave(paste0('annotation/tsne.featureplot.', y,'.jpg'), width = 4 * length(x) + 4, height = 4)
  } )



FeaturePlot(seurat, reduction = "tsne", features = markers$ISCs)

# resolution .6

seurat <- FindClusters(seurat, resolution = .8) 
DimPlot(seurat, reduction = "tsne", label = TRUE, pt.size = 0.5)

x = subset(seurat, orig.ident == "PBS")
p1 <- FeaturePlot(x, reduction = "umap", features = c("Lgr5","Ppara","Hmgcs2","Hopx"), ncol = 4)
x = subset(seurat, orig.ident == "BC")
p2 <- FeaturePlot(x, reduction = "umap", features = c("Lgr5","Ppara","Hmgcs2","Hopx"), ncol = 4)

p1 / p2
ggsave("featureplot_umap.pdf", width = 16, height = 7)

# FeaturePlot(seurat, reduction = "tsne", features = c("Dclk1","Trpm5","Chgb","Tac1","Defa24","Defa22","Lyz1",
#                                                      "Ptprc","Tff3","Agr2","Muc2","Ascl2","Fabp1","Alpi","Stmn1","Cd44"))

snn <- 0.8
seurat <- SetIdent(seurat, value = paste0("RNA_snn_res.", snn))

# 找到所有差异表达基因
allmarkers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.table(allmarkers, paste0("cluster/snn.", snn, ".allmarkers.diff.tsv"), sep = "\t", quote = F, row.names = F)

# 提取每个簇中具有最显著差异表达的前5个基因
top5_markers <- allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# 创建基因的Dot Plot
DotPlot(seurat, features = unique(top5_markers$gene)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
  NoLegend()
# ggsave(paste0("cluster/snn.", snn, ".cell.cluster.scatter.pdf"), width = 12, height = 5)

# 创建基因的Heatmap
DoHeatmap(seurat, features = unique(top5_markers$gene)) +
  theme(axis.text.y = ggplot2::element_text(size = 8)) + 
  NoLegend()
# ggsave(paste0("cluster/snn.", snn, ".cell.cluster.heatmap.pdf"), width = 12, height = 5)

# 选择snn 1.0

markergenes <- "Lgr5"
FeaturePlot(seurat, features = markergenes, reduction = "tsne")
ggsave("cluster/tSNE_lgr5.pdf", width = 6, height = 5)

allmarkers %>% 
  filter(gene == "Lgr5") %>%
  mutate(pval = p_val_adj) %>% 
  add_plab() %>% 
  ggbarplot("cluster", "avg_log2FC", fill = "cluster", sort.val = "desc", sort.by.groups = F,
            label = .$plab, lab.col = "red", palette = "Set2", legend = "none", xlab = "Cluster ID",
            ylab = "lgr5 avg_log2FC") +
  coord_fixed(ratio = .75)
ggsave("cluster/lgr5.diff.pdf", width = 4)

allmarkers %>% 
  filter(gene %in% c("Lgr5", "Hmgcs2", "Ppara")) %>% 
  mutate(pval = p_val_adj) %>% 
  add_plab() %>% 
  ggbarplot("gene", "avg_log2FC", fill = "gene", palette = "Set2", legend = "none",
            ylab = "avg_log2FC", x.text.angle = 90, xlab = "") +
  geom_hline(yintercept = 1, lty = "dashed", linewidth = .4) +
  facet_wrap(vars(cluster), nrow = 2) +
  coord_fixed(ratio = 1)
ggsave("cluster/lgr5.diff2.pdf", height = 5, width = 5)

# 细胞数量
names(seurat@meta.data)
DimPlot(seurat, split.by = "sample_title", reduction = "tsne")
ggsave("cluster/cell.count.tSNE.pdf", height = 5, width = 20)

seurat@meta.data %>% 
  select(sample = sample_title, cluster = RNA_snn_res.1) %>%
  group_by(sample, cluster) %>% 
  summarise(n = n()) %>%
  ggbarplot("sample", "n", fill = "sample", palette = "Set2", legend = "none",
            ylab = "cell count", xlab = "", x.text.angle = 90) +
  facet_wrap(vars(cluster), nrow = 3)
ggsave("cluster/cell.count.barplot.pdf", height = 9, width = 12)

# 细胞内差异基因分析
cell_cluster <- subset(seurat, RNA_snn_res.1 == "19")
table(seurat$RNA_snn_res.1)

seurat@meta.data$comparison <- paste0(seurat@meta.data$RNA_snn_res.1, "_", gsub(" ", "_", seurat@meta.data$sample_title))

Idents(seurat) <- "comparison"
x <- FindMarkers(seurat, ident.1 = "19_Crypts-Irradiated", ident.2 = "19_Crypts-Normal", verbose = F) %>% 
  rownames_to_column("name")
# x <- FindMarkers(seurat, ident.1 = "19_Whole_Epithelium-Irradiated", ident.2 = "19_Whole_Epithelium-Normal", verbose = F)


# We mapped receptor-ligand pairs onto cell minor types using CellphoneDB software tool to construct a putative cell-cell interaction network across disease states