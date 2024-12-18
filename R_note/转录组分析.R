# Jinxin Meng, 20241022 --------

# 获取基因长度，根据长度和rc计算基因的fpkm值 -----------------
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # 基因组信息
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb) # 基因长度获取
gene_lengths <- data.frame(ENTREZID = names(genes), length = width(genes))
metadata_gene <- select(org.Hs.eg.db, keys = rownames(profile), keytype = "SYMBOL", # profile是基因表达矩阵
                        columns = c("GENENAME", "ENTREZID", "GENETYPE", "SYMBOL")) %>% 
  filter(GENETYPE == "protein-coding") %>% 
  left_join(gene_lengths, by = "ENTREZID") %>% 
  filter(!is.na(length)) 
source("/code/R_func/transform_rc.R") # rc2fpkm
fpkm <- rc2fpkm(profile_HT, dplyr::select(metadata_gene, name = SYMBOL, len = length))
write.table(fpkm, "profile_HT.fpkm.tsv", sep = "\t", quote = F, row.names = F)

# 转录组差异分析 -------------------
# reads count profile
rc <- read.delim("rc.txt", row.names = 1, check.names = F) %>% 
  dplyr::select(all_of(samples))
rc <- rc[rowSums(rc > 5) >= 2, ]

group <- read.delim("map.txt") %>% 
  filter(sample %in% samples)
group_levels <- c("Pbs.1d","BC.1d")

metadata <- group %>% 
  dplyr::select(sample, group) %>% 
  mutate(sample = factor(sample, colnames(rc))) %>% 
  arrange(sample) %>% 
  column_to_rownames(var = "sample") %>% 
  mutate(group = factor(group, group_levels))

# DESeq_obj
dds <- DESeqDataSetFromMatrix(countData = rc, colData = metadata, design = ~ group)
DEseq_obj <- DESeq(dds)

gp <- c("BC.1d", "Pbs.1d")
diff <- results(DEseq_obj, contrast = c("group", gp)) %>% 
  data.frame(.) %>% 
  rownames_to_column(var = "name")

diff <- results(DEseq_obj, contrast = c("group", gp)) %>% 
  data.frame(.) %>% 
  rownames_to_column(var = "name") %>% 
  mutate(enriched = ifelse(log2FoldChange > 1 & pvalue < 0.05, "up", 
                           ifelse(log2FoldChange < -1 & pvalue < 0.05, "down", "none"))) %>% 
  filter(enriched == "up")

# 转录组差异基因富集分析 GO 和 GSEA (过表达基因富集分析) --------------
genes <- select(org.Mm.eg.db, diff$name, keytype = "SYMBOL", columns = "ENTREZID") # id转换
diff <- read.delim("data/Mock_CW3.vs.Mock.txt") %>%  # 差异分析
  dplyr::select(name, log2FC, pval = PValue, FDR) %>% 
  mutate(enriched = ifelse(log2FC > 0 & pval < 0.05, "up",
                           ifelse(log2FC < 0 & pval < 0.05, "down", "none"))) %>% 
  filter(enriched == "down")

# KEGG集分析
eKEGG <- enrichKEGG(gene = na.omit(genes$ENTREZID), organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1)
eKEGG_out <- data.frame(eKEGG) %>% 
  mutate(geneName = map_vec((strsplit(x = geneID, "/")), \(x) paste(genes$SYMBOL[match(x, genes$ENTREZID)], collapse = "/")), .after = geneID, 
         Description = gsub(" - Mus.*", "", Description))

# GO 富集分析
# library(org.Mm.eg.db)
eGO <- enrichGO(gene = na.omit(genes$ENTREZID), org.Mm.eg.db, pvalueCutoff = 1, qvalueCutoff = 1)
eGO_out <- data.frame(eGO) %>%
  mutate(geneName = map_vec((strsplit(x = geneID, "/")), \(x) paste(genes$SYMBOL[match(x, genes$ENTREZID)], collapse = "/")), .after = geneID)

# 富集分析作图
data <- read.csv("CBCs BC vs PBS eKEGG.csv") %>% 
  dplyr::select(all_of(c("Description","GeneRatio","pvalue"))) %>% 
  mutate(GeneRatio = as.numeric(stringr::str_split_i(GeneRatio, "/", 1))/
           as.numeric(stringr::str_split_i(GeneRatio, "/", 2))) %>% 
  arrange(desc(GeneRatio)) %>% 
  mutate(Description = factor(Description, (Description)))
ggbarplot(data, "Description", "GeneRatio", fill = "pvalue", rotate = T, width = .6, 
          ylab = "Gene ratio", xlab = "", legend = "right",
          title = "KEGG pathways enrichment for up-regulated\ngenes in stem cells after BC-colonization") +
  scale_fill_gradient(low = "#abd9e9", high = "#fdae61")

# clusterProfiler 获取一个物种信息KEGG的相关信息，通路和基因
download_KEGG("mmu")

# 转录组差异基因富集分析 GSEA --------------
# 写了一个python脚本；
# clusterprofiler中有函数做GSEA分析的，但是我没测试明白；
# clusterprofiler中支持GSEA富集分析，可以用函数 gseGO(), gseKEGG(), GSEA()
# 在我的机器上，运行不了GSEA富集分析，查到的原因可能如下：重新安装 fgsea和BiocParallel包；但是最终还是无法运行这个分析。后来发现要在运行代码前register(SerialParam())加这个参数，告诉BioParallel接下来的计算不进行并行处理，采用串行的方式。
gene_list = structure(diff$log2FoldChange, names = diff$name) %>% sort(decreasing = T)
register(SerialParam())
gsea_result <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", minGSSize = 1, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
genes <- select(org.Mm.eg.db, keys(org.Mm.eg.db), keytype = "ENTREZID", columns = "SYMBOL")
diff <- read.delim("data/Mock_CW3.vs.Mock.txt") %>% 
  dplyr::select(name, log2FC, pval = PValue, FDR) %>% 
  left_join(genes, by = c("name" = "SYMBOL")) %>% 
  filter(!is.na(ENTREZID))
gene_list <- structure(diff$log2FC, names = diff$ENTREZID) %>% sort(decreasing = T)
register(SerialParam())
gseKEGG <- gseKEGG(gene_list, organism = "mmu", minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1)
gseKEGG_out <- data.frame(gseKEGG, row.names = NULL) %>% 
  mutate(geneName = map_vec((strsplit(x = core_enrichment, "/")), \(x) paste(genes$SYMBOL[match(x, genes$ENTREZID)], collapse = "/")), .after = core_enrichment, 
         Description = gsub(" - Mus.*", "", Description))

# 自己构建基因集进行GSEA
lib <- data.frame(item = "ketone body metabolism", name = c("AACS","OXCT1","HMGCS2","HMGCS2","ACSS3","TYRP1","HMGCL"))
lib <- read.delim("/database/MSigDB/c7.all.v2024.1.Hs.symbols.lib", header = F, col.names = c("item", "name")) %>% head(300)
gsea_out = GSEA(gene_list, TERM2GENE = lib, pvalueCutoff = 1, maxGSSize = 1000, minGSSize = 11)

# GSEA作图的方法
library(enrichplot)
gseaplot(gsea_out, geneSetID = "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN")
gseaplot2(gsea_out, geneSetID = "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN")
gseaplot2(gsea_out, geneSetID = "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN", color = "red")
gseaplot2(gsea_out, geneSetID = "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN", color = "red", 
          base_size = 20, subplots = 1:2, title = "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN")

# ReactomeGSA通路富集分析 ---------------
library("ReactomeGSA")


# ReactomePA通路富集分析 ---------------
library("ReactomePA") # 也是Y叔写的，牛
genes <- AnnotationDbi::select(org.Mm.eg.db, keys = keys(org.Mm.eg.db), columns = c("ENSEMBL", "ENTREZID", "SYMBOL"))
gseReactome <- map(diff, \(x)
                   x = left_join(x, genes, by = c("name" = "ENSEMBL")) %>% 
                     filter(!is.na(ENTREZID)) %>% 
                     pull(log2FoldChange, ENTREZID) %>%
                     sort(decreasing = T) %>% 
                     gsePathway(organism = "mouse", pvalueCutoff = 1) %>% 
                     data.frame() %>% 
                     mutate(core_enrichment_name = lapply((strsplit(x = core_enrichment, "/")), 
                                                          \(x) genes$SYMBOL[match(x, genes$ENTREZID)]) %>% 
                              sapply(\(x) paste(x, collapse = "/")), .after = core_enrichment)  )
saveRDS(gseKEGG, "diff.gseReactome.rds")
openxlsx::write.xlsx(gseKEGG, "diff.gseReactome.xlsx")

# 批量快速完成差异分析和富集分析 ----------
dds <- DESeqDataSetFromMatrix(countData = profile, colData = metadata, design = ~ group) # 构建dds对象
des <- DESeq(dds)
gps <- list( # 构建分组比较对, "group"是group表的第二列名称
  c("group", "Dub", "PBS"),
  c("group", "Akk", "PBS"),
  c("group", "EF", "PBS"),
  c("group", "BC", "PBS") )
diff <- map(gps, \(x)
            results(des, contrast = x) %>% 
              data.frame() %>%
              rownames_to_column("name") %>% 
              mutate(enriched = ifelse(log2FoldChange > 1 & pvalue < 0.05, "up", 
                                       ifelse(log2FoldChange < -1 & pvalue < 0.05, "down", "none")))
) %>% 
  setNames(map_vec(gps, \(x) paste0(x[2:3], collapse = "_vs_")))
saveRDS(diff, "diff.rds")

BiocParallel::register(BiocParallel::SerialParam())
genes <- AnnotationDbi::select(org.Mm.eg.db, keys = keys(org.Mm.eg.db), columns = c("ENSEMBL", "ENTREZID", "SYMBOL"))
gseKEGG <- map(diff, \(x)
               x = left_join(x, genes, by = c("name" = "ENSEMBL")) %>% 
                 filter(!is.na(ENTREZID)) %>% 
                 pull(log2FoldChange, ENTREZID) %>%
                 sort(decreasing = T) %>% 
                 gseKEGG(organism = "mmu", pvalueCutoff = 1) %>% 
                 data.frame() %>% 
                 mutate(core_enrichment_name = lapply((strsplit(x = core_enrichment, "/")), 
                                                      \(x) genes$SYMBOL[match(x, genes$ENTREZID)]) %>% 
                          sapply(\(x) paste(x, collapse = "/")), .after = core_enrichment) %>% 
                 mutate(Description = gsub(" - Mus.*", "", Description))  )

saveRDS(gseKEGG, "diff.gseKEGG.rds")
openxlsx::write.xlsx(gseKEGG, "diff.gseKEGG.xlsx")

# 火山图 ---------------------
ggscatter(plot_data, "log2FC", "vip", color = "enriched", legend = "right",
          palette = c("#eeacec", "grey", "#21c1dc"), size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text(aes(label = label), size = 1, fontface = "italic")
ggscatter("log2FC", "pval2", color = "enriched", legend = "right",
                palette = c(up = "#fdae61", none = "grey", down = "#abd9e9"), size = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      ggrepel::geom_text_repel(aes(label = name), size = 4, fontface = "italic")

colors <- structure(c("#56B4E9", "#E69F00", "grey"), names = c("NBP", "CP", "none"))
ggplot(plot_data, aes(estimate, -log10(padj), color = enriched)) +
  geom_point(size = 2, alpha = .8) +
  geom_text(aes(label = label), color = "black", size = 1) +
  scale_color_manual(values = colors, labels = labels) +
  scale_y_continuous(expand = c(.01, .01)) +
  labs(x = "Estimate of Meta-analysis", y = "-Log10-transformed P.adj", color = "Enriched in") + 
  geom_hline(yintercept = -log10(0.05), lty = 2, lwd = .4, color = "black") +
  theme_pubr() +
  theme(aspect.ratio = 0.9)

ggplot(plotdat, aes(x = log2FC, y = -log10(padj), 
                    fill = factor(enriched, levels = c("Control", "None", "Disease")))) +
  geom_point(aes(size = lab), alpha = .7, shape = 21, color = "black", stroke = .2) + 
  geom_vline(xintercept = c(-2, 2), lty = "dashed", lwd = .4, color = "black") + 
  geom_hline(yintercept = -log10(0.05), lty = "dashed", lwd = .4, color = "black") +
  scale_fill_manual(values = c("#0c84c6", "#e0e0e0", "#f7ac00"),
                    labels = c(paste0("Control (", sum(diff$enriched == "Control"), ")"),
                               paste0("None(", sum(diff$enriched == "None"), ")"),
                               paste0("Disease"," (", sum(diff$enriched == "Disease"), ")"))) + 
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), 
                     labels = c(-6, -4, -2, 0, 2, 4, 6), expand = c(.02, .02)) +
  scale_y_continuous(expand = c(.02, .02)) +
  labs(y = "-log10(padj)", x = "log2FoldChange", fill = "Enriched in", size = "Avg. abundance") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black", size = 9),
        axis.title = element_text(color = "black", size = 10), 
        panel.background = element_rect(color = "black", size = .8),
        panel.grid = element_blank(),
        aspect.ratio = 3/4)

res %>% mutate(padj = ifelse(padj < 1e-10, 1e-10, padj), 
                log2FoldChange = ifelse(log2FoldChange < -2.5, -2.5, log2FoldChange),
                log2FoldChange = ifelse(log2FoldChange > 2.5, 2.5, log2FoldChange),
                enriched = factor(enriched, levels = c("Up", "None", "Down"))) %>% 
  ggplot(aes(log2FoldChange, -log10(padj), color = enriched)) +
  geom_point(size = .7) +
  scale_color_manual(values = c("#f46d43", "grey75", "#66bd63"),
                      labels = c(paste0("Up (n=", sum(res$enriched == "Up"), ")"),
                                paste0("None (n=", sum(res$enriched == "None"), ")"),
                                paste0("Down (n=", sum(res$enriched == "Down"), ")"))) + 
  scale_y_continuous(expand = c(.01, .01)) +
  labs(y = expression(-log[10]*P.adj), x = expression(log[2]*FoldChange), color = "",
        title = paste0("Volcano plot  ", paste0(gp[[i]][-1], collapse = " vs.")), 
        subtitle = 'Coefficient: '~~italic('padj')~'<'~'0.05') +
  geom_vline(xintercept = 0, lty = 2, lwd = .3, color = "black") + 
  geom_hline(yintercept = -log10(0.05), lty = 2, lwd = .3, color = "black") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.line = element_blank(),
        plot.title = element_text(size = 10, color = "black"),
        plot.subtitle = element_text(size = 8, color = "black"),
        panel.background = element_rect(linewidth = .4, color = "black"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.background = element_blank(),
        legend.key.height = unit(6, "mm"),
        legend.spacing.x = unit(.2, "mm"),
        legend.title.align = .5, 
        aspect.ratio = 3/4)

# GSVA
# To elucidate the mechanism underlying the effects of p-coumaric acid on cell death in the livers of BDL mice, single-sample gene
# set enrichment analysis (ssGSEA) was employed. Mice were oral administration with p-coumaric acid for 1 h followed by BDL surgery.
# Liver samples were collected on day 4 after BDL. The activity scores for six cell death signaling pathways were obtained using
# the ‘‘GSVA’’ package (v. 1.44.2)63 in R (v 4.2.0). The gene gets for the cell death signaling pathways were obtained from the MsigDB
# (https://www.gsea-msigdb.org/) and KEGG (https://www.kegg.jp/) databases.