library(ggClusterNet)
library(sna)
library(phyloseq)
library(igraph)
library(network)
library(ggplot2)

# 直接输入phyloseq格式的数据
data(ps)

# 按照丰度过滤微生物表格，并计算相关矩阵，按照指定的阈值挑选矩阵中展示的数值。
# 调用了psych包中的corr.test函数，使用三种相关方法。
# N参数提取丰度最高的150个OTU；method.scale参数确定微生物组数据的标准化方式，这里我们选用TMM方法标准化微生物数据。
# 提取丰度最高的指定数量的otu进行构建网络
result <- corMicro(ps = ps, N = 150, method.scale = "TMM", r.threshold = 0.8, p.threshold = 0.05, method = "spearman")

# 提取相关矩阵
cor <- result[[1]]
head(cor)

# 网络中包含的OTU的phyloseq文件提取
ps_net = result[[3]]

#-导出otu表格
otu_table = ps_net %>% vegan_otu() %>% t() %>% as.data.frame()

# 制作分组,我们模拟五个分组
# ggClusterNet中提供了多种优秀的网络可视化布局算法：
# 1. PolygonClusterG：环状模块，环状布局
# 2. PolygonRrClusterG：环状模块，模块半径正比于节点数量，环状布局
# 3. randomClusterG：随机布局，环状模块
# 4. ArtifCluster：环状模块，人工布局
# 5. randSNEClusterG：sna包中的布局按照模块布局
# 6. PolygonModsquareG 环状布局，顺序行列排布
# 7. PolygonRdmNodeCir 实心圈布局，环状布局，控制半径
# 8. model_Gephi.2：模仿Gephi布局
# 9. model_igraph:模仿igraph布局
# 10. model_maptree：按照maptree算法布局模块
# 11. model_maptree2：内聚算法改进离散点排布
# 12. model_filled_circle 实心圆布局，基于分组信息
# 等···········
# 这是网络布局的基础，无论是什么聚类布局，都需要制作一个分组文件，这个文件有两列，一列是节点，一列是分组信息，
# 这个分组信息名称为：group。这个文件信息就是用于对节点进行分组，然后按照分组对节点归类，使用包中可视化函数计算节点位置。
# 注意分组文件的格式，分为两列，第一列是网络中包含的OTU的名字，第二列是分组信息，同样的分组标记同样的字符。
# 人工构造分组信息：将网络中全部OTU分为五个部分，等分
netClu = data.frame(ID = row.names(otu_table), group = rep(1:5,length(row.names(otu_table)))[1:length(row.names(otu_table))])
netClu$group = as.factor(netClu$group)
head(netClu)

# PolygonClusterG 根据分组，计算布局位置坐标
# 不同的模块按照分组聚集成不同的圆，并且圆形的大小一样。如果一个分组只有一个点，则这个点坐落在圆心位置。
result2 = PolygonClusterG(cor = cor, nodeGroup = netClu)
node = result2[[1]]
head(node)
# nodeadd 节点注释的：用otu表格和分组文件进行注释
# nodeadd函数只是提供了简单的用注释函数，用户可以自己在node的表格后面添加各种注释信息。
tax_table = ps_net %>% vegan_tax() %>% as.data.frame()
nodes = nodeadd(plotcord = node, otu_table = otu_table, tax_table = tax_table)
head(nodes)
# 计算边
edge = edgeBuild(cor = cor,node = node)
head(edge)
# 出图
ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)), data = edge, linewidth = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean), pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(aspect.ratio = 1)


# 模拟不同的分组—可视化
# 模拟不同分组效果展示：1个分组
# 这是网络布局的基础，无论是什么聚类布局，都需要制作一个分组文件，这个文件有两列，一列是节点，一列是分组信息，这个分组信息名称必须设定为：group。
netClu = data.frame(ID = row.names(tax_table),group = rep(1,length(row.names(tax_table)))[1:length(row.names(tax_table))])
netClu$group = as.factor(netClu$group)
result2 = PolygonClusterG (cor = cor,nodeGroup =netClu)
node = result2[[1]]
head(node)
# node节点注释
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
# 计算边
edge = edgeBuild(cor = cor,node = node)
# 出图
ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)), data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

# 模拟不同的分组查看效果：8个分组
netClu = data.frame(ID = row.names(cor),group =rep(1:8,length(row.names(cor)))[1:length(row.names(cor))] )
netClu$group = as.factor(netClu$group)
result2 = PolygonClusterG (cor = cor,nodeGroup =netClu )
node = result2[[1]]
# 节点注释
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
# 计算边
edge = edgeBuild(cor = cor,node = node)
# 出图
ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)), data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())


# 按照微生物分类不同设定分组
# 计算相关
result = corMicro (ps = ps, N = 200, method.scale = "TMM", r.threshold=0.8, p.threshold=0.05, method = "spearman")
# 提取相关矩阵
cor = result[[1]]
head(cor)
# 网络中包含的OTU的phyloseq文件提取
ps_net = result[[3]]
# 导出otu表格
otu_table = ps_net %>% vegan_otu() %>% t() %>% as.data.frame()
tax = ps_net %>% vegan_tax() %>% as.data.frame()
tax$filed = tax$Phylum
group2 <- data.frame(ID = row.names(tax),group = tax$Phylum)
group2$group  =as.factor(group2$group)
result2 = PolygonClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]
# 节点注释
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
# 计算边
edge = edgeBuild(cor = cor,node = node)
# 出图
ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)), data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

# 结果发现这些高丰度OTU大部分属于放线菌门和变形菌门，其他比较少。所以下面我们按照OTU数量的多少，对每个模块的大小进行重新调整。
result2 = PolygonRrClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]
# 节点注释
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
# 计算边
edge = edgeBuild(cor = cor,node = node)
# 出图
ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)), data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())


# 用实心点作为每个模块的布局方式1
set.seed(12)
# 实心圆
result2 = model_filled_circle(cor = cor, culxy =TRUE,
                              da = NULL,# 数据框，包含x,和y列
                              nodeGroup = group2,
                              mi.size = 1,# 最小圆圈的半径，越大半径越大
                              zoom = 0.3# 不同模块之间距离
)
node = result2[[1]]
# 节点注释
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
# 计算边
edge = edgeBuild(cor = cor,node = node)
# 出图
ggplot() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)), data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

# 用实心点作为每个模块布局方式2
set.seed(12)
# 实心圆2
result2 = model_maptree_group(cor = cor, nodeGroup = group2)
node = result2[[1]]
# node节点注释
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
# 计算边
edge = edgeBuild(cor = cor,node = node)
# 出图
ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)), data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

# 按照网络模块分析定义分组
netClu  = modulGroup(cor = cor,cut = NULL,method = "cluster_fast_greedy" )
result2 = model_maptree_group(cor = cor, nodeGroup = group2)
node = result2[[1]]
# node节点注释
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
nodes2 = nodes %>% inner_join(netClu,by = c("elements" = "ID"))
nodes2$group = paste("Model_",nodes2$group,sep = "")
# 计算边
edge = edgeBuild(cor = cor,node = node)
# 出图
ggplot() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)), data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

# 使用升级的model_maptree2：不在可以将每个模块独立区分，而是将模块聚拢，并在整体布局上将离散的点同这些模块一同绘制到同心圆内。
set.seed(2)
result2 = model_maptree2(cor = cor, method = "cluster_fast_greedy")
node = result2[[1]]
# node节点注释
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
nodes2 = nodes %>% inner_join(netClu,by = c("elements" = "ID"))
nodes2$group = paste("Model_",nodes2$group,sep = "")
# 计算边
edge = edgeBuild(cor = cor,node = node)
# 出图
ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)), data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())


# model_igraph布局
result = cor_Big_micro(ps = ps, N = 1000,r.threshold=0.6,p.threshold=0.05,method = "spearman")
# 提取相关矩阵
cor = result[[1]]
dim(cor)
result2 <- model_igraph(cor = cor, method = "cluster_fast_greedy", seed = 12)
node = result2[[1]]
head(node)
dat = result2[[2]]
head(dat)
tem = data.frame(mod = dat$model,col = dat$color) %>%  dplyr::distinct( mod, .keep_all = TRUE)  
col = tem$col
names(col) = tem$mod
# node节点注释
otu_table = as.data.frame(t(vegan_otu(ps)))
tax_table = as.data.frame(vegan_tax(ps))
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
# 计算边
edge = edgeBuild(cor = cor,node = node)
colnames(edge)[8] = "cor"
head(edge)
tem2 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color)
head(tem2)
tem3 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color)
head(tem3)

tem4 = tem2 %>%inner_join(tem3)
head(tem4)

edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
                        manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1"))
head(edge2)
col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>% 
  select(color,manual)
col0 = col_edge$manual
names(col0) = col_edge$color

library(ggnewscale)

p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                              data = edge2, size = 1) +
  scale_colour_manual(values = col0) 
p2 = p1 +
  new_scale_color() +
  geom_point(aes(X1, X2,color =model), data = dat,size = 4) +
  scale_colour_manual(values = col) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p2

# 节点模块化计算和可视化
result4 = nodeEdge(cor = cor)
# 提取边文件
edge = result4[[1]]
# 提取节点文件
node = result4[[2]]
igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
p <- res[[1]]
p

# 模块内连通度（Zi）和模块间连通度（Pi）的定义
# 依据节点的的拓扑特征可将节点属性分为4种类型，包括：
# Module hubs（模块中心点，在模块内部具有高连通度的节点，Zi > 2.5且Pi < 0.62），
# Connectors（连接节点，在两个模块之间具有高连通度的节点，Zi < 2.5且Pi > 0.62），
# Network hubs（网路中心点，在整个网络中具有高连通度的节点，Zi > 2.5且Pi > 0.62）、
# Peripherals（外围节点，在模块内部和模块之间均不具有高连通度的节点，Zi < 2.5且Pi <0.62）。

# 基于模块内连通度（Zi）和模块间连通度（Pi）的网络核心节点判别方法得到了广泛的应用。
# 例如，微生物共发生网络一般可以被划分成多个模块，模块是网络中高度连接的区域，
# 模块可能反应了栖息地的异质性、系统发育上亲缘关系较近物种的聚集、生态位的重叠和物种的共进化，
# 被认为是系统发育、进化或功能上独立的单元（Olesen et al. 2007）。
# 在生态网络模块中识别的关键节点，往往代表了在维持微生物群落结构稳定性上可能起重要作用关键物种。
# Shi等（2016）在野燕麦根际微生物共发生网络中通过模块内连通度（Zi）和模块间连通度（Pi）的概念寻找核心微生物物种就是依据此原理。
# (摘抄自生信小白鱼)

# 网络性质计算，22年6月升级后版本包括了16项网络属性，包括周集中老师21年NCC文章中全部属性
dat = net_properties(igraph)
head(dat)
# 升级后包含的网络属性更多
dat = net_properties.2(igraph, n.hub = T)
head(dat, n = 16)


# 节点性质计算
nodepro = node_properties(igraph)
head(nodepro)

# 扩展-关键OTU挑选
# Hub节点是在网络中与其他节点连接较多的节点，Hub微生物就是与其他微生物联系较为紧密的微生物，可以称之为关键微生物（keystone）
hub = hub_score(igraph)$vector %>%
  sort(decreasing = TRUE) %>%
  head(5) %>%
  as.data.frame()
colnames(hub) = "hub_sca"
ggplot(hub) +
  geom_bar(aes(x = hub_sca,y = reorder(row.names(hub),hub_sca)),stat = "identity",fill = "#4DAF4A")

# 对应随机网络构建和网络参数比对
result = random_Net_compate(igraph = igraph, type = "gnm", step = 100, netName = layout)
p1 = result[[1]]
sum_net = result[[4]]
p1
head(sum_net)


# 微生物组网络pipeline分析
# 微生物网络-大网络
# 使用model_maptree2布局计算微生物大网络。
# 大网络运算时间会比较长，这里我没有计算zipi，用时5min完成全部运行。
# N=0，代表用全部的OTU进行计算。
# 3000个OTU不计算zipi全套需要18min。
data("ps16s")
library(WGCNA)
path = "./result_big_1000/"
dir.create(path)
result = network.2(ps = ps16s,
                   N = 1000,
                   big = TRUE,
                   maxnode = 5,
                   select_layout = TRUE,
                   layout_net = "model_maptree2",
                   r.threshold=0.4,
                   p.threshold=0.01,
                   label = FALSE,
                   path = path,
                   zipi = FALSE)

# 多组网络绘制到一个面板
p = result[[1]]
# 全部样本网络参数比对
data = result[[2]]
num= 3
# plotname1 = paste(path,"/network_all.jpg",sep = "")
# ggsave(plotname1, p,width = 16*num,height = 16,dpi = 72)

plotname1 = paste(path,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 10*num,height = 10,limitsize = FALSE)

tablename <- paste(path,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)