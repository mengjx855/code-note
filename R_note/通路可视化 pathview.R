# Jinxin Meng, 20240801, 20240829 ---------------------
# https://academic.oup.com/nar/article/45/W1/W501/3804420
# https://bioconductor.org/packages/release/bioc/html/pathview.html
# https://cloud.tencent.com/developer/article/1539928
# Pathview主要用于可视化pathway图上的数据。
# pathview可以生成KEGG视图和Graphviz视图，前者将用户数据呈现在原生KEGG pathway图上，更自然，更易于阅读。
# 后者使用Graphviz引擎对pathway图进行布局，可以更好地控制节点或边缘属性和pathway拓扑。

library(pathview)
data(gse16873.d)
head(gse16873.d)

# pathway上色图例的限制默认是1，使用limit函数修改
limit = list(gene = 5)


# 先看下Pathview最常见的一种用法：将某个样本的表达量（读入的数据需要是归一化后的表达量）；
# 其实也可以任何列数据，不仅仅是表达量数据，也可以是fold change, p-value, 组蛋白修饰数据等，
# 人为指定的数值型数据也行 (关键是要懂需要展示什么数据、说明什么问题，原理最重要，就像GSEA基因集富集分析也是一样)；
# 最后以color bar的形式在KEGG通路图上的对应节点 (一定注意节点名字的匹配)展示；
# 如下例子所示，我们通过指定gene.data和pathway.id来观察单个样本在典型信号通路细胞周期上的表达变化。
# 该基因芯片是在人体组织上进行的，因此species=“hsa”。
pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa", out.suffix = "gse16873", kegg.native = T)

# 该例子中的图只有一个单层，在原始图层修改节点颜色，保留原始KEGG节点标签 (节点名)。
# 这样输出的文件大小与原始的KEGG PNG文件一样小，但是计算时间相对较长。
# 如果我们想要一个快速的视图，并且不介意将输出文件大小，我们可以通过same.layer = F使用两层图。
# 通过这种方式，节点颜色和标签被添加到原始KEGG的额外图层上。原来的KEGG基因标签（或EC编号）被替换为官方基因符号。
pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
         out.suffix = "gse16873.2layer", kegg.native = T, same.layer = F)

# 上面的两个例子中，我们查看了KEGG pathway图的数据，在KEGG图上我们可以得到所有注释和meta-data，因此数据更具可读性和解释性。
# 然而输出的是PNG格式的栅格图像。我们也可以使用Graphviz engine重新绘制pathway图来查看数据，
# 这样我们对节点和边缘属性能有更多的控制，更重要的是可以保存为PDF格式的矢量图像。
pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
         out.suffix = "gse16873", kegg.native = F, sign.pos="bottomleft")

# 该图的主图和图例都在一个图层或者说一个页面中，我们只列出KEGG边的类型，忽略图例中的节点类型，以节省空间。
# 如果我们想要完整的图例，我们可以使用两个层来创建Graphviz视图: 第1页是主图，第2页是图例。
# 注意，对于Graphviz视图 (PDF文件)，“层”的概念与KEGG视图 (PNG文件)略有不同。
# 在这两种情况下，我们都为两层图设置参数same.layer=F。
pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
         out.suffix = "gse16873.2layer", kegg.native = F, sign.pos="bottomleft", same.layer = F)

# 图形布局样式调整
# 在Graphviz视图中，我们对图形布局有更多的控制，比如可以将节点组拆分为独立的节点，甚至可以将多基因节点扩展为单个基因。
# 分裂的节点或扩展的基因可能从未分裂的组或未扩展的节点继承边。这样我们就可以得到一个基因/蛋白-基因/蛋白相互作用网络。
# 可以更好地查看网络特性(模块化等)和基因方面(而不是节点方面)的数据。
# 注意在KEGG视图中，一个基因节点可能代表多个功能相似或重复的基因/蛋白。成员基因的数量从1到几十不等。
# 为了更好的清晰度和可读性，一般将它们作为路径图上的单个节点放在一起。
# 因此，默认情况下，我们不分割节点并单独标记每个成员基因。
# 但是，我们可以通过总结基因数据来可视化节点数据，用户可以使用node.sum参数指定摘要方法。
pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
         out.suffix = "gse16873.split", kegg.native = F, sign.pos="bottomleft",
         split.group = T)

# 下面的图更复杂了，对简单通路适用，复杂通路就头秃了！！
pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
         out.suffix = "gse16873.split.expanded", kegg.native = F, sign.pos="bottomleft",
         split.group = T, expand.node = T)

# 数据整合
# Pathview为数据集成提供了强大的支持。
# 它可以用来整合、分析和可视化各种各样的生物数据：基因表达、遗传关联、代谢产物、基因组数据、文献和其他可映射到通路的数据类型。
# 当数据映射到KEGG ortholog pathways时，它可以直接用于宏基因组、微生物组或未知物种的数据。

# 化合物和基因集同时绘制在通路上
# 在上面的例子中，我们查看了具有典型的信号通路的基因数据。
# 有时候我们也想研究代谢通路。除了基因节点外，这些通路还有复合节点。
# 因此，我们可以将基因数据和化合物数据与代谢途径进行整合或可视化。
# 这里的基因数据是一个广泛的概念，包括基因、转录本、蛋白质、酶及其表达、修饰和任何可测量的属性。
# 化合物数据也是如此，包括代谢物、药物、它们的测量值和属性。这里我们仍然使用乳腺癌微阵列数据集作为基因数据。
# 然后生成模拟的化合物或代谢组数据，并加载适当的化合物ID类型(具有足够数量的惟一条目)进行演示。
sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000)
data(cpd.simtypes)
head(sim.cpd.data)

# 我们生成了一个包含基因数据和化合物数据的KEGG视图。
# pathview生成的代谢通路图与原始KEGG图相同，只是为了更好地查看颜色，将复合节点放大。
pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
         pathway.id ="00640", species = "hsa", out.suffix = "gse16873.cpd",
         keys.align = "y", kegg.native = T, key.pos = "topright")

# 我们还生成了相同pathway和数据的Graphviz视图。
# Graphviz视图更好地显示了层次结构。
# 对于代谢通路，解析xml文件中的反应条目，并将其转换为基因和复合节点之间的关系。
# 对复合节点使用省略号。标签是从CHEMBL数据库中检索到的标准化合物名称 (KEGG在pathway数据库文件中没有提供它)。
# 化学名称是长字符串，我们需要对它们进行换行，以使其符合图上指定的宽度。
pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data, pathway.id ="00640",
         species = "hsa", out.suffix = "gse16873.cpd", keys.align = "y", kegg.native = F,
         key.pos = "bottomright", sign.pos ="topright", cpd.lab.offset =-0.8)

# 多状态或样本同时或分开绘制
# 在前面的所有示例中，我们查看了单个样本数据，这些数据要么是向量，要么是单列矩阵。
# Pathview还可以处理多个样本数据，用于为每个样本生成图形。
# 从1.1.6版开始，Pathview就可以整合并绘制多状态或样本到一个图中。
set.seed(10)
sim.cpd.data2 = matrix(sample(sim.cpd.data, 18000,  replace = T), ncol = 6)
rownames(sim.cpd.data2) = names(sim.cpd.data)
colnames(sim.cpd.data2) = paste("exp", 1:6, sep = "")
head(sim.cpd.data2, 3)

# KEGG view
p <- pathview(gene.data = gse16873.d[, 1:3],  cpd.data = sim.cpd.data2[, 1:2],
              pathway.id = "00640",  species = "hsa", out.suffix = "gse16873.cpd.3-2s",
              keys.align = "y", kegg.native = T, match.data = F, multi.state = T,
              same.layer = T)

head(p$plot.data.cpd)
#基因表达和化合物都展示前3个样品，一一对应。
p <- pathview(gene.data = gse16873.d[, 1:3],  cpd.data = sim.cpd.data2[, 1:2],
              pathway.id = "00640",  species = "hsa", out.suffix = "gse16873.cpd.3-2s.match",
              keys.align = "y", kegg.native = T, match.data = T, multi.state = T, same.layer = T)

# 同样，也可以使用graphviz view。
#graphviz view
p <- pathview(gene.data = gse16873.d[, 1:3],  cpd.data = sim.cpd.data2[, 1:2],
              pathway.id ="00640",  species = "hsa", out.suffix = "gse16873.cpd.3-2s",
              keys.align = "y", kegg.native = F, match.data = F, multi.state = T,
              same.layer = T,  key.pos = "bottomright", sign.pos = "topright")

# 同样，我们可以选择分别绘制样本，即每个图形一个样本。
# 请注意，在这种情况下，必须匹配gene.data和cpd.data中的样本，以便将其分配给同一图表。
# 因此，参数match.data在这里并不是很有效。(图3就没有化合物的映射了)
#plot samples/states separately
p <- pathview(gene.data = gse16873.d[, 1:3],  cpd.data = sim.cpd.data2[, 1:2],
              pathway.id = "00640", species = "hsa", out.suffix = "gse16873.cpd.3-2s",
              keys.align = "y", kegg.native = T, match.data = F, multi.state = F, same.layer = T)

# 如上所述，同一层上的KEGG视图可能需要一些时间。
# 同样，如果我们不介意丢失原始的KEGG基因标签(或EC编号)，我们可以选择使用两层的KEGG视图来加快这个过程。
p <- pathview(gene.data = gse16873.d[, 1:3],  cpd.data = sim.cpd.data2[, 1:2],
              pathway.id = "00640",  species = "hsa", out.suffix = "gse16873.cpd.3-2s.2layer",
              keys.align = "y", kegg.native = T, match.data = F, multi.state = T,  same.layer = F)

# 离散数据标记上下调或是否存在
# 到目前为止，我们一直在处理连续数据。但我们也经常处理离散数据。
# 例如，我们根据一些统计数据(p值、折叠变化等)选择显著基因或化合物列表。
# 输入数据可以命名为两个层次的向量，1或0(显着或不显着)，也可以是一个更短的显着基因/化合物名称列表。
# 在接下来的两个例子中，我们只使gene.data和cpd.data或gene.data离散。

require(org.Hs.eg.db)
gse16873.t <- apply(gse16873.d, 1, function(x) t.test(x, alternative = "two.sided")$p.value)
sel.genes <- names(gse16873.t)[gse16873.t < 0.1]
sel.cpds <- names(sim.cpd.data)[abs(sim.cpd.data) > 0.5]

# 我们分别查看下sel.genes和sel.cpds的数据结构：
# 选择高亮的基因
head(sel.genes)
head(sel.cpds)
pathview(gene.data = sel.genes, cpd.data = sel.cpds,  pathway.id ="00640",
         species = "hsa", out.suffix = "sel.genes.sel.cpd", keys.align = "y", kegg.native = T,
         key.pos = "topright",  limit = list(gene = 1, cpd = 1), bins = list(gene = 1, cpd = 1),
         na.col = "gray", discrete = list(gene = T, cpd = T))

pathview(gene.data = sel.genes, cpd.data = sim.cpd.data,  pathway.id = "00640",
         species = "hsa", out.suffix = "sel.genes.cpd",  keys.align = "y", kegg.native = T,
         key.pos = "topright",  limit = list(gene = 1, cpd = 1), bins = list(gene = 1, cpd = 10),
         na.col = "gray", discrete = list(gene = T, cpd = F))

# 不同来源的ID的转换和映射
# pathview的一个显著特点是它强大的ID映射能力。
# 集成的Mapper模块将10多种基因或蛋白ID、20多种化合物或代谢物ID映射到标准KEGG基因或化合物ID。
# 换句话说，使用这些不同ID类型命名的用户数据可以精确地映射到目标KEGG路径。
# Pathview适用于大约4800个物种的路径，物种可以以多种格式指定:KEGG代码、科学名称或常用名称。
cpd.cas <- sim.mol.data(mol.type = "cpd", id.type = cpd.simtypes[2],  nmol = 10000)
gene.ensprot <- sim.mol.data(mol.type = "gene", id.type = gene.idtype.list[4], nmol = 50000)

head(cpd.simtypes)
head(cpd.cas)
head(gene.ensprot)
# 注意参数中的gene.idtype和cpd.idtype，用来指定基因和化合物的ID类型。
pathview(gene.data = gene.ensprot, cpd.data = cpd.cas,  gene.idtype = gene.idtype.list[4],
         cpd.idtype = cpd.simtypes[2],  pathway.id = "00640", species = "hsa", same.layer = T,
         out.suffix = "gene.ensprot.cpd.cas", keys.align = "y", kegg.native = T,
         key.pos ="bottomright", sign.pos = "topright", limit = list(gene = 3, cpd = 3),
         bins = list(gene = 6, cpd = 6))

# 对于不在自动映射列表中的外部ID，我们可以使用mol.sum函数 (也是Mapper模块的一部分)显式地进行ID和数据映射。
# 这里我们需要提供id.map (外部ID和KEGG标准ID之间的映射矩阵)。
# 我们使用ID映射函数id2eg和cpdidmap等来得到id.map矩阵。注意，这些ID映射函数可以独立于pathview 主函数使用。
id.map.cas <- cpdidmap(in.ids = names(cpd.cas), in.type = cpd.simtypes[2],
                       out.type = "KEGG COMPOUND accession")
cpd.kc <- mol.sum(mol.data = cpd.cas, id.map = id.map.cas)
id.map.ensprot <- id2eg(ids = names(gene.ensprot), category = gene.idtype.list[4], org = "Hs")
gene.entrez <- mol.sum(mol.data = gene.ensprot, id.map = id.map.ensprot)
p <- pathview(gene.data = gene.entrez, cpd.data = cpd.kc, pathway.id = "00640",
              species = "hsa", same.layer = T, out.suffix = "gene.entrez.cpd.kc", keys.align = "y",
              kegg.native = T,  key.pos ="bottomright", sign.pos = "topright",
              limit = list(gene = 3, cpd = 3), bins = list(gene = 6, cpd = 6))

# 不同物种使用时名称的处理
# 当对KEGG处理时，物种是一个棘手但重要的问题。
# KEGG拥有自己的物种专用命名法和数据库，其中包括大约4800个基因组完整的生物体。
# 换句话说，这些生物体的基因数据都可以通过pathview进行映射、可视化和分析。
# 这种全面的物种覆盖是pathview数据集成能力的一个突出特点。
# 然而，KEGG并不以同样的方式对待所有这些生物体/基因组。
# KEGG使用NCBI GeneID(或Entrez基因)作为最常见的模型动物的默认ID，包括人类、小鼠、大鼠等，以及一些作物，如大豆、葡萄和玉米。
# 另一方面，KEGG对所有其他生物体使用 Locus 标记和其他id，包括动物、植物、真菌、原生生物以及细菌和古生菌。
# Pathview带有一个数据矩阵korg，其中包括支持的KEGG物种数据和默认基因ID的完整列表。
# 让我们探索korg数据矩阵，以便对KEGG物种及其Gene ID的使用有所了解。
data(korg)
head(korg)
#number of species which use Entrez Gene as the default ID
sum(korg[,"entrez.gnodes"]=="1",na.rm=T)  #204

#number of species which use other ID types or none as the default ID
sum(korg[,"entrez.gnodes"]=="0",na.rm=T)  #5041

#new from 2017: most species which do not have Entrez Gene annotation any more
na.idx=is.na(korg[,"ncbi.geneid"])
sum(na.idx)  #4674

# 从上面的korg的探索中，我们看到4800个KEGG物种中，只有少数没有NCBI（Entrez）基因ID或基因ID（注释）其中的一个。
# 几乎所有物种都具有默认的KEGG基因ID（通常是Locus标签）和Entrez Gene ID注释。
# 因此，pathview接受所有这些物种的gene.idtype =“kegg”或“Entrez”（不区分大小写）。
# 用户需要确保正确的gene.idtype是特定的，因为除了35种常见物种外，KEGG和Entrez基因ID不同。
# 对于19种，Bioconductor提供基因注释包。
# 用户可以自由地输入gene.data和其他gene.idtype。有关这些Bioconductor注释物种和额外基因ID类型的列表允许：
data(bods)
bods
data(gene.idtype.list)
gene.idtype.list

# 所有先前的例子用了人类数据，其中Entrez Gene是KEGG的默认基因ID。
# Pathview现在（从版本1.1.5开始）显式处理使用Locus标签或其他基因ID作为KEGG默认ID的物种。
# 以下是大肠杆菌菌株K12数据的几个例子。
# 首先，我们使用大肠杆菌K12的默认KEGG ID（基因座标签）处理基因数据。
eco.dat.kegg <- sim.mol.data(mol.type="gene",id.type="kegg",species="eco",nmol=3000)
head(eco.dat.kegg)

p <- pathview(gene.data = eco.dat.kegg, gene.idtype="kegg",  pathway.id = "00640",
              species = "eco", out.suffix = "eco.kegg", kegg.native = T, same.layer=T)

# 我们也可以用与人类相同的方法对大肠杆菌K12的Entrez Gene ID进行基因数据处理。
eco.dat.entrez <- sim.mol.data(mol.type="gene",id.type="entrez",species="eco",nmol=3000)

# 查看部分eco.dat.entrez数据：
head(eco.dat.entrez)
p <- pathview(gene.data = eco.dat.entrez, gene.idtype="entrez",  pathway.id = "00640",
              species = "eco", out.suffix = "eco.entrez", kegg.native = T, same.layer=T)
# 基于上述bods数据，大肠杆菌K12在Bioconductor中有注释。
# 因此，我们可以进一步研究其基因数据与其他ID类型，例如，官方基因符号。
# 在调用pathview时，首先需要将这些数据映射到Entrez Gene ID（如果还没有），然后再映射到默认的KEGG ID（Locus标签）。
# 因此，它比上一个例子花费更长的时间。
egid.eco=eg2id(names(eco.dat.entrez), category="symbol", pkg="org.EcK12.eg.db")
eco.dat.symbol <- eco.dat.entrez
names(eco.dat.symbol) <- egid.eco[,2]
head(eco.dat.kegg)
p <- pathview(gene.data = eco.dat.symbol, gene.idtype="symbol", pathway.id = "00640", species = "eco", out.suffix = "eco.symbol.2layer", kegg.native = T, same.layer=F)

# 未注释物种的处理 （直接用于宏基因组或微生物组数据）
# 重要的是，当数据被映射到KEGG ortholog pathways时，pathview可以直接用于宏基因组或微生物组数据。
# 来自KEGG中未注释和未包含的任何新物种（非KEGG物种）的数据也可以通过pathview用同样的方法映射到KEGG ortholog pathways中进行分析和可视化。
# 在下一个例子中，我们首先模拟映射的KEGG ortholog基因数据。
# 然后将数据作为gene.data输入，其中species =“ko”。
ko.data=sim.mol.data(mol.type="gene.ko", nmol=5000)
head(ko.data)
p <- pathview(gene.data = ko.data, pathway.id = "04112", species = "ko", out.suffix = "ko.data", kegg.native = T)