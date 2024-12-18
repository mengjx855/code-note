## MaAslin2
- MaAsLin2包可以有效的确定`表型`、`环境`、`暴露`、`协变量`和`微生物宏组学特征`之间的多变量关联。MaAsLin2依靠“通用线性模型”来适应大多数现代流行病学研究设计，包括横断面和纵向的研究。并提供各种数据探索、规范化和转换方法。
- 引用: Multivariable Association Discovery in Population-scale Meta-omics Studies. 2021, PLoS Computational Biology
- 官网：https://huttenhower.sph.harvard.edu/maaslin  

### 1. Bioconductor安装
```R {.line-numbers}
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Maaslin2")
```

### 2.输入文件
MaAsLin2需要输入两个文件，微生物物种和功能的特征表和样本的metadata。不要求两个文件要对应起来，样本可以是行可以是列；软件会自动合并两个表共有的样本。以下官方给的例子，就是后边说的HMP2数据。
```R {.line-numbers}
input_data = system.file(
    "extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file
input_data
input_metadata = system.file(
    "extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file
input_metadata
```

### 运行MaAsLin2
使用上述导入的数据运行MaAsLin2s，建立多变量回归模型去检测微生物物种的丰度和IBD诊断以及菌群紊乱得分之间的关联(`fixed_effects = c("diagnosis", "dysbiosis")`)。输出在当名前工作目录下为`demo_Output`的文件夹中。
```R {.line-numbers}
fit_data = Maaslin2(
    input_data = input_data, 
    input_metadata = input_metadata, 
    output = "demo_output", 
    fixed_effects = c("diagnosis", "dysbiosis"))
```
几种不同类型的统计模型可用于 MaAsLin 中的关联测试（例如简单线性、零膨胀、基于计数等）。我们在相关手稿中评估了它们的优缺点，虽然默认值通常适合大多数分析，但用户可能希望在某些情况下选择不同的模型。 对于许多不同的微生物群落数据类型（分类或功能profile）、环境（人类或其他）和测量（计数或相对比例）来说都是如此，只要使用替代模型和适当修改的标准化/转换方案。
- 对于模型，如果您的输入是计数，则可以使用 NEGBIN 和 ZINB，而对于非计数（例如百分比、CPM 或相对丰度）输入，您可以使用 LM 和 CPLM。 
- 在 MaAsLin2 中实现的标准化方法中，TMM 和 CSS 仅适用于计数，并且与 TSS 和 CLR 不同，它们还返回标准化计数。 因此，如果您的输入是计数，则可以使用上述两种标准化（即 TMM、CSS 或 NONE（如果数据已标准化）），而无需进一步转换（即，transform = NONE）。
- 在非计数模型中，CPLM 要求数据为正。 因此，任何产生负值的转换通常不适用于 CPLM。
- 所有非 LM 模型都使用内部的基于log转换，因为它们与 GLM 密切相关，并且建议在 Transform = NONE 的情况下运行。
- 除此之外，LM 是唯一同时适用于正值和负值（遵循标准化/转换）的模型，并且（根据我们的手稿），它通常对参数变化更加稳健（这通常仅限于非 LM 楷模）。 关于是否使用 LM、CPLM 或其他模型，直观上看，CPLM 或零膨胀替代方案在零存在的情况下应该表现更好，但根据我们的基准测试，我们没有证据表明 CPLM 在实践中明显优于 LM 。

| Model (-m ANALYSIS_METHOD)|Data type|Normalization (-n NORMALIZATION)|Transformation (-t TRANSFORM)|
| - | - | - | - |
|LM|count and non-count|TSS,CLR, NONE|LOG, LOGIT, AST,NONE|
|CPLM|count and non-count|TSS, TMM, CSS, NONE|NONE|
|NEGBIN|count|TMM, CSS, NONE|NONE|
|ZINB|count|TMM, CSS, NONE|NONE|
*[TMM]: trimmed mean of M-values; [TSS]: Total Sum Normalization：通过将样本中每个OTU的读数除以该样本中的总读数来将数据转换成比例;LOGIT:Logistic标准化是借助Logistic函数来对 x 进行非线性映射,经过变换后的f取值区间为：[0,1];*


### MaAsLin2输出
这个例子的是参考知乎大佬的帖子。[宏基因组下游分析之Maaslin2差异分析：你总有机会用的到(另有超大宏基因组数据彩蛋)](https://zhuanlan.zhihu.com/p/582619073)
1. 显著关联：查看`significant_results.tsv`
- value列：对于分类特征，报告的是我们选择的特征的因子的系数和显著性。
举个例子，`Y ~ β0 + β1 * X1 + β2 * X2 + υ`，这个是一个多元线性回归模型的表达式；如果β2是分类变量，样本点∈{男，女}。解释这个模型：当我们选择基数为男时，β0、β1、υ不变的情况下，男的Y多β2个单位。
- coef列：模型的系数值（效应值）；分类变量的系数表示值指定的类别与参考类别之间的对比。默认情况下，maaslin2将字母顺序中的第一个类别设置为引用。
- stderr列：模型的标准误差。
- N：用于构建模型的样本数量。
- N.not.0：非0的样本数量。
- pval：模型的显著性。
- qval：校正后的p值，默认使用BH矫正方法。
2. MaAsLin2也会生成显著模型的可视化图标（分类变量生成箱线图，连续变量生成散点图）。 
3. 分类变量默认用的是分类变量的第一个因子，也可以在数据中修改因子的顺序，也可以用`reference`设置参考水平。
4. MaAsLin2不支持相互作用的变量。这种情况下生成子自变量。
```R {.line-numbers}
df_input_metadata$CD_dysbiosis = (df_input_metadata$diagnosis_modified == "CD") *
                                 df_input_metadata$dysbiosis
df_input_metadata$UC_dysbiosis = (df_input_metadata$diagnosis_modified == "UC") *
                                 df_input_metadata$dysbiosis
fit_data5 = Maaslin2(
    input_data = df_input_data, 
    input_metadata = df_input_metadata, 
    output = "demo_output5", 
    fixed_effects = c("diagnosis_modified", "dysbiosis", "CD_dysbiosis", "UC_dysbiosis"))
```
5. 随机效应（`random effect`）去解决样本间非独立性问题。==不太明白，等听完线性回归补上这块内容==
6. 最小丰度（`min_abundance`）设置feature的最小丰度，默认为0。
7. 最小流行率（`min_prevalence`）设置feature的流行率，默认为0.1。
8. 最小显著性（`max_significance`）设置qval阈值，默认为0.25。
9. 归一化（`normalization`）设置归一化方法，默认为TSS。一般我们用的相对丰度的就不用标准化，用"None"。
10. 数据转化(`transform`）设置数据转化的方法，默认为log。
11. 分析方法（`analysis_method`）设置模型方法，默认是线性模型。
12. 矫正方法（`correction`）设置q值矫正方法，默认是BH。
13. 标准化（`standardize`）默认设置Z值标准化。
14. 热图（`plot_heatmap`）
15. 热图展示的特征数量（`heatmap_first_n`）
16. 散点图（`plot_scatter`）
17. 并行数（`cores`）

### 其他例子
```R {.line-numbers}
library("Maaslin2")

# 按照疾病、样本位置筛选，同时去除没有BMI指数和没有性别的行
metadata <- curatedMetagenomicData::sampleMetadata %>% 
  filter(str_detect(disease, "T2D") &
             body_site == "stool" & 
             study_condition == "T2D" & 
             is.na(BMI) != TRUE &
             is.na(gender) != TRUE)

metadata <- metadata %>% 
  group_by(study_name, subject_id) %>% 
  filter(row_number() == 1) %>% 
  ungroup()

datasets_tokeep <- metadata %>%
  select(study_name, gender) %>%
  group_by(study_name) %>%
  summarise(n_males = sum(gender=="male"),
            n_females= sum(gender=="female"),
            N=n()) %>%
  mutate(keep = (pmin(n_males,n_females) >= 40) & (n_females/N >= 0.25) & (n_males/N >= 0.25)) %>%
  filter(keep == TRUE)

datasets_tokeep <- datasets_tokeep$study_name

metadata <- metadata %>%
  filter(study_name %in% datasets_tokeep)

# Retrieve queried samples
# 获取菌群丰度数据，该包是先设置metadata，然后根据metadata下载你的样本菌群丰度数据，数据类型和行名是根据自己的需要进行修改
# 行名可以选择long和short，short的话行名就是菌种的名称，long会包含界门纲目科属种的分类信息。
tse <- curatedMetagenomicData::returnSamples(metadata, dataType = "relative_abundance", rownames = "short")

metadata_T2D <- data.frame(colData(tse))
metadata_T2D <- metadata_T2D %>%
  rownames_to_column(var = "sample")

write_tsv(metadata_T2D, "Data/MaAslin2/metadata_T2D.tsv")

relative_abundance_T2D <- data.frame(t(assay(tse)))
relative_abundance_T2D <- relative_abundance_T2D %>% 
  rownames_to_column(var = "sample")

write_tsv(relative_abundance_T2D, "relative_abundance_T2D.tsv")

# Maaslin2差异分析
# 需要两个文件一个是菌群丰度文件一个是metadata文件
metadata <- read.table("Data/MaAslin2/metadata_T2D.tsv", header = T, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
data <- read.table("Data/MaAslin2/relative_abundance_T2D.tsv", header = T, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

fit_data2 <- Maaslin2(
  input_data = data, 
  input_metadata = metadata, 
  min_prevalence = 0,
  normalization = "NONE",
  output = "Data/MaAslin2/demo_output", 
  fixed_effects = c("gender"),
  random_effects = c("BMI"),
  reference = c("gender, female"))
```
