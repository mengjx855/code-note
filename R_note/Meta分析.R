####info####
# Creator：Jinxin Meng
# Created Date:2022-7-19
# Modified Date:2022-9-2


# arcsine-square root-transformed
otu <- asin(sqrt(otu))
otu <- t(otu) %>% data.frame(check.names = F)

# meta input data
# meta-analysis for each family among all study
# use measure method of standardized mean difference
# random effect model applied for meta-analysis
smd_meta_dt <- rbind()
smd_meta_model_dt <- rbind()
for (i in seq_len(ncol(otu))) {
  # prepare data. 准备数据
  taxa <- colnames(otu)[i]
  
  otu_i <- otu %>% 
    select(taxa = all_of(taxa)) %>%
    rownames_to_column(var = "sample") %>% 
    merge(x = ., y = group, by = "sample") %>% 
    select(-sample, -case)
    
  # generate meta-analysis input. 产生输入meta分析的数据
  in_case <- otu_i %>%  
    subset(group%in%"Disease") %>% 
    group_by(study, group) %>% 
    summarise(Mean = mean(taxa), Sd = sd(taxa), N = n()) %>% 
    ungroup() %>% 
    select(-group) %>% 
    rename(d_Mean = Mean, d_Sd = Sd, d_N = N)
  
  in_ctr <- otu_i %>%  
    subset(group%in%"Control") %>% 
    group_by(study, group) %>% 
    summarise(Mean = mean(taxa), Sd = sd(taxa), N = n()) %>% 
    ungroup() %>% 
    select(-group) %>% 
    rename(c_Mean = Mean, c_Sd = Sd, c_N = N)
  
  in_meta <- merge(in_case, in_ctr, by = "study")
  
  # Calculate effect size and variance in each project. 计算效应值和案例内方差
  # We select the method of standardized mean difference provided by Hedges. 使用Hedges提供的SMD方法计算效应值（yi）和案例内方差（vi)
  smd_meta <- escalc(measure = "SMD", data = in_meta, append = T,
                     m1i = d_Mean, m2i = c_Mean, 
                     sd1i = d_Sd, sd2i = c_Sd, 
                     n1i = d_N, n2i = c_N)
  
  # Calculate cumulative effect size using Random-effect model. 计算累积效应值，使用随机效应模型，随机效应模型除了随机因素引起的误差外，还考虑一些案例间的差异
  # We calculate between-case variance using REML method (restricted maximum likelihood estimator) 使用REML方法计算案例间方差。 
  # tau^2：是案例间方差，认为不同研究之间有一些其他因素导致的差异。
  # I^2：去判断案例间差异大小占总差异的指标之一，但是I2不可以作为选择哪种模型（固定vs.随机）的依据。
  # Qt：效应值总体的异质性，是评价效应值的差异程度，表示效应值偏离均值的程度。
  # Qt越大，则效应值越离散，暗示我们有些因素对效应值有强烈的影响，我们可以去寻找一些因素，例如年龄性别，收集数据进行下一步分析。
  # Qt的优势是可以进行显著性检验的。如果p值不显著，那么我们认为案例间的差异是随机因素造成的，这种情况下就不需要往下继续进行Meta分析。
  # estimate 累积效应值，到底是大于0还是小于0，就能知道某种处理下，feature多了还是少了，是否显著，疾病下是否对feature影响明显呢。
  # ci.lb ci.ub 置信区间
  smd_rma <- rma(yi, vi, method = "REML", data = smd_meta)
  
  # merge each data
  smd_meta <- smd_rma$data %>% 
    add_column(taxa = taxa) %>% 
    as_tibble()
    
  smd_meta_model <- tibble(measure_method = "Standardized mean difference", # 效应值的计算方法
                           meta_model = "Random effect model", # 累积效应值计算模型
                           tau2_method = "REML", # 随机效应模型估计案例内方差（Tau^2）的计算方法
                           tau2_val = as.numeric(smd_rma$tau2), # 案例内方差的值
                           I2 = paste0(smd_rma$I2, "%"), # 案例间差异大小占总差异的比例
                           Qm = smd_rma$QM, # 某一因素引起的异质性
                           Qm_p = smd_rma$QMp, 
                           Qe = smd_rma$QE, # 残差的检验
                           Qe_p = smd_rma$QEp,
                           estimate = as.numeric(smd_rma$beta), # 模型累积效应值
                           ci_lb = smd_rma$ci.lb, # 累积效应值的上置信区间
                           ci_ub = smd_rma$ci.ub, # 累积效应值的下置信区间
                           pval = smd_rma$pval, # 显著性
                           taxa = taxa)
  
  smd_meta_dt <- rbind(smd_meta, smd_meta_dt)
  smd_meta_model_dt <- rbind(smd_meta_model, smd_meta_model_dt)
}

# 对接meta.metafor函数的结果 --------------
source("/share/data1/mjx/script/metafor.R")
meta <- meta.metafor(index, group = "group", group_pair = c("Disease", "Control"), proj = "study")
write.table(meta, "diversity_meta.txt", sep = "\t", quote = F, row.names = F)

# ggplot2
colors <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#FFD320","#ff7f00","#a65628","#f781bf","#999999")
ggplot() +
  geom_hline(yintercept = 0, color = "grey30", lwd = .5, lty = 2) +
  geom_point(data = meta, aes(x = feature, y = estimate), size = 3.5, shape = 15, inherit.aes = F) + 
  geom_errorbar(data = meta, aes(feature, ymin = ci_lb, ymax = ci_ub), width = .3, lwd = .5, inherit.aes = F) +
  geom_text(data = meta, aes(y = -1.22, x = feature, 
                             label = ifelse(pval < 0.05 | pval > 0.01, "*",
                                            ifelse(pval <= 0.01 | pval > 0.001, "**",
                                                   ifelse(p <= 0.001, "***","")))), color = "red") +
  geom_point(data = meta, aes(x = feature, y = yi, color = proj), size = 2.5, inherit.aes = F) +
  scale_color_manual(values = colors) +
  scale_y_continuous(limits = c(-1.28, 1.15), expand = c(0,0)) +
  labs(x = "", y = "Standradized Mean Difference (Random Effect Model)", color = "Study") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(color = "black"),
        axis.line.x = element_line(size = .5),
        axis.ticks.x = element_line(size = .5),
        aspect.ratio = .2)
ggsave("diversity_meta.pdf", width = 6, height = 3)
