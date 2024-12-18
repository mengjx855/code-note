####info####
# Encoding: utf-8
# modifier：Jinxin Meng
# Modified date：2023-1-5
# 理解adonis分析
# 参考：https://mp.weixin.qq.com/s?src=11&timestamp=1672907711&ver=4269&signature=9b5yeoMbHQOR*gwQoOMcmjS5SnIzAbcPAVsMB6UBopaCw5V9JmUqzYxev*8yv5N1jJbsXTo4pRvQ5*nDP1HqiqSy1RLfC*5zgdHWs2ysEYtiFzi1ptTHzJVzCTK1ccQv&new=1
# 参考：【你的adonis用对了吗？不同因素的顺序竟然对结果有很大影响】生信宝典

####PERMANOVA####
source("R2adjust.R")
group <- read.delim("../00.Data/sample_group.txt", sep = "\t")
load("../00.Data/profile.tpm.RData")
otu <- tpm[group$sample]

# distance
distance <- vegdist(t(otu), method = "bray")

group <- group[match(rownames(as.matrix(distance)), group$sample),] %>% 
  as.data.frame(row.names = NULL)  # 整体水平比较,保证样本在两个数据集中对应

# by = "margin"时，意思应该是控制其他因素的情况下，某个因素对菌群的影响。
# 举个例子，在控制除了年龄之外的因素干扰，年龄对于菌群时没有显著影响的。
adonis <- adonis2(distance ~ group + Age + BMI + SUA + SCr + Urea + CRP + ESR + eGFR, data = group, by = "margin", permutations = 999, parallel = 4)
# ermutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# adonis2(formula = distance ~ group + Age + BMI + SUA + SCr + Urea + CRP + ESR + eGFR, data = group, permutations = 999, by = "margin", parallel = 4)
#           Df SumOfSqs      R2      F Pr(>F)  
# group      1    0.484 0.01396 1.9930  0.029 *
# Age        1    0.175 0.00505 0.7216  0.765
# BMI        1    0.278 0.00803 1.1468  0.241
# SUA        1    0.195 0.00561 0.8013  0.669
# SCr        1    0.453 0.01308 1.8674  0.037 *
# Urea       1    0.280 0.00807 1.1528  0.263
# CRP        1    0.548 0.01582 2.2586  0.014 *
# ESR        1    0.213 0.00614 0.8773  0.532
# eGFR       1    0.179 0.00517 0.7375  0.736
# Residual 130   31.566 0.91056                
# Total    139   34.667 1.00000                
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# by = "terms"时，应该是根据顺序控制靠前的因素的情况下，下一个因素对菌群的影响。
# 举个例子，在控制除了疾病状态和年龄的干扰，BMI对于菌群时没有显著影响的；控制疾病状态，年龄，BMI，SUA的情况下，SCr这个临床因子对菌群有显著的影响。
adonis <- adonis2(distance ~ group + Age + BMI + SUA + SCr + Urea + CRP + ESR + eGFR, data = group, by = "margin", permutations = 999, parallel = 4)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(formula = distance ~ group + Age + BMI + SUA + SCr + Urea + CRP + ESR + eGFR, data = group, permutations = 999, by = "terms", parallel = 4)
#           Df SumOfSqs      R2      F Pr(>F)   
# group      1    0.592 0.01707 2.4367  0.011 * 
# Age        1    0.167 0.00483 0.6894  0.839   
# BMI        1    0.274 0.00791 1.1288  0.267   
# SUA        1    0.186 0.00536 0.7650  0.724   
# SCr        1    0.653 0.01884 2.6900  0.007 **
# Urea       1    0.217 0.00625 0.8923  0.511   
# CRP        1    0.620 0.01788 2.5528  0.013 * 
# ESR        1    0.213 0.00615 0.8775  0.534   
# eGFR       1    0.179 0.00517 0.7375  0.748   
# Residual 130   31.566 0.91056                 
# Total    139   34.667 1.00000                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# 分别计算每个影响因素对菌群的影响程度
# 写了个循环，之后考虑写个function
dat <- rbind()
for (i in c("group", "Age", "BMI", "SUA", "SCr", "Urea", "CRP", "ESR", "eGFR")) {
  adonis <- adonis2(formula = eval(parse(text = paste0("distance ~ ", i))), 
                    group, by = "margin", permutations = 999, parallel = 4)
  dat <- rbind(dat, data.frame(term = i, 
                               R2 = round(adonis[1,3], digits = 4),
                               R2adjust = get.adjusted.r2(adonis) %>% round(., digits = 4),
                               pval = adonis[1,5]))
}
write.table(dat, "adonis_for_single_term.txt", sep = "\t", quote = F, row.names = F)

# 可视化
ggplot(dat, aes(factor(term, levels = c("group", "Age", "BMI", "SUA", "SCr", "Urea", "CRP", "ESR", "eGFR")), R2*100)) +
  geom_bar(stat = "identity", fill = "grey60", width = .8, size = .5) +
  geom_text(aes(label = paste0("p = ", pval)), nudge_y = .05, color = "red", size = 2.5) +
  scale_x_discrete(expand = c(.06, .06)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "R2(%)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        aspect.ratio = .5)
ggsave("adonis_for_single_term.pdf", width = 5, height = 3)

