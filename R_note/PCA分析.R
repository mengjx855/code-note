# 如何比较PCA分析不同样本的差异？

# 1. PCA 得分比较
pca_result <- prcomp(profile_x[apply(profile_x, 1, \(x) sum(x)) != 0,] %>% t, scale. = T)
# 提取主成分得分
scores <- pca_result$x
# 添加组信息
scores_df <- scores %>% 
  data.frame %>% 
  rownames_to_column("sample") %>% 
  left_join(group_x %>% select(sample, group) %>% distinct(), by = "sample")
# 比较主成分得分
ggplot(scores_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  stat_ellipse()

# 2. 使用ANOVA比较主成分得分的组间差异
anova_result <- aov(PC1 ~ group, data = scores_df)
summary(anova_result)

# 3. 使用MANOVA比较多主成分得分
manova_result <- manova(cbind(PC1, PC2) ~ group, data = scores_df)
summary(manova_result)

# 4. 使用K均值聚类进行分析
set.seed(123)
kmeans_result <- kmeans(scores[, 1:2], centers = 4)
# 可视化聚类结果
ggplot(scores_df, aes(x = PC1, y = PC2, color = factor(kmeans_result$cluster))) +
  geom_point() +
  geom_text(aes(label = sample)) +
  labs(color = "")

# 5. 计算组中心之间的欧氏距离
library(cluster)
# 计算组中心
group_centers <- aggregate(scores[, 1:2], by = list(scores_df$group), FUN = mean)
# 计算欧氏距离
dist(group_centers %>% column_to_rownames("Group.1"))