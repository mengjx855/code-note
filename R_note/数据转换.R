####  clr transform ####
clr_transform <- function(data) {
  # 计算每个样本的几何平均数
  gmean <- apply(data, 2, \(x) exp(compositions::geometricmean(as.matrix(log(otu_table)))))
  
  # 计算每个观测值与其对应样本的几何平均数的比值
  ratio <- sweep(data, 2, gmean, "/")
  
  # 对比值取对数
  log_ratio <- log(ratio)
  
  return(log_ratio)
}

# 使用示例
sample1 = c(100, 200, 300, 400)
sample2 = c(200, 400, 500, 800)
sample3 = c(100, 500, 600, 1000)
otu_table = data.frame(sample1, sample2, sample3)
rownames(otu_table) = c("genus1", "genus2", "genus3")

clr_data = clr_transform(otu_table)
