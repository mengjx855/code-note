# Jinxin Meng, 20241024, 20241024 --------

# 核函数, 中文, 适用范围, 参数
# linear, 线性核, 线性数据, 无
# polynomial, 多项式核, 偏向于线性数据, gamma/degree/coef0
# radial basis, RBF或高斯径向基核, 偏向于非线性数据, gamma
# sigmoid, sigmoid核, 非线性数据, gamma/coef0

# 使用默认的tune.svm()调整超参数，我们就用常见的径向基核为例进行演示，我们同时调整2个超参数：gamma和cost。
# tune_model <- tune.svm(Species ~ ., data = iris, cost = 10^(-1:3), # 设置cost的值
#                        gamma = 10^(-3:1), # 设置gamma的值
#                        # 重抽样方法选择自助法，次数选择100次
#                        tunecontrol=tune.control(sampling = "bootstrap", nboot = 100 )) #自助法
# Bootstrap 自助法 放回抽样生成多个样本集
## - best parameters:
##  gamma cost
##   0.01  100

library(dplyr)
library(tibble)
library(tidyr)
library(e1071)
library(pROC)
source("/code/R_func/plot_roc.R")

# Holdout留出法建模 ---------
svm_base <- function(profile, group, rep = 5, seed = 2024, train_perc = .8, kernel = "radial", 
                     scale = T, probability = T) {
  start_time <- Sys.time()
  profile <- t(profile)
  if(rep <= 1) seqs = 1
  if(rep > 1) seqs = seq_len(rep)
  
  auc = 0
  pb <- txtProgressBar(style = 3, width = 50, char = "#")
  for (i in seqs) {
    set.seed(seed + i - 1)
    sample_i <- caret::createDataPartition(group$group, p = train_perc)[[1]]
    train_group <- group_x[sample_i, ]
    train_data <- profile_x[train_group$sample, ]
    test_group <- group_x[-sample_i, ]
    test_data <- profile_x[test_group$sample, ]
    
    tune_model <- tune.svm(train_data, factor(train_group$group), cost = 10^(-1:3), gamma = 10^(-3:1),
                           tunecontrol = tune.control(sampling = "bootstrap", nboot = 10))
    svm_model <- svm(train_data, factor(train_group$group), kernel = kernel, scale = scale, probability = probability, 
                     gamma = tune_model$best.parameters[1, "gamma"], cost = tune_model$best.parameters[1, "cost"])
    
    pred <- predict(svm_model, test_data, probability = T)
    pred <- attr(pred, "probabilities") %>% 
      data.frame() %>% 
      rownames_to_column("sample") %>% 
      left_join(group, by = "sample")
    roc <- roc(pred$group, pred$SLE, quiet = T)
    
    if(roc$auc > auc) {
      auc <- roc$auc
      out = list(
        seed = seed + i - 1,
        train_data = train_data,
        train_group = train_group,
        test_data = test_data,
        test_group = test_group,
        tune_model = tune_model,
        svm_model = svm_model,
        kernel = kernel,
        pred = pred,
        roc = roc,
        roc_plot = plot_roc(roc) )
    }
    
    setTxtProgressBar(pb, i/rep) # 进度计算
  }
  close(pb)
  end_time <- Sys.time()
  run_time <- difftime(end_time, start_time, units = "sec")
  cat(paste0(" Consumption of time: ", round(as.numeric(run_time), 3), " sec\n"))
  return(out)
}

# K折svm建模 --------
# 输出一个预测数据框
# group文件是样本的分组信息，分组的列名为sample和group
# k为几折检验
svm_Kfold <- function(profile, group, k = 5, seed = 2024, kernel = "radial",
                      scale = T, probability = T){
  start_time <- Sys.time()
  group <- group[match(colnames(profile), group$sample),]
  profile <- t(data.frame(profile, check.names = F))
  
  sample_result <- caret::createFolds(group$sample, k = k)
  
  pred <- rbind()
  pb <-  txtProgressBar(style = 3, width = 50, char = "#")
  for(i in 1: k){ # 划分测试集和训练集
    sample_i <- sample_result[[i]]
    test_data <- profile[sample_i, ]
    test_group <- group[sample_i, ]
    train_data <- profile[-sample_i, ]
    train_group <- group[-sample_i, ]
    
    set.seed(seed)
    tune_model <- tune.svm(train_data, factor(train_group$group), cost = 10^(-1:3), gamma = 10^(-3:1),
                           tunecontrol = tune.control(sampling = "bootstrap", nboot = 10))
    svm_model <- svm(train_data, factor(train_group$group), kernel = kernel, scale = scale, probability = probability, 
                     gamma = tune_model$best.parameters[1, "gamma"], cost = tune_model$best.parameters[1, "cost"])
    
    pred_i <- predict(svm_model, test_data, probability = T)
    pred_i <- attr(pred_i, "probabilities") %>% 
      data.frame() %>% 
      rownames_to_column("sample")
    pred <- rbind(pred, pred_i)
    setTxtProgressBar(pb, i / k) # 进度计算
  }
  close(pb)

  end_time <- Sys.time()
  run_time <- difftime(end_time, start_time, units = "sec")
  cat(paste0(" Consumption of time: ", round(as.numeric(run_time), 3), " sec\n"))
  return(pred)
}

# profile_x 建模，profile_y 验证 -----------------
# 两个otu，第一个是发现集，第二个为验证集
# group文件是样本的分组信息，分组的列至少有sample和group,第一列必须为"sample"，第二列必须为"group"
# 返回一个预测结果的数据框
svm_next_vaildate <- function(profile_x, profile_y, group_x, group_y, seed = 2024, 
                              kernel = "radial", scale = T, probability = T,
                              label = "otu_x for modeling and otu_y for predicting"){
  start_time <- Sys.time()
  
  group_x <- group_x[match(colnames(profile_x), group_x$sample),]
  profile_y <- profile_y[rownames(profile_x),]
  
  # discovery
  profile_x <- t(profile_x)
  set.seed(seed)
  tune_model <- tune.svm(profile_x, factor(group_x$group), cost = 10^(-1:3), gamma = 10^(-3:1),
                         tunecontrol = tune.control(sampling = "bootstrap", nboot = 20))
  svm_model <- svm(profile_x, factor(group_x$group), kernel = kernel, scale = scale, probability = probability, 
                   gamma = tune_model$best.parameters[1, "gamma"], cost = tune_model$best.parameters[1, "cost"])
  
  # validation
  group_y <- group_y[match(colnames(profile_y), group_y$sample),]
  profile_y <- t(profile_y)
  
  pred <- predict(svm_model, profile_y, probability = T) %>% 
    attr("probabilities") %>% 
    data.frame() %>% 
    rownames_to_column(var = "sample") %>% 
    left_join(group_y, by = "sample")
  roc <- roc(pred$group, pred[,2])
  
  # model_evaluation
  confusion_matrix <- table(pred = predict(svm_model, profile_y), actu = pred$group)
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix) # 准确率
  sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) # 灵敏度和特异性
  specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[1, ])

  out <- list(tune_model = tune_model,
              model = svm_model,
              pred = pred,
              roc = roc,
              roc_plot = plot_roc(roc),
              method = "svm",
              confusion_matrix = confusion_matrix,
              accuracy = accuracy,
              sensitivity = sensitivity, 
              specificity = specificity)
  
  end_time <- Sys.time()
  run_time <- difftime(end_time, start_time, units = "sec")
  cat(paste0(" Consumption of time: ", round(as.numeric(run_time), 3), " sec\n"))
  return(out)
}