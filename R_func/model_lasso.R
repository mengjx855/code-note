# Jinxin Meng, 20241023, 20241023 --------

# cv_fit <- cv.glmnet(profile_x, group_x$group, family = family, alpha = 1)
# 这里alpha=1为LASSO回归，如果等于0就是岭回归
# 参数 family 规定了回归模型的类型：
# family="gaussian" 适用于一维连续因变量（univariate）
# family="mgaussian" 适用于多维连续因变量（multivariate）
# family="poisson" 适用于非负次数因变量（count）
# family="binomial" 适用于二元离散因变量（binary）
# family="multinomial" 适用于多元离散因变量（category）
# 我们这里结局指标是2分类变量，所以使用binomial

# data <- predict(lasso_model, profile_y, type = 'class') 
# 在做分类预测的时候，如果预测值为概率，则type = "response" 给出具体的预测概率，而 type = "class" 按规定的阙值给出分类
# link: 返回线性预测器的值，即Xβ。
# response: 返回模型的预测概率（对于二分类问题），即应用了对数几率函数或逻辑回归函数的值。
# coefficients: 返回模型的系数值。
# class: 返回分类标签（仅适用于分类模型）。
# nonzero: 返回非零系数的索引

library(dplyr)
library(tibble)
library(tidyr)
library(glmnet)

# K折lasso建模 --------
# 输出一个预测数据框
# group文件是样本的分组信息，分组的列名为sample和group
# k为几折检验
lasso_Kfold <- function(profile, group, k, seed = 2024, family = "binomial"){
  start_time <- Sys.time()
  group <- group[match(colnames(profile), group$sample),]
  profile <- t(data.frame(profile, check.names = F))
  
  # 随机按照K折采样
  n_sample <- nrow(profile)
  size <- round(n_sample/k, 0)
  count <- seq_len(n_sample)
  sample_result <- list()
  for (i in 1:(k-1)) {
    set.seed(seed)
    sample_i <- sample(x = count, size = size, replace = F)
    sample_result[[i]] <- sample_i
    count <- setdiff(count, sample_i)
  }
  sample_result[[i+1]] <- count
  
  # 训练和测试
  pred <- rbind()
  pb <-  txtProgressBar(style = 3, width = 50, char = "#")
  for(i in 1:length(sample_result)){ # 划分测试集和训练集
    sample_i <- sample_result[[i]]
    test_data <- profile[sample_i, ]
    test_group <- group[sample_i, ]
    train_data <- profile[-sample_i, ]
    train_group <- group[-sample_i, ]
    
    set.seed(seed)
    cv_fit <- cv.glmnet(train_data, train_group$group, family = family, alpha = 1)
    lasso_model <- glmnet(train_data, train_group$group, family = family, alpha = 1, lambda = cv_fit$lambda.min)
    pred_i <- data.frame(predict(lasso_model, test_data, type = "response"))
    pred <- rbind(pred, pred_i)
    setTxtProgressBar(pb, i/length(sample_result)) # 进度计算
  }
  close(pb)
  
  pred <- dplyr::rename(pred, pred = 1) %>%
    rownames_to_column("sample")
  end_time <- Sys.time()
  run_time <- end_time - start_time
  cat(paste0(" Consumption of time: ", round(as.numeric(run_time), 3), " sec\n"))
  return(pred)
}

# profile_x 建模，profile_y 验证 -----------------
# 两个otu，第一个是发现集，第二个为验证集
# group文件是样本的分组信息，分组的列至少有sample和group,第一列必须为"sample"，第二列必须为"group"
# 返回一个预测结果的数据框
lasso_next_vaildate <- function(profile_x, profile_y, group_x, group_y, seed = 2024, 
                                family = "binomial",
                                label = "otu_x for modeling and otu_y for predicting"){
  group_x <- group_x[match(colnames(profile_x), group_x$sample),]
  profile_y <- profile_y[rownames(profile_x),]
  
  # discovery
  profile_x <- t(profile_x)
  cv_fit <- cv.glmnet(profile_x, group_x$group, family = family, alpha = 1)
  # 这里alpha=1为LASSO回归，如果等于0就是岭回归
  # 参数 family 规定了回归模型的类型：
  # family="gaussian" 适用于一维连续因变量（univariate）
  # family="mgaussian" 适用于多维连续因变量（multivariate）
  # family="poisson" 适用于非负次数因变量（count）
  # family="binomial" 适用于二元离散因变量（binary）
  # family="multinomial" 适用于多元离散因变量（category）
  # 我们这里结局指标是2分类变量，所以使用binomial
  best_lambda <- cv_fit$lambda.min
  
  lasso_model <- glmnet(profile_x, group_x$group, family = family, alpha = 1, lambda = best_lambda)

  # validation
  group_y <- group_y[match(colnames(profile_y), group_y$sample),]
  profile_y <- t(profile_y)
  
  pred <- predict(lasso_model, profile_y, type = 'response') %>% 
    data.frame() %>% 
    dplyr::rename(pred = 1) %>% 
    rownames_to_column(var = "sample")
  
  data <- predict(lasso_model, profile_y, type = 'class') %>% 
    data.frame() %>% 
    dplyr::rename(pred = 1) %>% 
    rownames_to_column(var = "sample") %>% 
    left_join(group_y, by = "sample")
  
  # 在做分类预测的时候，如果预测值为概率，则type = "response" 给出具体的预测概率，而 type = "class" 按规定的阙值给出分类
  # link: 返回线性预测器的值，即Xβ。
  # response: 返回模型的预测概率（对于二分类问题），即应用了对数几率函数或逻辑回归函数的值。
  # coefficients: 返回模型的系数值。
  # class: 返回分类标签（仅适用于分类模型）。
  # nonzero: 返回非零系数的索引
  
  # model_evaluation
  confusion_matrix <- table(pred = data$pred, actu = data$group)
  
  # 准确率
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  # 灵敏度和特异性
  sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
  specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[1, ])

  out <- list(model = lasso_model,
              pred = pred,
              method = "lasso",
              confusion_matrix = confusion_matrix,
              accuracy = accuracy,
              sensitivity = sensitivity, 
              specificity = specificity)
  
  return(out)
}