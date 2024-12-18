# Jinxin Meng, 20241023, 20241023 ---------------------

# https://blog.csdn.net/dege857/article/details/111693504
# LASSO 回归也叫套索回归，是通过生成一个惩罚函数是回归模型中的变量系数进行压缩，达到防止过度拟合，解决严重共线性的问题，LASSO 回归最先由英国人Robert Tibshirani提出，目前在预测模型中应用非常广泛
# 在新格兰文献中，有大牛提出，对于变量过多而且变量数较少的模型拟合，首先要考虑使用LASSO 惩罚函数。今天我们来讲讲怎么使用R语言通过LASSO 回归构造预测模型。
# 首先我们要下载R的glmnet包，由 LASSO 回归的发明人，斯坦福统计学家 Trevor Hastie 领衔开发。

# LASSO，全称Least absolute shrinkage and selection operator，是一种数据挖掘方法，即在常用的多元线性回归中，添加惩罚函数，不断压缩系数，从而达到精简模型的目的，以避免共线性和过拟合。
# 当系数为0时，同时达到筛选变量的效果。（以下是一个不严谨的示意图）
# Y = β0 + β1X1 + β2X2 + β3X3 + β4X4 + βnXn
#          |
# 引入惩罚分数，压缩系数，当β1=0时，X1被剔除方程中
#          |
# Y = β0 + ———— + β2X2 + β3X3 + β4X4 + βnXn

# 所以，LASSO回归高效解决了筛选变量的难题：区别于传统的逐步回归stepwise前进、后退变量筛选方法，LASSO回归可以利用较少样本量，高效筛选较多变量。
# 比如在基因组学、影像学、以及其他小样本分析中，LASSO回归都可以派上大用场。

# 1. 安装并加载必要的包
# 首先，确保你已经安装了 glmnet 包。
install.packages("glmnet")
library(glmnet)

# 2. 数据准备
# 假设你已经有两个数据框 microbiome_data 和 labels，其中 microbiome_data 包含样本的特征，labels 包含相应的标签（病人或健康人）。
# 示例数据
set.seed(123)
microbiome_data <- matrix(rnorm(200 * 100), ncol = 100)  # 200 样本，100 特征
labels <- sample(c(0, 1), 200, replace = TRUE)  # 0 表示健康，1 表示病人

# 3. 数据标准化
# LASSO 回归对特征的尺度敏感，建议进行标准化。
normalized_data <- scale(microbiome_data)

# 4. 拆分训练集和测试集
# 将数据分成训练集和测试集。
set.seed(456)
train_index <- sample(1:nrow(normalized_data), 0.7 * nrow(normalized_data))  # 70% 训练集
train_data <- normalized_data[train_index, ]
test_data <- normalized_data[-train_index, ]
train_labels <- labels[train_index]
test_labels <- labels[-train_index]

# 5. 模型训练
# 使用交叉验证来选择最优的惩罚参数 lambda。
cv_fit <- cv.glmnet(train_data, train_labels, family = "binomial", alpha = 1)
#这里alpha=1为LASSO回归，如果等于0就是岭回归
#参数 family 规定了回归模型的类型：
#family="gaussian" 适用于一维连续因变量（univariate）
#family="mgaussian" 适用于多维连续因变量（multivariate）
#family="poisson" 适用于非负次数因变量（count）
#family="binomial" 适用于二元离散因变量（binary）
#family="multinomial" 适用于多元离散因变量（category）
#我们这里结局指标是2分类变量，所以使用binomial

# 获取最优的 lambda 值
best_lambda <- cv_fit$lambda.min
# 使用最优 lambda 重新训练模型
lasso_model <- glmnet(train_data, train_labels, family = "binomial", alpha = 1, lambda = best_lambda)

# 6. 模型预测
# 使用训练好的模型进行预测。
# 对测试集进行预测
predicted_probabilities <- predict(lasso_model, newx = test_data, type = "response")
# 在做分类预测的时候，如果预测值为概率，则type = "response" 给出具体的预测概率，而 type = "class" 按规定的阙值给出分类
# link: 返回线性预测器的值，即Xβ。
# response: 返回模型的预测概率（对于二分类问题），即应用了对数几率函数或逻辑回归函数的值。
# coefficients: 返回模型的系数值。
# class: 返回分类标签（仅适用于分类模型）。
# nonzero: 返回非零系数的索引
predicted_labels <- ifelse(predicted_probabilities > 0.5, 1, 0)

# 7. 评估模型性能
# 使用混淆矩阵、准确率、灵敏度、特异性等指标评估模型性能。
# 混淆矩阵
confusion_matrix <- table(Predicted = predicted_labels, Actual = test_labels)
print(confusion_matrix)
# 准确率
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(paste("Accuracy:", accuracy))
# 灵敏度和特异性
sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[1, ])
print(paste("Sensitivity:", sensitivity))
print(paste("Specificity:", specificity))
