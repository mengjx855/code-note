# Jinxin Meng, 20241024, 20241024 ---------------------

# https://blog.csdn.net/qq_53123067/article/details/136060974
# 支持向量机（Support Vector Machine，SVM）是一种经典的监督学习算法，用于解决二分类和多分类问题。其核心思想是通过在特征空间中找到一个最优的超平面来进行分类，并且间隔最大。
# SVM能够执行线性或非线性分类、回归，甚至是异常值检测任务。它是机器学习领域最受欢迎的模型之一。SVM特别适用于中小型复杂数据集的分类。

# 超平面最大间隔介绍：决策边界，不仅分离了两个类别，且尽可能远离最近的训练实例。
# 硬间隔分类：在上面我们使用超平面进行分割数据的过程中，如果我们严格地让所有实例都不在最大间隔之间，并且位于正确的一边，这就是硬间隔分类。硬间隔分类有两个问题，首先，它只在数据是线性可分离的时候才有效；其次，它对异常值非常敏感。当有一个额外异常值的鸢尾花数据：左图的数据根本找不出硬间隔，而右图最终显示的决策边界与我们之前所看到的无异常值时的决策边界也大不相同，可能无法很好地泛化。
# 软间隔分类：要避免这些问题，最好使用更灵活的模型。目标是尽可能在保持最大间隔宽阔和限制间隔违例（即位于最大间隔之上，甚至在错误的一边的实例）之间找到良好的平衡，这就是软间隔分类。

# 在Scikit-Learn的SVM类中，可以通过超参数C来控制这个平衡：C值越小，则间隔越宽，但是间隔违例也会越多。上图显示了在一个非线性可分离数据集上，两个软间隔SVM分类器各自的决策边界和间隔。
# 左边使用了高C值，分类器的错误样本（间隔违例）较少，但是间隔也较小。
# 右边使用了低C值，间隔大了很多，但是位于间隔上的实例也更多。看起来第二个分类器的泛化效果更好，因为大多数间隔违例实际上都位于决策边界正确的一边，所以即便是在该训练集上，它做出的错误预测也会更少。

# SVM算法的优点：
# 1）SVM方法既可以用于分类（二/多分类），也可用于回归和异常值检测。
# 2）SVM具有良好的鲁棒性，对未知数据拥有很强的泛化能力，特别是在数据量较少的情况下，相较其他传统机器学习算法具有更优的性能。
# 使用SVM作为模型时，通常采用如下流程：
# 1）对样本数据进行归一化
# 2）应用核函数对样本进行映射（最常采用和核函数是RBF和Linear，在样本线性可分时，Linear效果要比RBF好)
# 3）用cross-validation和grid-search对超参数进行优选
# 4）用最优参数调练得到模型
# 5）测试

# https://ayueme.github.io/machine_learning_base_r/%E6%94%AF%E6%8C%81%E5%90%91%E9%87%8F%E6%9C%BA.html
# 支持向量机的原理非常复杂，我这里尽量用简单的话说明。还是用前一章中的二分类问题为例进行说明。
# 下图有两个特征可以用来预测肿瘤是”良性”还是”恶性”，SVM的关键就是找到一条完美的线（或者平面）把它分成两类，下图中展示的实线和虚线都可以做到这一点，但这只是无数条线其中的两条而已。
# 如果数据的维度多于2，那就不是一条线能解决的了，那就需要一个平面，这个平面就被称为超平面，
# 如果一个超平面能够把数据分为两类，每类只包含一种类别，那么这个超平面就是数据的决策边界。比如上图中的实线和虚线就是两个决策边界。
# 如上图所示，我们可以找到多个决策边界把数据分为两类，但是哪一个才是最好的呢？
# 如下左图所示，我们把其中一条决策边界（比如实线）往左右两边移动，当碰到第一个数据点时停下来，此时两条实线之间的距离就被称为（中间）实线这条决策边界的边际（margin）。

# 支持向量机就是要找到使得边际最大的决策边界来实现对数据的分类。如上右图所示的边际就明显小于作图所示的边际。
# 如果这条决策边界非常完美，能够把所有的样本完美分成两类，那么此时模型的误差就是0，因为完全分类正确，但是此时模型对新数据的预测可能会出现较大的误差，也就是说出现了过拟合。
# 上图中在左右边界上的样本被称为支持向量（support vector）。
# 实际生活中很多数据并不是完全线性可分的（上面这个例子是线性可分的），为了能够适用于非线性的数据，支持向量机可以将原始数据进行各种转换，这就是核技巧（kernel trick），或者叫核函数（kernel function）。通过使用的不同的核技巧，支持向量机可以适用于多种不同类型的数据。

# SVM有多种实现方法，最著名的就是台湾林智仁教授开发的LIBSVM了，python中的sklearn中的支持向量机就是使用的这种方法，在R语言中是通过e1071这个包实现的。
# 该包还可以实现缺失值插补、异常值检测、计算偏度/峰度、多种算法的超参数调优（决策树、随机森林、knn、神经网络等）等，功能非常强大。

library(e1071)
#example("svm")

# 参数解释
# 支持向量机在e1071中主要是通过svm()这个函数实现的。既可以使用公式语法，也可以使用x/y的形式，此时如果因变量y是因子型，则函数默认这是一个分类任务，如果y不是因子型则默认为回归任务，如果不提供y，则默认是异常值检测。
svm(x, y = NULL, scale = TRUE, type = NULL, kernel = "radial", degree = 3, 
    gamma = if (is.vector(x)) 1 else 1 / ncol(x),
    coef0 = 0, cost = 1, nu = 0.5, class.weights = NULL, cachesize = 40, 
    tolerance = 0.001, epsilon = 0.1, shrinking = TRUE, cross = 0, 
    probability = FALSE, fitted = TRUE, ..., subset, na.action = na.omit)

# 部分参数介绍：
# scale：是否进行标准化，默认是TRUE
# type：任务类型，比如回归还是分类，具体可选以下几种：C-classification，nu-classification，one-classification (for novelty detection)，eps-regression，nu-regression
# kernel：e1071中的核函数有4种，每个核函数都有最适合的情况，但通常都会选择径向基核函数。除了线性核之外，其余3种都有属于自己的超参数：
# linear：线性核函数
# polynomial：多项式核函数，(3个超参数:gamma, degree, coef0)
# radial basis：RBF，径向基核函数，(1个超参数: gamma)
# sigmoid： (2个超参数:gamma, coef0)
# degree：多项式核的阶数，默认是3，如果是1就类似于线性核
# gamma：默认是1/特征数量，取值范围是0到正无穷，一般在优化这个参数时会进行log转换，gamma越大，通常导致支持向量越多
# coef0：多项式核的参数，默认是0，需要整数
# cost：正则化程度，或者表示犯错的成本，类似于sklearn中的C，正整数。一个较大的成本意味着模型对误差的惩罚更大，从而将生成一个更复杂的分类边界，对应的训练集中的误差也会更小（也就是边际大），但也意味着可能存在过拟合问题，即对新样本单元的预测误差可能很大。较小的成本意味着分类边界更平滑，但可能会导致欠拟合
# class.weight：因变量的权重，比如svm(x, y, class.weights = c(A = 0.3, B = 0.7))，指定inverse会直接使用和类别比例相反的类别权重
# cachesize：使用的内存大小，默认40MB
# cross：交叉验证的折数，默认是0，不进行交叉验证。如果使用了交叉验证，结果中会多出一些额外信息，回归任务会多出：
# accuracies：每一折的准确率；
# tot.accuracy:平均准确率。 分类任务会多出：
# MSE:每一折的MSE；
# tot.MSE：平均MSE；
# scorrcoef:相关系数的平方
# probability:是否生成各个类别的概率
# fitted：是否生成分类结果

# 这个核函数的概念非常重要，所以我做了一个表格：
# 核函数           中文             适用范围            参数
# linear         线性核             线性数据            无
# polynomial    多项式核            偏向于线性数据      gamma/degree/coef0
# radial basis  RBF,高斯径向基核    偏向于非线性数据      gamma
# sigmoid       sigmoid核           非线性数据          gamma/coef0

# amma/degree/coef0这3个参数对SVM的影响非常复杂，因为并不是完全的正相关或者负相关关系，尤其是对多项式核、sigmoid核来说，多个参数会共同影响，关系就更加复杂了，所以通常此时我们会使用网格搜索等方法确定这些参数的值。对于RBF核来说，它只有一个参数gamma，此时我们可以用学习曲线来确定最佳的gamma值。

# 建立模型
# 演示数据就用R语言自带的鸢尾花数据集，这个数据集是一个3分类的数据，其中Species是结果变量，为3分类，其余变量是预测变量。
rm(list = ls())
data(iris)

# 公式形式
model <- svm(Species ~ ., data = iris, probability = T)

# 或者x/y形式
x <- subset(iris, select = -Species)
y <- iris$Species
#model <- svm(x, y, probability = T) 

print(model)
## 
## Call:
## svm(formula = Species ~ ., data = iris, probability = T)
## 
## 
## Parameters:
##    SVM-Type:  C-classification 
##  SVM-Kernel:  radial 
##        cost:  1 
## 
## Number of Support Vectors:  51
summary(model)
## 
## Call:
## svm(formula = Species ~ ., data = iris, probability = T)
## 
## 
## Parameters:
##    SVM-Type:  C-classification 
##  SVM-Kernel:  radial 
##        cost:  1 
## 
## Number of Support Vectors:  51
## 
##  ( 8 22 21 )
## 
## 
## Number of Classes:  3 
## 
## Levels: 
##  setosa versicolor virginica

# 查看结果
# 查看支持向量：
# 51个支持向量，这是预处理(比如na.omit)后的序号
model$index
##  [1]   9  14  16  21  23  24  26  42  51  53  54  55  57  58  60  61  64  67  69
## [20]  71  73  77  78  79  84  85  86  87  88  99 107 109 111 117 119 120 122 124
## [39] 126 127 128 130 132 134 135 138 139 143 147 149 150

# 查看预测类别和预测概率，只有在建模时使用了probability=T，才能在predict()时使用此参数：

# 获取训练集预测结果
# newdata写测试集就是测试集结果
pred <- predict(model, newdata = x, probability = T)

# 查看结果
head(pred)
##      1      2      3      4      5      6 
## setosa setosa setosa setosa setosa setosa 
## Levels: setosa versicolor virginica
# 查看预测概率
head(attr(pred,"probabilities"))
##      setosa versicolor   virginica
## 1 0.9816955 0.01075790 0.007546567
## 2 0.9746218 0.01725092 0.008127315
## 3 0.9804885 0.01135611 0.008155367
## 4 0.9766379 0.01459413 0.008767961
## 5 0.9809263 0.01108326 0.007990473
## 6 0.9756922 0.01600990 0.008297935
# 查看混淆矩阵：
# 混淆矩阵
table(pred, y)
##             y
## pred         setosa versicolor virginica
##   setosa         50          0         0
##   versicolor      0         48         2
##   virginica       0          2        48

# 或者caret版混淆矩阵
#caret::confusionMatrix(table(pred, y))
caret::confusionMatrix(pred, y)
## Confusion Matrix and Statistics
## 
##             Reference
## Prediction   setosa versicolor virginica
##   setosa         50          0         0
##   versicolor      0         48         2
##   virginica       0          2        48
## 
## Overall Statistics
##                                           
##                Accuracy : 0.9733          
##                  95% CI : (0.9331, 0.9927)
##     No Information Rate : 0.3333          
##     P-Value [Acc > NIR] : < 2.2e-16       
##                                           
##                   Kappa : 0.96            
##                                           
##  Mcnemar's Test P-Value : NA              
## 
## Statistics by Class:
## 
##                      Class: setosa Class: versicolor Class: virginica
## Sensitivity                 1.0000            0.9600           0.9600
## Specificity                 1.0000            0.9800           0.9800
## Pos Pred Value              1.0000            0.9600           0.9600
## Neg Pred Value              1.0000            0.9800           0.9800
## Prevalence                  0.3333            0.3333           0.3333
## Detection Rate              0.3333            0.3200           0.3200
## Detection Prevalence        0.3333            0.3333           0.3333
## Balanced Accuracy           1.0000            0.9700           0.9700

# 三分类也是可以绘制ROC曲线的，只不过稍微复杂一点，大家有需要的可以参考ROC曲线绘制合集。
# https://mp.weixin.qq.com/mp/appmsgalbum?__biz=MzUzOTQzNzU0NA==&action=getalbum&album_id=3157190853568921605&uin=&key=&devicetype=Windows+11+x64&version=63090819&lang=zh_CN&ascene=0

# 可视化决策边界
# 如果预测变量超过1个就要使用formula，因为只能在两个维度上画出这个图，slice=list(Sepal.Width=3,Sepal.Length=4)表示要把Sepal.Width这一列都当成3，把Sepal.Length这一列都变成4，这个数值具体怎么影响这幅图也没找到合适的解释，希望有大佬能指点迷津。
plot(model, data = iris, formula = Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4),
     svSymbol = "x" # 支持向量的形状
     )
# 图中的x表示支持向量，o表示数据点。

# 不同核函数比较
# 我们把4种核函数的结果都画出来，放在一起比较下决策边界。
model_rbf <- svm(Species ~ ., data = iris, kernel = "radial")#径向基核
model_linear <- svm(Species ~ ., data = iris, kernel = "linear")#线性核
model_ploy <- svm(Species ~ ., data = iris, kernel = "polynomial")#多项式核
model_sig <- svm(Species ~ ., data = iris, kernel = "sigmoid")#sigmoid核

par(mfrow=c(2,2))
plot(model_rbf, data = iris, formula = Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))
plot(model_linear, data = iris, formula = Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))
plot(model_ploy, data = iris, formula = Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))
plot(model_sig, data = iris, formula = Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))

# 超参数调优
# 使用默认的tune.svm()调整超参数，我们就用常见的径向基核为例进行演示，我们同时调整2个超参数：gamma和cost。
set.seed(123)
tune_model <- tune.svm(Species ~ ., data = iris,
                       cost = 10^(-1:3), # 设置cost的值
                       gamma = 10^(-3:1), # 设置gamma的值
                       
                       # 重抽样方法选择自助法，次数选择100次
                       tunecontrol=tune.control(sampling = "bootstrap",
                                                nboot = 100
                                                )#自助法
                       )

tune_model
## 
## Parameter tuning of 'svm':
## 
## - sampling method: bootstrapping 
## 
## - best parameters:
##  gamma cost
##   0.01  100
## 
## - best performance: 0.03924088

# 错误率最低是0.03924088，此时gamma=0.01，cost=100。
# 使用这个超参数重新拟合模型：

model_final <- svm(Species ~ ., data = iris, cost=100, gamma=0.01)
print(model_final)
## 
## Call:
## svm(formula = Species ~ ., data = iris, cost = 100, gamma = 0.01)
## 
## 
## Parameters:
##    SVM-Type:  C-classification 
##  SVM-Kernel:  radial 
##        cost:  100 
## 
## Number of Support Vectors:  24

# 预测新数据
test_df <- head(iris)
pred <- predict(model, newdata = test_df)
caret::confusionMatrix(pred, test_df$Species)
## Confusion Matrix and Statistics
## 
##             Reference
## Prediction   setosa versicolor virginica
##   setosa          6          0         0
##   versicolor      0          0         0
##   virginica       0          0         0
## 
## Overall Statistics
##                                      
##                Accuracy : 1          
##                  95% CI : (0.5407, 1)
##     No Information Rate : 1          
##     P-Value [Acc > NIR] : 1          
##                                      
##                   Kappa : NaN        
##                                      
##  Mcnemar's Test P-Value : NA         
## 
## Statistics by Class:
## 
##                      Class: setosa Class: versicolor Class: virginica
## Sensitivity                      1                NA               NA
## Specificity                     NA                 1                1
## Pos Pred Value                  NA                NA               NA
## Neg Pred Value                  NA                NA               NA
## Prevalence                       1                 0                0
## Detection Rate                   1                 0                0
## Detection Prevalence             1                 0                0
## Balanced Accuracy               NA                NA               NA

# 可视化决策边界：
plot(model_final, data = iris, formula = Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))

# 使用建议
# SVM对超参数值很敏感，建议多试几次；
# 对于分类任务，首选C-classification和RBF(default)，表现好且超参数比较少，只有cost和gamma。libsvm的作者建议先对cost用一个很大的范围，比如1~1000，然后用交叉验证选择比较好的，最后用这几个选中的cost值再尝试不同的gamma；
# 可以使用网格搜索实现超参数调优，tune.svm()；
# 对数据进行标准化会提高模型表现，强烈建议，svm()默认会对数据进行标准化。

# 除了支持向量机外，e1071包还支持对其他模型进行超参数调优，包括：
# tune.gknn()：generalized k-nearest neighbors
# tune.knn()：knn
# tune.nnet()：nnet，神经网络
# tune.randomForest()：随机森林
# tune.rpart()：决策树

# 这几个函数都是tune()的变体。每种算法可以调优的超参数都在帮助文档中有详细的说明，作为一个轻量化的调优R包来用是非常不错的选择，虽然我们有更好的方法。
# 支持向量机的递归特征消除，可参考推文：递归特征消除

# 参考资料
# example(“svm”)
# https://stackoverflow.com/questions/49877752/slice-option-in-plot-function-of-package-e1071
# https://c3h3notes.wordpress.com/2010/10/20/r%E4%B8%8A%E7%9A%84libsvm-package-e1071
# https://www.datacamp.com/tutorial/support-vector-machines-r
# https://rpubs.com/cliex159/865583
# 菜菜的sklearn课堂

# 下面给大家演示下R语言做支持向量机的例子，并且比较下在不进行调参的默认情况下，4种核函数的表现情况。分别是：线性核，多项式核，高斯径向基核，sigmoid核。
# 支持向量机非常强，应用非产广泛，不管是分类还是回归都能用，万金油一样的算法。不过它的理论知识比随机森林复杂了非常多，但是实现起来并不难哈，我们就直接调包即可。
# 加载数据和R包
# 使用e1071包做演示。数据使用印第安人糖尿病数据集。
library(e1071)
library(pROC)
library(dplyr)
rm(list = ls())
load(file = "datasets/pimadiabetes.rdata")

# 数据划分
# 划分训练集和测试集，经典7：3分：
# 划分是随机的，设置种子数可以让结果复现
set.seed(123)
ind <- sample(1:nrow(pimadiabetes), size = 0.7*nrow(pimadiabetes))
# 训练集、测试集
train <- pimadiabetes[ind,]
test <- pimadiabetes[-ind, ]

# 训练集建模
# e1071使用起来非常简单，直接一个函数搞定，也是使用R语言经典的formula写法，二分类数据我们通常希望获得预测概率，所以加上probability = TRUE
# 然后kernel参数就是分别用4种核函数。
set.seed(123)
svmLinear <- svm(diabetes ~ ., data = train,
                 probability = TRUE,
                 kernel="linear"
                 )

svmPoly <- svm(diabetes ~ ., data = train,
               probability = TRUE,
               kernel="polynomial"
               )

svmRadial <- svm(diabetes ~ ., data = train,
                 probability = TRUE,
                 kernel="radial"
                 )

svmSigmoid <- svm(diabetes ~., data = train,
                  probability = TRUE,
                  kernel="sigmoid"
                  )

# 接下来就是查看模型在训练集中的表现，我们为了少写几行代码，先定义一个函数，可以自定帮我们提取训练结果，并组成1个数据框，内含原始数据的结果变量，预测结果，预测概率。

# 定义函数
getres <- function(svmfunc, dataset){
  data_pred <- predict(svmfunc, newdata=dataset, probability = T)
  data_pred_df <- dataset %>% dplyr::select(diabetes) %>% 
  dplyr::bind_cols(status_pred = data_pred) %>% 
  dplyr::bind_cols(attr(data_pred, "probabilities"))
}

# 接下来提取数据即可，我们先提取1个看看：
Linear_train_pred_df <- getres(svmLinear, train)
head(Linear_train_pred_df)
##     diabetes status_pred        neg       pos
## 415      neg         pos 0.36393796 0.6360620
## 463      pos         pos 0.16223183 0.8377682
## 179      pos         neg 0.73073668 0.2692633
## 526      pos         pos 0.04261936 0.9573806
## 195      pos         pos 0.08214236 0.9178576
## 118      pos         pos 0.12508113 0.8749189

# 上面这个是：线性核函数，训练集，的结果，看起来没什么问题，第一列是真实结果变量，第2列是预测结果类别，第3和4列是预测的类别概率。
# 如果你想看看混淆矩阵，可以借助caret包实现：
caret::confusionMatrix(Linear_train_pred_df$diabetes,
                       Linear_train_pred_df$status_pred,
                       mode = "everything")
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction pos neg
##        pos 312  38
##        neg  82 105
##                                           
##                Accuracy : 0.7765          
##                  95% CI : (0.7389, 0.8111)
##     No Information Rate : 0.7337          
##     P-Value [Acc > NIR] : 0.01292         
##                                           
##                   Kappa : 0.4792          
##                                           
##  Mcnemar's Test P-Value : 8.661e-05       
##                                           
##             Sensitivity : 0.7919          
##             Specificity : 0.7343          
##          Pos Pred Value : 0.8914          
##          Neg Pred Value : 0.5615          
##               Precision : 0.8914          
##                  Recall : 0.7919          
##                      F1 : 0.8387          
##              Prevalence : 0.7337          
##          Detection Rate : 0.5810          
##    Detection Prevalence : 0.6518          
##       Balanced Accuracy : 0.7631          
##                                           
##        'Positive' Class : pos             
## 

# 内容非常全面，我们就不解读了。
# 我们直接把剩下的核函数在训练集、测试集中的结果都提取出来，方便接下来使用。

# 提取4种核函数分别在训练集、测试集的结果
Linear_test_pred_df <- getres(svmLinear, test)
Poly_train_pred_df <- getres(svmPoly, train)
Poly_test_pred_df <- getres(svmPoly, test)

Radial_train_pred_df <- getres(svmRadial, train)
Radial_test_pred_df <- getres(svmRadial, test)

Sigmoid_train_pred_df <- getres(svmSigmoid, train)
Sigmoid_test_pred_df <- getres(svmSigmoid, test)

# 接下来又是大家喜闻乐见的画图环节，就选大家最喜欢的ROC曲线吧。
# 关于这个ROC曲线，我一共写了十几篇推文，应该是全面覆盖了，大家还不会的去翻历史推文吧。
# 其实这里你也可以写个函数哈，大神们都说只要重复超过3遍的都建议写函数实现…

# 首先构建训练集中4个ROC对象
roc_train_linear <- roc(Linear_train_pred_df$diabetes, 
                        Linear_train_pred_df$pos,
                        auc=T
                        )
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases

roc_train_Poly <- roc(Poly_train_pred_df$diabetes, 
                      Poly_train_pred_df$pos,
                      auc=T
                        )
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases

roc_train_Radial <- roc(Radial_train_pred_df$diabetes, 
                        Radial_train_pred_df$pos,
                        auc=T
                        )
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases

roc_train_Sigmoid <- roc(Sigmoid_train_pred_df$diabetes, 
                        Sigmoid_train_pred_df$pos,
                        auc=T
                        )
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases

# 然后我们准备4种颜色，这种小代码，建议大家记住，因为使用很高频，它可以直接给你十六进制颜色代码，复制粘贴就可以使用了！
RColorBrewer::brewer.pal(4,"Set1")
## [1] "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3"

# 然后就是把训练集中，4种核函数的ROC曲线画在1张图上：
plot.roc(Linear_train_pred_df$diabetes, 
         Linear_train_pred_df$pos,
         col="#1c61b6",legacy=T,lwd=2)
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
lines.roc(Poly_train_pred_df$diabetes, 
          Poly_train_pred_df$pos, col="#008600")
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
lines.roc(Radial_train_pred_df$diabetes, 
          Radial_train_pred_df$pos, col="#E41A1C")
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
lines.roc(Sigmoid_train_pred_df$diabetes, 
          Sigmoid_train_pred_df$pos, col="#984EA3")
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases

legend("bottomright", 
       legend=c(paste0("train_svmLinear AUC: ",round(roc_train_linear[["auc"]],3)),
                paste0("train_svmPoly AUC: ",round(roc_train_Poly[["auc"]],3)),
                paste0("train_svmRadial AUC: ",round(roc_train_Radial[["auc"]],3)),
                paste0("train_svmSigmoid AUC: ",round(roc_train_Sigmoid[["auc"]],3))
                ),
       col=c("#1c61b6", "#008600","#E41A1C","#984EA3"), 
       lwd=2)

# easy！看着还行。果然是高斯径向基核函数最牛逼，堪称万金油！

# 测试集的数据已经提取好了，直接用即可。还是写个函数吧….
# 构建测试集中4个ROC对象
roc_test <- lapply(list(Linear_test_pred_df,Poly_test_pred_df,
                        Radial_test_pred_df,Sigmoid_test_pred_df), function(x){
                          roc_res <- roc(x$diabetes, x$pos,auc=T)
                          }
                   )
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
roc_test[[1]]
## 
## Call:
## roc.default(response = x$diabetes, predictor = x$pos, auc = T)
## 
## Data: x$pos in 150 controls (x$diabetes pos) > 81 cases (x$diabetes neg).
## Area under the curve: 0.8414

# 然后把测试集中，4种核函数的ROC曲线画在一起：
plot.roc(Linear_test_pred_df$diabetes, 
         Linear_test_pred_df$pos,
         col="#1c61b6",legacy=T)
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
lines.roc(Poly_test_pred_df$diabetes, 
          Poly_test_pred_df$pos, col="#008600")
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
lines.roc(Radial_test_pred_df$diabetes, 
          Radial_test_pred_df$pos, col="#E41A1C")
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases
lines.roc(Sigmoid_test_pred_df$diabetes, 
          Sigmoid_test_pred_df$pos, col="#984EA3")
## Setting levels: control = pos, case = neg
## Setting direction: controls > cases

legend("bottomright", 
       legend=c(paste0("test_svmLinear AUC: ",round(roc_test[[1]][["auc"]],3)),
                paste0("test_svmPoly AUC: ",round(roc_test[[2]][["auc"]],3)),
                paste0("test_svmRadial AUC: ",round(roc_test[[3]][["auc"]],3)),
                paste0("test_svmSigmoid AUC: ",round(roc_test[[4]][["auc"]],3))
                ),
       col=c("#1c61b6", "#008600","#E41A1C","#984EA3"), 
       lwd=2)

# 结果看起来还不错哦。