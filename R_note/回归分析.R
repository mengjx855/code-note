#### 线性拟合 ####
# 有这样的一个数据，每个组中有多个观测值，使用ggplot2拟合出一条线来，脚本如下
# sample  val   group
# s1      0.12  g1
# s2      0.32  g2
# s3      0.11  g3
# s4      0.22  g4
# ....
ggplot(dat, aes(group, val)) +
  geom_boxplot(position = position_dodge(), fill = "transparent", width = .6, color = "#000000", linewidth = .3, outlier.shape = NA) +
  geom_jitter(aes(color = group), size = 1.5, width = .3, show.legend = F) + 
  stat_smooth(aes(group = 1), level = 0.95, method = "lm", color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 3, label.x.npc = 0, label.y.npc = 1, color = "#8da0cb")
# 这时，拟合曲线，置信区间以及直线的方程，rr和pval都出来了。怎么算的呢？怎么能自己写脚本计算这个值
# 例如这个组中的group有四个组，线性拟合需要把这个分组的列也变成数值型的值，不能是字符串或者因子，如下
group_order <- c("g1", "g2", "g3", "g4")
dat <- mutate(dat, group = factor(group, group_order, labels = seq_len(9)) %>% as.numeric())
lm_res <- lm(val ~ group, data = dat) # formula: y ~ x
# Call:
# lm(formula = val ~ group, data = dat)
# Coefficients:
# (Intercept)        group  
#    0.470772     0.002134  
# Coefficients表示拟合后的系数，Intercep为拟合直线的截距，右边的为斜率，也就是说该拟合的直线方程为y = 0.002134x+0.470772
# 如何知道rr和pval呢？
lm_sum <- summary(lm_res)
# Call:
# lm(formula = val ~ group, data = dat)
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.27175 -0.04747  0.01743  0.05443  0.15444 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.470772   0.017650  26.673   <2e-16 ***
# group       0.002134   0.003096   0.689    0.492    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.08049 on 103 degrees of freedom
# Multiple R-squared:  0.004589,	Adjusted R-squared:  -0.005075 
# F-statistic: 0.4748 on 1 and 103 DF,  p-value: 0.4923
# summary 结果解析：
# Residuals为残差，结果给出了参数的四分位数以及极值的信息。
# Coefficients 除了给出了截距和斜率值，外还提供了很多信息。由数据可以看出，斜率的T检验的值为12.525 (这个值大于t0.05)，
# 对应的p值为1.25e-08小于0.05，从而拒绝原假设，表明x对y的预测贡献了信息。后面的星号越多，线性关系越显著。
# R-squared 为决定系数(rr)，为0.004589
# 如何获得回归曲线的p值？需要使用pf函数
pf(lm_sum$fstatistic[1],lm_sum$fstatistic[2],lm_sum$fstatistic[3], lower.tail = F)
# 这个时候就得到所有的数据了
# 虽然通过summary函数可以得到线性模型的很多信息，但是有的时候信息太多，反而会扰乱我们的分析，我们只需要某一个具体的值，例如决定系数等等。可以参考下表中总结了对拟合线性模型非常有用的函数
# summary|展示拟合的详细结果
# coefficients|列出拟合模型的模型参数截距项和斜率
# fitted|列出拟合模型的预测值
# residuals|列出拟合模型的残差值，生成一个拟合模型的方差分析,或者比较两个或更多拟合模型的
# anova|方差分析表
# vcov|列出模型参数的协方差矩阵
# AIC|输出赤池信息统计量
# plot|生成评价拟合模型的诊断图
# predict|用拟合模型对新的数据集预测响应变量值
# confint|给定回归模型中参数的置信区间，默认为95%区间

#### 简单线性回归 ####
data <- read.csv("F:/Git/Regression-model-master/Regression-model-master/data/student_data.csv")
data <- data[,-1]
lm_data <- na.omit(data) # 删除列中有缺失值的数据

attach(lm_data) # 我感觉就是把一个数据搞成数据框，每一列是一个向量，列名就是向量的名称，可以在后头的使用中直接引用
plot(height, weight, main = "scatter plot") # 散点图可视化身高体重线性关系
abline(h = mean(weight), col = "red")
abline(v = 170, col = "red")

model <- lm(weight ~ height, data = lm_data) # 最小二乘法拟合简单线性回归
# Call:
# lm(formula = weight ~ height, data = lm_data)
# Coefficients:
# (Intercept)       height
#     -86.318        0.903
# β0 = -86.318 当人的身高为0时，体重为-86.318
# β1 = 0.903 身高每升高1cm，体重增加0.903kg

model$Coefficients # β0和β1
model$residuals # 残差，多少个样本对饮多少残差，残差=真实值-估计值
model$rank # β的数量，这里是两个，一个是β0，一个是β1
model$fitted.values # 样本预测值
model$df.residual # 残差的自由度

abline(model, col="red") # 将模型在图中画出

summary(model)
# Call:
# lm(formula = weight ~ height, data = lm_data)
# Residuals:
#     Min      1Q  Median      3Q     Max
# -23.029 -11.243  -5.614   3.178  87.068
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -86.318     22.819  -3.783 0.000193 ***
# height         0.903      0.134   6.740 1.02e-10 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Residual standard error: 18.63 on 260 degrees of freedom
# Multiple R-squared:  0.1487,    Adjusted R-squared:  0.1455
# F-statistic: 45.43 on 1 and 260 DF,  p-value: 1.018e-10
# Estimate这列表示的是求出来的参数β0和β1；
# Residuals部分显示的是所有样本的残差的五分位数，理想状态下，median为0，一分位点和三分位点的绝对值应该是相等的；
# RSE表示随机误差项，有公式计算；
# 260表示残差的自由度；
# Std. Error表示对应参数的标准差，β0和β1的标准差是由公式计算出来的，做假设检验的时候会用；
# 假设检验用的三种方法：区间估计（算出95%的把握，区间要是包含0，支持原假设<原假设就是X和Y没有线性关系>），p值法估计（），F检验法；

# 假设检验：身高和体重是否呈线性关系
# 1） 区间估计法，公式中用的标准误已经算出来，在summary表格里
t = -qt(0.025, df = model$df.residual) # 都是有公式计算的
c(0.903 - t*0.134, 0.903 + t*0.134)
# 2）p值法，认为身高服从正态分布，但是数据不能代表全部，所以最终认为是服从t分布，公式中用的t值已经算出来，在summary表格里
(1 - pt(6.740, df = model$df.residual))*2 # 都是有公式计算的
# 3）F值法
anova(model) # 先根据anova表找出SSE和SST
# Analysis of Variance Table
# Response: weight
#           Df Sum Sq Mean Sq F value    Pr(>F)    
# height      1  15771 15770.7  45.426 1.018e-10 ***
# Residuals 260  90265   347.2
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# [height, Sum Sq]是SSR，SSR除以[height, Df](自身自由度)就是MSR；
# [Residuals, Sum Sq]是SSE，SSE除以[Residuals, Df](自身自由度)就是MSE；
# MSR除以MSE就是F value；
# 然后根据F分布，检验F值的p值；
1 - pf(45.426, df1 = 1, df2 = 260)

# 估计：β
# 1）点估计
model$coefficients # 模型的β0和β1约等于实际中的β0和β1
# 2）区间估计
summary(model)
t = -qt(0.025, df = 260) # 估计斜率
c(0.903 - t*0.134, 0.903 + t*0.134) # 有公式计算
# [1] 0.6391366 1.1668634 # 有95%的把握认为斜率是出现在这个区间内
confint(model) # 这个命令能直接进行区间估计
#                    2.5 %     97.5 %
# (Intercept) -131.2515986 -41.383962 # β0的区间估计情况
# height         0.6391812   1.166827 # β1的区间估计情况

# 估计：y, 身高为165中国人的体重
# 1）点估计
model$coefficients
-86.3177802+0.903*165
# 2）区间估计
x = data.frame(height = c(165,180))
predict(model, newdata = x , interval = "predict") # 预测区间
predict(model, newdata = x , interval = "confidence") # 置信区间

# 在散点图中绘制出回归模型的置信区间和预测区间
plot(height, weight, main = "scatter plot")
abline(model, col = "red")

x = data.frame(height = seq(150, 200, by = 1)) # 预测150-200的体重
conf_interval = predict(model, newdata = x, interval = "confidence")
conf_interval
lines(x$height, conf_interval[,2], col = "blue", lty = 2)
lines(x$height, conf_interval[,3], col = "blue", lty = 2)


####模型诊断####
data <- read.csv("F:/Git/Regression-model-master/Regression-model-master/data/student_data.csv")
data <- data[,-1]
lm_data <- na.omit(data) # 删除列中有缺失值的数据

attach(lm_data) # 我感觉就是把一个数据搞成数据框，每一列是一个向量，列名就是向量的名称，可以在后头的使用中直接引用
model = lm(weight ~ height, data = lm_data) # 最小二乘法拟合简单线性回归

plot(height, weight, main = "scatter plot") # 散点图可视化身高体重线性关系
abline(model, col = "red")

# residual plot（残差图）
# x轴为样本的预测值，y为样本的残差，残差要满足随机性，正态性和等方差性
par(mfrow = c(1,2))
plot(x = model$fitted.values, y = model$residuals)
abline(h = 0, col = "red") # 满足随机性，但是不满足等方差性
plot(model, which = 1) # plot线性模型会生成四张图，Residuals vs Fitted, Normal QQ, Scale-Location, Residuals vs Leverage. which选择输出的图是第几个

# Scale-Location Graph 位置尺度图
plot(model$fitted.values, sqrt(abs(scale(model$residuals)))) # 只是y轴的残差标准化后的绝对值再开方
plot(model, which = 3)

# QQ图，可以判断数据的分布情况
hist(weight) # 直方图
hist(log(weight)) # 直方图
qqnorm(weight) # qq图
qqline(weight)
plot(model,which = 2)

# 异常点处理
plot(model,which=4)
plot(model,which=5)


####非线性回归####
nonlinear.student.data <- read.csv("F:/Git/Regression-model-master/Regression-model-master/data/student_data.csv", row.names = 1)
nonlinear.student.data <- na.omit(nonlinear.student.data)
head(nonlinear.student.data)

attach(nonlinear.student.data) # 散点图可视化身高体重线性关系
plot(height, weight, main = "scatter plot")

# 1）多项式回归模型建模
plot(height, weight, main = "scatter plot")

# 第一种方法：
poly.model = lm(weight ~ height + I(height^2) + I(height^3)) # 我也没闹清楚I()函数是咋用，可能是把内容作为整体的意思
poly.model
# 第二种方法：生成多项式
poly.model = lm(weight ~ poly(height, 100, raw = TRUE)) # raw要是T，否则建立出来的是正交多项式

x <- seq(150, 200, by = 0.1) # 可视化多项式回归
lines(x, predict(poly.model, newdata = data.frame(height = x)), col = "red")

conf_interval = predict(poly.model, data.frame(height = x), interval = "confidence") # 加置信区间
lines(x, conf_interval[,2], col = "blue", lty = 2)
lines(x, conf_interval[,3], col = "blue", lty = 2)

predict_interval = predict(poly.model, data.frame(height = x), interval = "predict") # 加预测区间
lines(x, predict_interval[,2], col = "green", lty = 2)
lines(x, predict_interval[,3], col = "green", lty = 2)

# 如何判断polynomial degree
anova(poly.model) # 发现简单线性回归就够了，二次项和三次项均不显著

# 2）把身高看成分类型变量回归
cut(height, 3)
step.regression = lm(weight ~ cut(height, 3))
step.regression

plot(height, weight ,main = "scatter plot") # 可视化分类型变量回归
lines(sort(height), sort(step.regression$fitted.values), col = "red", lty = 2, lwd = 2)

# 3）多分类变量回归
plot(height, weight, main = "scatter plot")
multi.factor.lm = lm(weight ~ cut(height, 3) + gender)
multi.factor.lm

# 4）分类型数值型混合回归
mix.model = lm(weight ~ height + gender)
mix.model

plot(height, weight, main = "scatter plot") # 可视化混合回归
x <- seq(150, 200, by = 0.1)
lines(x, predict(mix.model,newdata = data.frame(height = x, gender = "Male")), col = "red", lty = 2, lwd = 2)
lines(x, predict(mix.model,newdata = data.frame(height = x, gender = "Female")), col = "Blue", lty = 2, lwd = 2)
legend(180, 140, legend = c("Male", "Female"), col = c("red", "blue"), lty = c(2, 2))

# 5）分类型数值型混合回归（含交互作用）
# 第一种方法：
mix.model = lm(weight ~ height*gender)
mix.model
# 第二种方法：
mix.model = lm(weight ~ height+gender+height:gender)
mix.model

plot(height, weight, main = "scatter plot") # 可视化混合回归
x <- seq(150, 200, by = 0.1)
lines(x, predict(mix.model, newdata = data.frame(height = x,gender = "Male")), col = "red", lty = 2, lwd = 2)
lines(x, predict(mix.model, newdata = data.frame(height = x,gender = "Female")), col = "Blue", lty = 2, lwd = 2)
legend(180, 140, legend = c("Male", "Female"), col = c("red", "blue"), lty = c(2, 2))

####多元线性回归####
housing_data <- read.csv("F:/Git/Regression-model-master/Regression-model-master/data/housing_data.csv", row.names = 1)
multi.data <- housing_data[,-1]

plot(multi.data) # 可视化数据，希望x和y呈线性关系，而不同的x之间不要有线性关系，否则引起共线性关系
cor(multi.data) # 相关性同样如此

# install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
chart.Correlation(multi.data, histogram = TRUE, pch = 19) # 这个函数可以检查x和y之间的关系

model <- lm(price ~ sqft.living+bedrooms+bathrooms, data = multi.data) # 建立多元线性回归模型
model

par(mfrow=c(2,2))
plot(model) # 回归诊断，发现数据呈现偏态，最好的办法就是y变量取log

model <- lm(log(price) ~ sqft.living+bedrooms+bathrooms, data = multi.data) # 建立新的多元线性回归模型
model

# 解释多元线性回归模型
model
exp(12.32446) # β0，斜距

par(mfrow=c(1,1)) # 解释其他β，在两项不变的情况下，第三项x增加一个单位，y增加一个项系数的单位 
bathrooms <- 1:10
bedrooms <- rep(mean(multi.data$bedrooms), 10) 
sqft.living <- rep(mean(multi.data$sqft.living), 10)
dat.new <- data.frame(bedrooms, bathrooms, sqft.living)
plot(bathrooms, exp(predict(model, newdata = dat.new)), type = "l", lwd = 3, ylab = "price", xlab = "bedrooms", cex.lab = 2)

# 研究X变量是否都显著需要看ANOVA表
summary(model)
# Call:
#   lm(formula = log(price) ~ sqft.living + bedrooms + bathrooms, 
#      data = multi.data)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.47558 -0.29821  0.00659  0.26738  1.23097 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    1.232e+01  1.855e-02 664.539  < 2e-16 ***
#   sqft.living  4.000e-04  8.495e-06  47.082  < 2e-16 ***
#   bedrooms    -5.816e-02  6.439e-03  -9.032  < 2e-16 ***
#   bathrooms    5.070e-02  9.367e-03   5.412 6.42e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3867 on 6976 degrees of freedom
# Multiple R-squared:  0.4625,	Adjusted R-squared:  0.4622 
# F-statistic:  2001 on 3 and 6976 DF,  p-value: < 2.2e-16
# 这里表的最后一列Pr(>|t|)假设检验的p值表示，在其他两项存在的情况下，该项的显著性，等于ANOVA的第三个假设
# F-statistic是F总检验

anova(model) 
# Analysis of Variance Table
# Response: log(price)
#              Df  Sum Sq Mean Sq  F value    Pr(>F)    
# sqft.living   1  882.53  882.53 5902.754 < 2.2e-16 ***
#   bedrooms    1   10.40   10.40   69.577 < 2.2e-16 ***
#   bathrooms   1    4.38    4.38   29.294 6.425e-08 ***
#   Residuals   6976 1042.99    0.15                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 其实就是三个假设检验，但是是有顺序的
# 第一个假设检验H0就是平米数对于log（价格）不显著，Ha为平米数对于log（价格）显著
# 第二个假设检验H0认为在平米数这个变量存在的情况下，bedrooms数对log（价格）不显著，Ha为bedrooms数对于log（价格）显著
# 第三个假设检验H0认为在平米数和bedrooms数这两个变量存在的情况下，bathrooms数对log（价格）不显著，Ha为bathrooms数对于log（价格）显著
# F = MSR/MSE

# overall F test:
MSR = (882.53+10.40+4.38)/3
MSE = 1042.99/6976
F = MSR/MSE
1 - pf(F, df1 = 3, df2 = 6976)

# 预测
model
x = data.frame(sqft.living = 1000, bedrooms = 3, bathrooms = 3)
exp(predict(model, x)) # 点估计
exp(predict(model, x, interval = "predict")) # 区间估计，预测区间
exp(predict(model, x, interval = "confidence")) # 区间估计，置信区间

# 1）pearson 相关系数 处理多重共线性
cor(multi.data)

# 2）VIF 方差膨胀系数(variance inflation factor)
partial.regrs <- lm(sqft.living ~ bedrooms+bathrooms, data = multi.data ) # 辅助回归
summary(partial.regrs) # 取出R2
VIF = 1/(1-0.6281) # 如果大于10就认为是sqft.living和bedrooms、bathrooms相关，就要去掉sqft.living避免多重共线性

# install.packages("car")
library(car)
vif(model) # 直接计算VIF

# 3) 偏残差图
#判断X变量和Y的线性关系
par(mfrow=c(1,3))
plot(multi.data$bedrooms,log(multi.data$price))
plot(multi.data$bathrooms,log(multi.data$price))
plot(multi.data$sqft.living,log(multi.data$price))

# 判断在bedrooms变量和bathrooms存在的情况下，sqft.living和log(price)是否呈线性关系
# sqft.living为X轴，不含sqft.living变量拟合出回归模型的残差为Y轴，做散点图
par(mfrow=c(1,1))
partial.model = lm(log(price) ~ bathrooms+bedrooms,data = multi.data)
plot(multi.data$sqft.living, partial.model$residuals)

# 偏残差图的包
#install.packages("car")
library(car)
model = lm(log(price) ~ sqft.living+bedrooms, data = multi.data)
crPlots(model)


####异常点处理####
#建立回归模型
mod <- lm(SavingsRate ~ Pop15+Pop75+DPIgrowth, data = savings)
anova(mod)

# 1）高杠杆点 high leverage point
hatv <- hatvalues(mod)
plot(hatv, type="h", ylab="Leverages")
# 标出界限
p <- 4
n <- nrow(savings)
abline(h=(2*p/n), lwd=2, col="red") 
# 标出高杠杆点
countries <- row.names(savings)
text(c(1:50)[hatv>(2*p/n)], hatv[hatv>(2*p/n)],savings$Country[hatv>(2*p/n)] , col="dodgerblue", cex=1.5)

# 2）异常点 outlier
par(mfrow=c(1,2))
plot(mod,which=1)
plot(mod$fitted.values,rstandard(mod))
abline(h = 0, col="red")
# 标出 outlier 
countries <- row.names(savings)
text(mod$fitted.values[abs(rstandard(mod))>2], rstandard(mod)[abs(rstandard(mod))>2],savings$Country[abs(rstandard(mod))>2] , col="dodgerblue", cex=1)
# 根据标准QQ图找出异常点
par(mfrow=c(1,1))
plot(mod,which=2)

# 3）cook distance找强影响点（high influential point）
cooks.distance(model)
plot(cooks.distance(model), type="h", ylab="cook.distance")
# 标出影响点
countries <- row.names(savings)
text(c(1:50)[which.max(cooks.distance(model))],cooks.distance(model)[which.max(cooks.distance(model))], savings$Country[which.max(cooks.distance(model))] , col="dodgerblue", cex=1)
# 用R语言内置函数画cook distance
plot(model,which=4)

# 在一幅图中展示前三个结果
plot(model,which=5)

# 4）斜率的改变率 找强影响点（high influential point）
par(mfrow=c(2,2))
# intercept change percentage
intercept.change = dfbeta(mod)[,1]/coef(mod)[1]*100
plot(intercept.change, type = "h", ylab="intercept", pch=16, cex.lab=1.5)
abline(h=0)
text(c(1:50)[which.max(abs(intercept.change))],intercept.change[which.max(abs(intercept.change))],savings$Country[which.max(abs(intercept.change))], col="dodgerblue", cex=1)

# Pop15 change percentage
Pop15.change = dfbeta(mod)[,2]/coef(mod)[2]*100
plot(Pop15.change, type = "h", ylab="Pop15", pch=16, cex.lab=1.5)
abline(h=0)
text(c(1:50)[which.max(abs(Pop15.change))],Pop15.change[which.max(abs(Pop15.change))],savings$Country[which.max(abs(Pop15.change))], col="dodgerblue", cex=1)

# Pop75 change percentage
Pop75.change = dfbeta(mod)[,3]/coef(mod)[3]*100
plot(Pop75.change, type = "h", ylab="Pop75", pch=16, cex.lab=1.5)
abline(h=0)
text( c(1:50)[which.max(abs(Pop75.change))],Pop75.change[which.max(abs(Pop75.change))],savings$Country[which.max(abs(Pop75.change))], col="dodgerblue", cex=1)

# DPIgrowth change percentage
DPIgrowth.change = dfbeta(mod)[,4]/coef(mod)[4]*100
plot(DPIgrowth.change, type = "h", ylab="DPIgrowth", pch=16, cex.lab=1.5)
abline(h=0)
text(c(1:50)[which.max(abs(DPIgrowth.change))],DPIgrowth.change[which.max(abs(DPIgrowth.change))],savings$Country[which.max(abs(DPIgrowth.change))], col="dodgerblue", cex=1)


####线性模型特征选择####
#数据从这个ISLR包中获取
install.packages("ISLR")
library("ISLR")
View(Credit)
attach(Credit)

#建立回归模型
model = lm(Balance~.,data = Credit)

#1）逐步回归
#install.packages("leaps")
library(leaps)
regfit.full=regsubsets(Balance~.,data=Credit,method = "exhaustive",nvmax = 12)
reg.summary = summary(regfit.full)
summary(regfit.full)

#2）adjusted R^2
par(mfrow=c(1,2))
plot(reg.summary$rsq,xlab="Number of Variables",type="b",ylab="RSq")
plot(reg.summary$adjr2,xlab="Number of Variables", ylab="Adjusted RSq",type="b",ylim=c(0.953,0.954))
par(mfrow=c(1,1))
plot(regfit.full,scale="adjr2")
coef(regfit.full,id =7)

#3） BIC and Mellow'CP
plot(regfit.full,scale="bic")
coef(regfit.full,id =4)
which.min(reg.summary$cp)
coef(regfit.full,id =6)


# 4）Validation selection
cv_set=sample(1:nrow(Credit),nrow(Credit)*3/4)
regfit.best=regsubsets(Balance~.,data=Credit[cv_set,],nvmax=12)
test.mat=model.matrix(Balance~.,data=Credit[-cv_set,])
val.errors=rep(NA,11)

for(i in 1:12){
  coefi=coef(regfit.best,id=i)
  pred=test.mat[,names(coefi)]%*%coefi
  val.errors[i]=mean((Credit$Balance[-cv_set]-pred)^2)
}

which.min(val.errors)
coef(regfit.full,id=4)

#5）10-fold Validation selection
k=10
folds=sample(1:k,nrow(Credit),replace=TRUE)
out.cv = matrix(0, nrow = k, ncol = 12)
colnames(out.cv)=1:12

predict.regsubsets = function(object, newdata, id, ...){
  form = as.formula(object$call[[2]])
  coefi = coef(object, id)
  test.mat = model.matrix(form, newdata)
  pred = test.mat[ , names(coefi)]%*%coefi
  return(pred)
}

for (j in 1:k) {
  best.fit=regsubsets(Balance~.,data=Credit[folds!=j,],nvmax=12)
  for (i in 1:12) {
    pred = predict.regsubsets(best.fit, Credit[folds == j, ], i)
    MSE = mean((Credit$Balance[folds == j] - pred)^2)
    out.cv[j, i] = MSE
  }
}
out.cv
which.min(apply(out.cv, 2, mean))
coef(regfit.full,id=4)

# 6）Ridge and Lasso
gridLambda=10^seq(5,-2,length=100)
install.packages("glmnet")
library(glmnet)
train.validation.x=model.matrix(Balance~.,data=Credit)
train.validation.y=Credit$Balance
lasso.mod=cv.glmnet(train.validation.x,train.validation.y,alpha=0,lambda=gridLambda)
lasso.mod$cvm
# 给定最佳的lambda值，算coefficient
bestlam=lasso.mod$lambda.min
predict(lasso.mod,s=gridLambda[1],type="coefficients")
mean(Balance)

# 7）降维 PCR
library(pls)
pcr.fit <- pcr(Balance~.,data=Credit,scale=TRUE,validation="CV")
summary(pcr.fit)
validationplot(pcr.fit,val.type="MSEP")
predict(pcr.fit,newdata = Credit[1:5,],ncomp=11)

# 8）降维 PLS
library(pls)
set.seed(1)
pls.fit <- plsr(Balance~.,data=Credit,scale=TRUE,validation="CV")
summary(pls.fit)
validationplot(pls.fit,val.type="MSEP")
predict(pls.fit,newdata = Credit[1:10,],ncomp=4)





