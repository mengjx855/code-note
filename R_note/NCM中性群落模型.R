# 现代微生物生态学理论认为，微生物群落的形成过程中受到不同因素的影响，划分为确定性过程和随机性过程。如下图所示。
# 群落组装机制：
# 1. 随机性过程：
# a. 扩散限制：个体的迁移和/或在新地点建立个体(殖民化)受到限制，导致群落之间的结构更加不同。
# b. 漂移：由于固有的出生、死亡和繁殖的随机过程，一个群落内不同物种的相对丰度随着时间的推移而发生的关联物种特性的随机变化
# c. 均匀扩散：群落间非常高的分散率，使群落同质化
# 2. 确定性过程
# a. 异质选择：异质非生物和生物环境条件下的选择导致群落间结构的差异性增加
# b. 同质选择：同质的非生物和生物环境条件下的选择导致群落间更加相似
# Zhou et al.. 2017. Microbiology and Molecnlar Biology Reviews

# 容易理解，确定性过程主要指微生物在环境过滤和微生物种间作用影响下发生的群落改变。比如你往一块土壤里喷洒稀盐酸或放点蚯蚓进去，群落将向着确定性的方向演替。而随机性过程强调微生物的扩散和漂移，即微生物的随机扩散、定殖和正常的生老病死。比如存在一个与世隔绝的村庄，里面的人的生活很可能就属于随机过程；而如果有一天发洪水或大旱，又或者有贪婪的军队发现了他们，可想而知，他们的生活将变成确定性的。
# 群落组装机制是生态学中用来解释生物群落形成和演变的理论框架。它关注个体如何在特定环境中相互作用、竞争资源以及适应环境条件。群落的组装过程涉及多个生态学因素，包括生物间的相互作用、生物和环境之间的相互作用以及物种的迁移和扩散。
# 在群落组装的背景下，中性群落模型（neutral community model，NCM）成为生态学中一种重要的理论框架，强调了群落结构中的随机性和中性过程。中性群落模型假设微生物个体在生态学上是相同的，没有竞争优势，因此，群落的动态主要受到随机漂移和迁移的影响。这一模型在解释缺乏明显选择压力或环境梯度的生态系统中的微生物群体形成和多样性方面提供了一种理论基础。通过数学方程，中性群落模型模拟了微生物的随机分布和变化，突显了群体组装中的随机性和偶然性。这种理论框架促使我们重新审视生态系统中微生物群体动态的复杂性，强调了随机过程在群体组装中的重要性，为理解生物多样性和群体结构的形成提供了深刻的洞察。

# https://mp.weixin.qq.com/s?src=11&timestamp=1709778653&ver=5123&signature=7ugswtvtKj26ssERZ0IlVFbg0Xs7ApXsIQPsuVQ-*5fGz6Wg*HggE4hpLbtHwz74zclyjWPbcPbv420OOb6tN1lATlNeabWiX0zT9GbPN*7PNA5xLDbpYmdyz1jOMvje&new=1
# R2（R方）代表了中性群落模型的整体拟合优度，R2越高表明越接近中性模型，即群落的构建受随机性过程的影响越大，受确定性过程的影响越小。m量化了群落层面的迁移率（migration rate），该值对于每个群落成员都是统一的（与物种无关），m值越小说明整个群落中物种扩散越受限制，反之m值越高则表明物种受到扩散限制越低。蓝色实线表示NCM 的最佳拟合，蓝色虚线表示模型预测周围的 95% 置信区间。发生频率高于或低于 NCM 预测的 OTU 以不同的颜色显示。

# 生态位理论（niche-based theories）和中性理论（neutral-based theories）构成了理解微生物群落组装的两个重要且互补的机制。生态位过程认为微生物群落是由确定性的非生物因素（环境因素如pH、温度等）和生物因素（物种相互作用如竞争、捕食等）所决定，归因于微生物不同的栖息地偏好和适应能力。中性理论则认为，随机过程，如出生、死亡、迁移、物种形成和有限扩散，塑造了微生物群落结构，并假设微生物在类群的损失和收益之间呈现一种随机平衡。
# 尽管随机过程被认为在塑造微生物群落结构方面发挥着重要作用，但在以前，由于难以定义随机性和缺少用于表示随机性的方法，生态随机性在影响群落结构中的重要性却经常被忽视。随着微生物生态学的迅速发展，许多有关于微生物群落形成的概念性模型已经被提出。群落多样性和分布模式可以为群落组装的过程提供证据。本篇主要介绍Sloan等（2006）提出的中性群落模型（neutral community model，NCM），在量化中性过程的重要性时特别有用。

library(Hmisc)
library(minpack.lm)
library(stats4)

# 将 Sloan 等（2006）的中性模型拟合到一个物种或分类群的丰度表，并返回几个拟合统计数据。或者，将根据它们在元群落中的丰度返回每个分类群的预测出现频率
# 用非线性最小二乘法（Non-linear least squares，NLS）拟合模型参数
# 转置otu表格
spp <- t(spp)
# 计算总相对丰度的平均值
N <- mean(apply(spp, 1, sum))
# 计算每个物种的平均相对丰度
p.m <- apply(spp, 2, mean)
# 去除平均值为0的物种
p.m <- p.m[p.m != 0]
# 计算每个物种的相对丰度
p <- p.m/N
# 将原始数据二值化，表示物种的存在与否
spp.bi <- 1 * (spp > 0)
# 计算每个物种的出现频率
freq <- apply(spp.bi, 2, mean)
# 去除频率为0的物种
freq <- freq[freq != 0]
# 合并相对丰度和频率数据
C <- merge(p, freq, by = 0)
# 按照频率排序
C <- C[order(C[,2]),]
# 将结果转为数据框
C <- as.data.frame(C)
# # 去除包含0的行
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
# 提取相对丰度和频率的数据
p <- C.0[,2]
freq <- C.0[,3]
# 为数据命名
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
# # 计算d的值
d = 1/N
# 使用非线性最小二乘法（NLS）拟合模型参数 m（或 Nm）
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
# 输出拟合结果
m.fit
# 计算 m 的置信区间
m.ci <- confint(m.fit, 'm', level=0.95)
# 计算预测的频率值
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
# 使用 Wilson 方法计算置信区间
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
# 计算 R2 值
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
# Rsqr  #获取模型的 R2

# 输出 3 个统计结果数据表，包括各物种或分类群的平均相对丰度（p.csv）、出现频率（freq.csv）和预测的出现频率（freq.pred.csv）
# write.csv(p, file = "p.csv")
# write.csv(freq, file = "freq.csv")
# write.csv(freq.pred, file = "freq.pred.csv")
# p 是平均相对丰度（mean relative abundance）
# freq 是出现频率（occurrence frequency）的观测值
# freq.pred 是出现频率（occurrence frequency）的预测值，即中性模型的拟合值

# 绘制统计图
# 定义数据框 bacnlsALL，包含 p、freq、freq.pred、以及 pred.ci 的数据
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
# 定义点的颜色
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'#出现频率低于中性群落模型预测的部分
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'#出现频率高于中性群落模型预测的部分
# 加载 grid 包
library(grid)
# 创建新页面
grid.newpage()
# 设置视图
pushViewport(viewport(h = 0.6, w = 0.6))
# 设置数据视图范围
pushViewport(dataViewport(xData = range(log10(bacnlsALL$p)), yData = c(0, 1.02), extension = c(0.02, 0)))
# 绘制矩形
grid.rect()
# 绘制散点图
grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, gp = gpar(col = inter.col, cex = 0.7))
# 添加 y 轴和 x 轴
grid.yaxis()
grid.xaxis()
# 绘制预测的曲线
grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp = gpar(col = '#0a71b4', lwd = 2), default = 'native')
# 绘制置信区间
grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp = gpar(col = '#0a71b4', lwd = 2, lty = 2), default = 'native')
grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp = gpar(col = '#0a71b4', lwd = 2, lty = 2), default = 'native')
# 添加文字标签
grid.text(y = unit(0, 'npc') - unit(2.5, 'lines'), label = 'Mean Relative Abundance (log10)', gp = gpar(fontface = 2))
grid.text(x = unit(0, 'npc') - unit(3, 'lines'), label = 'Frequency of Occurrence', gp = gpar(fontface = 2), rot = 90)
# 定义绘制文字的函数
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=", round(Rsqr, 3), "\n", "Nm=", round(coef(m.fit) * N)), x = x[j], y = y[i], just = just)
}
# 添加文字标签
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

# 根据计算结果，R2=0.424。Nm=24417，而N=41373，相除得到m=0.59。表明该组随机性过程占比不低，其中扩散限制较低。
# R2代表了中性群落模型的整体拟合优度，R2越高表明越接近中性模型，即群落的构建受随机性过程的影响越大，受确定性过程的影响越小。
# N描述了元群落规模（metacommunity size），为每个样本中所有OTU的总丰度。
# m量化了群落层面的迁移率（migration rate），该值对于每个群落成员都是统一的（与物种无关），m值越小说明整个群落中物种扩散越受限制，反之m值越高则表明物种受到扩散限制越低。
# Nm是元群落规模（N）与迁移率（m）的乘积 (Nm = N*m)，量化了对群落之间扩散的估计，决定了发生频率和区域相对丰度之间的相关性。

plotdat <- data.frame(p, freq, freq.pred, pred.ci[,2:3]) %>% 
  dplyr::mutate(class = ifelse(freq <= Lower, "lower", ifelse(freq >= Upper, "upper", "inter")))
colors <- structure(c("#000000", "#a52a2a", "#29a6a6"), names = c("inter", "lower", "upper"))
subtitle <- paste0("Coefficient: Rsqr = ", round(Rsqr, 4), ", m = ", round(coef(m.fit), 4),
                   ", mN = ", round(coef(m.fit) * N, 4))
ggplot(plotdat) +
  geom_point(aes(x = log10(p), y = freq, color = class), size = 1.2, show.legend = F) +
  geom_line(aes(x = log10(p), y = freq.pred, group = 1), linewidth = .6, color = "blue") +
  geom_line(aes(x = log10(p), y = Upper, group = 1), linewidth = .6, color = "blue", lty = "dashed") +
  geom_line(aes(x = log10(p), y = Lower, group = 1), linewidth = .6, color = "blue", lty = "dashed") +
  scale_color_manual(values = colors) +
  labs(x = "Mean Relative Abundance (log10)", y = "Frequency of Occurance",
       title = "Neutral community model analysis", subtitle = subtitle) +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        axis.text = element_text(size = 8, color = "#000000"),
        axis.title = element_text(size = 8, color = "#000000"),
        plot.title = element_text(size = 10, color = "#000000"),
        plot.subtitle = element_text(size = 8, color = "#000000"),
        legend.text = element_text(size = 8, color = "#000000"),
        legend.title = element_text(size = 8, color = "#000000"),
        panel.grid = element_blank(),
        aspect.ratio = 3/4)
ggsave("xx.pdf", height = 4, width = 4)



