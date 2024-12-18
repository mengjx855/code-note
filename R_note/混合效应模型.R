# Jinxin Meng, 20240911, 20240920 ------------------


pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)
library(lme4)

# 混合效应模型
# 假设我们研究的是某种治疗药物对病人血压的影响。我们有多个病人，每个病人每天测量一次血压，持续一个月。
# 我们想分析药物剂量（固定效应）对血压的影响，同时考虑个体差异（随机效应）。

# 数据解读：
# 固定效应：药物剂量（低、中、高）
# 随机效应：病人个体差异
# 应变量：血压
# 
# 混合效应模型分析
# 构建模型：模型把药物剂量作为固定效应变量。病人作为随机效应变量，捕捉个体之间的血压差异。
# 模型公式：血压=β0 + β1 × 剂量 + ui + ϵ，其中，β0是截距，β1是剂量的系数，ui是病人的随机效应，ϵ是误差项。
# 模型解释：固定效应系数β1表示剂量增加对血压变化的平均影响。随机效应ui表示每个病人相对于平均水平的差异。
# 模型结论：如果β1显著为负，说明更高的药物剂量可能有效降低血压。随机效应的方差帮助我们了解病人之间的变异程度。

# 模拟 ---------
# 我们会创建一个数据集，包含病人的ID、药物剂量和对应的血压测量值。我们将有三种药物剂量（低、中、高），
# 并假设不同的病人有不同的基线血压水平（随机效应）。

# 设置随机数种子
set.seed(123)

# 模拟数据集
n_patients <- 100  # 病人数量
n_obs_per_patient <- 10  # 每个病人的观测次数
n_total_obs <- n_patients * n_obs_per_patient  # 总观测量

# 病人ID
patient_id <- rep(1:n_patients, each=n_obs_per_patient)

# 模拟剂量（低、中、高）作为固定效应变量
dose <- rep(c("low", "medium", "high"), length.out = n_total_obs)

# 随机生成病人的基线血压（每个病人一个随机基线）
baseline_bp <- rnorm(n_patients, mean = 120, sd = 10)
# 每个病人的观测血压基于他们的基线
blood_pressure <- baseline_bp[patient_id] + 
  ifelse(dose == "low", 0, ifelse(dose == "medium", -5, -10)) +  # 固定效应（剂量效应）
  rnorm(n_total_obs, mean = 0, sd = 5)  # 添加噪声

# 创建数据框
data <- data.frame(patient_id=as.factor(patient_id), dose=as.factor(dose), blood_pressure=blood_pressure)

# 查看前几行数据
head(data)

# 接下来，我们使用lme4包中的lmer函数来拟合模型。模型的公式如下：
# 固定效应：药物剂量对血压的影响。
# 随机效应：不同病人的基线血压差异。
# 拟合混合效应模型
model <- lmer(blood_pressure ~ dose + (1 | patient_id), data = data)

# 查看模型结果
summary(model)

# 在summary中，你会看到模型的结果，特别是以下两部分：
# Fixed Effects (固定效应)：表示药物剂量对血压的平均影响。
# 例如，如果中剂量和高剂量的系数为负，则意味着它们相较于低剂量可以显著降低血压。
# Random Effects (随机效应)：表示病人之间的基线差异。Std.Dev表示病人基线血压的标准差。

# 随机效应模型 函数 公式怎么写？
# 混合效应模型的公式通常写成如下形式：
# yij = β0 + β1 × Xij + ui + ϵij
# 其中：
# yij 是第i个个体的第j次观测的应变量。
# β0 是截距项。
# β1 是固定效应的系数（比如药物剂量的效应）。
# Xij 是固定效应变量（比如剂量）。
# ui 是随机效应，表示个体i的随机变异（个体差异）。
# ϵij 是误差项，表示不可解释的随机误差。

在R lmer 函数中，混合效应模型可以用以下公式表示：
lmer(response ~ fixed_effects + (random_effects_structure | group), data = dataset)
# response：应变量，即你要预测或分析的变量。
# fixed_effects：固定效应变量，可以是一个或多个，表示对所有观测单位有相同影响的变量。
# random_effects_structure：随机效应结构，指定哪些变量的效应在组间是随机的。
# group：分组因素，表示在哪个层次上引入随机效应。

# lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy) 
# lmer(response ~ fixed_effects + (random1 | group1) + (random2 | group2), data = dataset)

# random_effects_structure指的是模型中随机效应的组成方式。它描述了哪些变量的效应在不同的组之间是随机变化的。
# 有什么作用呢？
# 随机截距：允许每个组有不同的基线值。例如，(1 | group) 表示每个组都有自己的截距。
# 随机斜率：允许自变量在不同组的影响不同。例如，(Days | group) 表示每个组对 Days 这个自变量有不同的响应。
# 只有随机截距：表示每个 Subject 有不同的截距。
(1 | Subject)
# 随机截距和斜率：表示每个 Subject 对 Days 有不同的截距和斜率。
(Days | Subject)


