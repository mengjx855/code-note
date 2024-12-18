# Hierarchical Modelling of Species Communities (Hmsc) is a flexible framework for Joint Species Distribution Modelling (JSDMs). 
# The framework can be used to relate species occurrences or abundances to environmental covariates, species traits and phylogenetic relationships. 
# JSDMs are a special case of species distribution models (SDMs) that take into account the multivariate nature of communities which allows us to 
# estimate community level responses as well capture biotic interactions and the influence of missing covariates in residual species associations

# The Hmsc package contains functions to fit JSDMs, analyze the output and to generate predictions with these JSDMs. 
# The obligatory data for a HMSC analysis includes a matrix of species occurrences or abundances and a matrix of environmental covariates.
# Optionally, the user can include information species traits, phylogenetic relationships and information on the spatiotemporal context of the
# sampling design to account for dependencies among the sampling units.

# 物种群落层次模型（Hmsc）是联合物种分布模型（JSDM）的一个灵活框架。该框架可用于将物种出现或丰度与环境协变量、物种特征和系统发育关系联系起来。
# JSDM是物种分布模型（SDM）的一个特例，它考虑了群落的多元性质，使我们能够估计群落水平的反应，并捕捉生物相互作用和剩余物种关联中缺失协变量的影响。

# https://github.com/hmsc-r/HMSC/


# 1. Hmsc 的基本概念
# Hierarchical Modelling of Species Communities (Hmsc) 是一个灵活的框架，用于联合物种分布建模（Joint Species Distribution Modelling, JSDM）。
# JSDMs 是一种特殊的物种分布模型（Species Distribution Models, SDMs），它考虑了物种群落的多变量特性。
# 2. JSDM 的特点
# JSDMs 允许我们将物种的出现或丰度与环境变量（environmental covariates）、物种特征（species traits）和系统发育关系（phylogenetic relationships）联系起来。
# 这种模型特别之处在于，它能够估计社区层面的响应（community level responses），并捕捉生物相互作用（biotic interactions）以及在残差物种关联中缺失的协变量的影响。
# 3. Hmsc 包的功能
# Hmsc 包含了一些函数，用于：
# 拟合 JSDMs（建立模型）
# 分析模型输出（理解和解释结果）
# 生成基于这些 JSDMs 的预测（做出未来的推测）
# 4. 数据要求
# 进行 HMSC 分析所必需的数据包括：
# 一个物种出现或丰度的矩阵（matrix of species occurrences or abundances）：这是一个表格，其中行通常代表不同的采样单位（如地点或时间），列代表不同的物种，表中的值表示物种的出现（是/否）或丰度（数量）。
# 一个环境协变量的矩阵（matrix of environmental covariates）：这是另一个表格，行与前面的矩阵相同，列包含不同的环境因素（如温度、湿度、土壤类型等）。
# 5. 可选数据
# 用户还可以选择包含以下信息来改善模型：
# 物种特征（species traits）：例如，物种的生物学特性、生态习性等。
# 系统发育关系（phylogenetic relationships）：物种之间的进化关系，通常由系统发育树表示。
# 时空背景信息（spatiotemporal context）：有关采样设计的时空信息，以考虑采样单元之间的依赖性。
# 总结
# Hmsc 提供了一种强大的工具来建模和分析物种的空间分布及其与环境的关系。通过使用 JSDMs，研究人员能够更好地理解物种群落的动态，并考虑到物种之间的相互作用和环境因素的综合影响。

# 联合物种分布建模（Joint Species Distribution Modelling, JSDM）是一种强大而灵活的方法，可以用于解决多种生态学和生物多样性研究中的关键问题。以下是一些主要的应用领域和可以解决的问题：
# 1. 物种分布预测
# 预测物种的分布：JSDM可以用于预测在不同环境条件下物种的潜在分布，帮助研究人员了解在气候变化、栖息地破坏等情境下物种的未来分布模式。
# 2. 生物多样性评估
# 评估物种多样性和丰富度：通过将多种物种的分布结合在一起，JSDM可以提供对物种多样性和丰富度的综合评估，从而帮助制定保护策略。
# 3. 生物相互作用分析
# 捕捉生物相互作用：JSDM能够考虑物种之间的相互作用（如竞争、捕食和共生关系），从而提供更全面的生态网络视角。
# 4. 环境影响分析
# 评估环境因素的影响：通过将环境协变量纳入模型，研究人员可以量化不同环境因素（如土壤类型、气候条件、植被类型等）对物种分布的影响。
# 5. 处理缺失数据
# 应对缺失协变量：JSDM可以有效地处理数据中缺失的环境变量或物种信息，这在生态学研究中是一个常见的问题。
# 6. 系统发育关系的影响
# 考虑系统发育关系：通过整合系统发育信息，JSDM能够评估物种之间的亲缘关系如何影响其生态特征和分布模式。
# 7. 时空动态分析
# 研究物种的时空动态：通过加入时空背景信息，JSDM可以帮助研究人员理解物种分布的动态变化，例如季节性分布变化或长期趋势。
# 8. 保护优先级设定
# 制定保护策略：通过识别关键的物种和栖息地，JSDM可以帮助保护者制定优先级，集中资源于最需要保护的区域或物种，以有效维护生物多样性。
# 9. 影响评估
# 评估人类活动的影响：JSDM可以用于评估人类活动（如城市化、农业扩展、气候变化等）对物种分布和生态系统健康的影响，为生态管理提供科学依据。
# 通过联合物种分布建模，研究人员能够更全面地理解生态系统的复杂性和物种之间的关系，为生物多样性的保护和生态管理提供数据支持及科学指导。




