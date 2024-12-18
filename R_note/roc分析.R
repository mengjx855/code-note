# roc analysis --------

# anlaysis --------
# 建模型获得预测值
# pred
# sample  group_a group_b
# A0002	0.3704	0.6296
# A0005	0.1418	0.8582
# A0006	0.1184	0.8816
# A0010	0.1282	0.8718
pacman::p_load(randomForest, pROC)
rf_model <- randomForest(group ~ ., data = train, ntree = 1000, importance = F, proximity = T)
pred <- data.frame(predict(rf_model, test, type = 'prob'))
pred$group <- group$group[match(pred$sample, group$sample)] # 加一列真实的分组
roc <- roc(pred$group, pred$Disease) # roc分析

# plot ------------------
# 如果plot单个roc的结果
auc = round(roc$auc, digits = 4)
ggroc(roc, legacy.axes = T, color = "#377eb8") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black", lty = 2, lwd = .2) +
  labs(color = "", x = "1 - Specificity", y = "Sensitivity") +
  annotate("text" ,x = 0.75, y = 0.125 , label = auc) +
  theme_bw() +
  theme(aspect.ratio = 1)

# 如果plot多个roc的结果
# 各个roc要放到一个list
roc_list <- list(roc_a = roc_a, roc_b = roc_b)
ggroc(roc_l, legacy.axes = T, lwd = .35) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black", lty = 2, lwd = .4) +
  scale_color_manual(values = colors, breaks = study, labels = labels) +
  labs(color = "", x = "1 - Specificity", y = "Sensitivity") +
  theme_bw() +
  theme(aspect.ratio = 1)

# 多个roc的结果，弄一个label出来，展示auc的结果
labels <- lapply(roc_list, \(x) round(x$auc, digits = 4)) %>% 
  unlist() %>% as.numeric() %>% paste0(names(roc_list), ": ", .)
# 结果是这样的
roc_a: 0.901

# 计算auc置信区间
ci.auc(x) # 三个值: auc下限，auc，auc上限
# roc_list可以这么算
labels <- lapply(roc_list, \(x) paste0(round(x$auc, digits = 3), "\n(95% Cl: ", paste(round(ci.auc(x), digits = 3)[c(1,3)], collapse = " ~ "), ")")) %>% 
  unlist() %>% paste0(names(roc_list), " AUC: ", .)

# 置信区间展示在图中
# 如果是一个roc的话
roc_se <- ci.se(roc, specificities = seq(0, 1, 0.01), conf.level = 0.95) %>% 
      data.frame(check.names = F) %>% 
      dplyr::rename(lower = all_of("2.5%"), upper = all_of("97.5%"), median = all_of("50%")) %>% 
      rownames_to_column(var = "spec") %>% 
      mutate(spec = as.numeric(spec))
# 加在图中
p <- p + geom_ribbon(data = roc_se, aes(x = 1 - spec, ymin = lower, ymax = upper), fill = "#238443", alpha = .1)
# 如果是roc_list
roc_se <- lapply(roc_list, \(x) ci.se(x, specificities = seq(0, 1, 0.01), conf.level = 0.95) %>% 
                       data.frame(check.names = F) %>% 
                       rename(lower = all_of("2.5%"), upper = all_of("97.5%"), median = all_of("50%")) %>% 
                       rownames_to_column(var = "spec") %>% 
                       mutate(spec = as.numeric(spec))) %>% 
      purrr::map2_dfr(., names(.), \(x, y) add_column(x, class = y))
# 加在图中
p <- p + geom_ribbon(data = roc_se, aes(x = 1 - spec, ymin = lower, ymax = upper, fill = class), alpha = .1, inherit.aes = F, show.legend = F) +
      scale_fill_manual(values = colors)
