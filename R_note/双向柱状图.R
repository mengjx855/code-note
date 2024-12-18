####info####
# Encoding: utf-8
# Creator: 小梦
# Modified date：20222-5-28

####Option####
setwd("")
options(stringsAsFactors = F)
pacman::p_load(dplyr,tidyr,tibble,ggplot2)

####code####
dat <- read.delim("dat.txt", sep = "\t", row.names = NULL)
dat$item <- factor(dat$item, levels = unique(dat$item))
dat <- dat %>% mutate(percent = event/n*100)
dat$percent <- round(dat$percent, digits = 2)
dat$label[!is.nan(dat$percent)] <- paste0(dat$percent[!is.nan(dat$percent)], "%")

ggplot(dat) +
  geom_bar(aes(x = n, y = item, fill = species), stat = "identity") +
  scale_fill_manual(values = c("#d8b365","#5ab4ac")) +
  geom_bar(aes(x = event, y = item, fill = species), stat = "identity", fill = "grey99", alpha = .6) +
  geom_vline(xintercept = 0, size = .2) +
  scale_x_continuous(breaks = seq(-4500, 3000, by = 1500), 
                     labels = as.character(abs(seq(-4500, 3000, by = 1500))), 
                     limits = c(-5500, 3500)) +
  geom_text(data = subset(dat, species=="Cow"), aes(x = -5400, y = item, label = label), size = 2) +
  geom_text(data = subset(dat, species=="Cattle"), aes(x = 3400, y = item, label = label), size = 2) +
  labs(x = "Number of Samples", y = "Type", fill = "Species") +
  theme_bw() +
  theme(line = element_line(color = "black"),
        axis.text = element_text(color = "black"), 
        axis.ticks = element_line(size = .2))
ggsave(filename = "plot.pdf", width = 6, height = 4.5)
