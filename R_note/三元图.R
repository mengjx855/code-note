# Jinxin Meng, 20240801, 20240814 ---------------------

metadata <- read.delim("metadata.txt")
profile <- read.delim("profile.txt", row.names = 1, check.names = F) %>% 
  profile_transRA()
group <- read.delim("sample_group.txt")

plot_data <- profile_smp2grp_2(profile, group, method = "mean") %>% 
  rowwise() %>% 
  mutate(avg_ab = (CD + UC + nonIBD) / 3)

library(ggtern)
ggtern(plot_data, aes(CD, UC, nonIBD)) + 
  geom_point(aes(size = avg_ab), color = "#fc8d62") +
  scale_size_continuous(range = c(1, 6)) +
  theme_bw() +  
  theme_minimal()