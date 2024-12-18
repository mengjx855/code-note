library(ggridges)
library(ggplot2)

ggplot(diamonds, aes(x = price, y = cut, fill = cut)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

# library
library(ggridges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Load dataset from github
data <- read.table("https://raw.githubusercontent.com/zonination/perceptions/master/probly.csv", header=TRUE, sep=",")
data <- data %>% 
  gather(key="text", value="value") %>%
  mutate(text = gsub("\\.", " ",text)) %>%
  mutate(value = round(as.numeric(value),0)) %>%
  filter(text %in% c("Almost Certainly","Very Good Chance","We Believe","Likely","About Even", "Little Chance", "Chances Are Slight", "Almost No Chance"))

# Plot
data %>%
  mutate(text = fct_reorder(text, value)) %>%
  ggplot(aes(y=text, x=value,  fill=text)) +
  geom_density_ridges(alpha=0.6, stat="binline", bins=20) +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("") +
  ylab("Assigned Probability (%)")

library(ggridges)
library(ggplot2)

ggdensity(data, "value", fill = "group", palette = c("#6893ca", "#d0614f"), alpha = .1) +
  facet_wrap(~factor(time), ncol = 1)

ggplot(data, aes(value, y = time, fill = group)) +
  geom_density_ridges(scale = 0.9, alpha = 0.7) +
  scale_fill_manual(values = c("#6893ca", "#d0614f")) +
  theme_ridges() +
  theme(axis.title.y = element_blank()) +
  labs(x = "Scale abundance", y = "Visit time") +
  theme(axis.text.y = element_text(size = 12),
    legend.position = "bottom")
