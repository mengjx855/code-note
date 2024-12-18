df <- data.frame(
  kingdom = c("Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria"),
  phylum = c("Firmicutes", "Firmicutes", "Firmicutes", "Proteobacteria", "Proteobacteria"),
  class = c("Bacilli", "Bacilli", "Clostridia", "Gammaproteobacteria", "Alphaproteobacteria"),
  order = c("Lactobacillales", "Bacillales", "Clostridiales", "Enterobacterales", "Rhizobiales"),
  family = c("Lactobacillaceae", "Bacillaceae", "Clostridiaceae", "Enterobacteriaceae", "Rhizobiaceae"),
  genus = c("Lactobacillus", "Bacillus", "Clostridium", "Escherichia", "Rhizobium"),
  value = c(100, 200, 150, 300, 250)
)


df_long <- df %>%
  gather(key = "level", value = "name", -value) %>%
  arrange(level, name) %>%
  group_by(level, name) %>%
  summarise(value = sum(value)) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(
    ymax = cumsum(value),
    ymin = lag(ymax, default = 0)
  )


df_long <- df_long %>%
  mutate(level = factor(level, levels = c("kingdom", "phylum", "class", "order", "family", "genus")))

ggplot() +
  geom_rect(data = df_long %>% filter(level == "kingdom"), aes(fill = name, ymax = ymax, ymin = ymin, xmax = 6, xmin = 5), color = "white") +
  ggnewscale::new_scale_fill() +
  geom_rect(data = df_long %>% filter(level == "phylum"), aes(fill = name, ymax = ymax, ymin = ymin, xmax = 5, xmin = 4), color = "white") +
  ggnewscale::new_scale_fill() +
  geom_rect(data = df_long %>% filter(level == "class"), aes(fill = name, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3), color = "white") +
  ggnewscale::new_scale_fill() +
  geom_rect(data = df_long %>% filter(level == "order"), aes(fill = name, ymax = ymax, ymin = ymin, xmax = 3, xmin = 2), color = "white") +
  ggnewscale::new_scale_fill() +
  geom_rect(data = df_long %>% filter(level == "family"), aes(fill = name, ymax = ymax, ymin = ymin, xmax = 2, xmin = 1), color = "white") +
  ggnewscale::new_scale_fill() +
  geom_rect(data = df_long %>% filter(level == "genus"), aes(fill = name, ymax = ymax, ymin = ymin, xmax = 1, xmin = 0), color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "none")