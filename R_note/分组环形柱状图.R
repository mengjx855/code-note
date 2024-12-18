#### info ####
# Jinxin Meng, 20231216, 20231216
pacman::p_load(tidyr, dplyr, tibble, purrr)
pacman::p_load(ggplot2)

#### 有gap的环形柱状图 ####
# 逻辑上大概就是添加几个空的x轴变量，空开距离然后作图的时候使用极坐标转。
data <- rbind(data.frame(individual = paste( "Mister ", seq(1, 60), sep = ""), value = sample(seq(10, 100), 60, replace = T)),
              data.frame(individual = rep(NA, 10), value = rep(NA, 10))) %>% 
  add_column(id = seq_len(nrow(.)))

# 发散状的label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

ggplot(data, aes(x=as.factor(id), y=value)) +
  geom_bar(stat = "identity", fill = "#66c2a5") +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar(start = 0) + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

#### 分组环形柱状图 ####
# Create dataset
data <- data.frame(
  individual=paste( "Mister ", seq(1,60), sep=""),
  group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
  value=sample( seq(10,100), 60, replace=T)
)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p


# circular barplot
col_phy <- readRDS("../bacterial_genome/r.phylum_color.rds")
taxa <- read.delim("../bacterial_genome/dat.strains.metadata.txt")
BGCs <- read.delim("Network_Annotations_Full.tsv", check.names = F)
dat <- BGCs %>% 
  select(BGC, class = "BiG-SCAPE class") %>% 
  filter(!grepl("BGC", BGC)) %>% 
  mutate(genome = ifelse(grepl("BGC", x = BGC), "", stringr::str_extract(BGC, "[A-Za-z0-9.]+"))) %>% 
  left_join(taxa %>% select(name, phylum), by = c("genome" = "name")) %>% 
  mutate(phy = ifelse(phylum %in% names(col_phy), phylum, "p__Other")) %>% 
  select(phy, class) %>% 
  group_by(phy, class) %>% 
  summarise(n = n()) %>% 
  mutate(n_val = log10(n), 
         class = factor(class),
         n_lab = prettyNum(n, big.mark = ","))

# plot_dat
empty_bar <- 2
empty_df <- data.frame(matrix(NA, empty_bar*nlevels(dat$class), ncol(dat))) %>% 
  rename(any_of(structure(colnames(.), names = colnames(dat))))
empty_df$class <- rep(levels(dat$class), each = empty_bar)
plot_dat <- rbind(dat, empty_df) %>% 
  arrange(class, phy) %>% 
  add_column(id = 1:nrow(.))

# plot_lab
plot_lab <- plot_dat %>% 
  mutate(cols = nrow(.),
         angle = 90 - 360 * (id - 0.5) / cols,
         hjust = ifelse(angle < -90, 1, 0),
         angle = ifelse(angle < -90, angle + 180, angle))

plot_title <- plot_dat %>% 
  group_by(class) %>% 
  summarize(min = min(id), max = max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(mean = mean(c(min, max)))

plot_grid <- plot_dat %>% 
  filter(is.na(n)) %>% 
  select(class, id) %>% 
  group_by(class) %>% 
  summarise(min = min(id), max = max(id))

ggplot(plot_dat, aes(x = as.factor(id), y = n_val, fill = class)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  ylim(-3, max(na.omit(plot_dat$n_val)) + 1) +
  coord_polar() +
  theme_void() +
  geom_segment(data = plot_grid, aes(x = min, y = 1, xend = max, yend = 1), color = "grey", alpha = 1, size = 0.3, inherit.aes = F) +
  geom_segment(data = plot_grid, aes(x = min, y = 2, xend = max, yend = 2), color = "grey", alpha = 1, size = 0.3, inherit.aes = F) +
  geom_segment(data = plot_grid, aes(x = min, y = 3, xend = max, yend = 3), color = "grey", alpha = 1, size = 0.3, inherit.aes = F) +
  annotate("text", x = rep(.5, 3), y = c(1, 2, 3), label = c("1", "2", "3"), color = "grey", size = 3, angle = 0, fontface = "bold", hjust = 1) +
  geom_text(data = plot_lab, aes(x = id, y = n_val + .2, label = sub("p__", "", phy), hjust = hjust, angle = angle), color = "black", fontface = "bold", alpha = 0.6, size = 2.5, inherit.aes = F) +
  geom_text(data = plot_lab, aes(x = id, y = n_val - .1, label = n_lab, hjust = hjust, angle = angle), color="black", fontface = "bold", alpha = 0.6, size = 2.5, inherit.aes = F) +
  geom_segment(data = plot_title, aes(x = min, y = -.3, xend = max, yend = -.3), color = "black", alpha = 0.8, size = 0.8 , inherit.aes = F) +
  geom_text(data = plot_title, aes(x = mean, y = -.5, label = class, color = class), hjust = .5, vjust = .5, alpha = 0.8, size = 4, fontface = "bold", inherit.aes = F) +
  scale_fill_manual(values = c("#58BBD0","#EB9256","#9B99CB","#68A7BE","#AD9F2A","#87CCBA","#FEDFB2","#EE7E77")) +
  scale_color_manual(values = c("#58BBD0","#EB9256","#9B99CB","#68A7BE","#AD9F2A","#87CCBA","#FEDFB2","#EE7E77"))

