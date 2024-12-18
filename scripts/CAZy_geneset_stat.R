#!/home/mengjx/.conda/envs/mengjx/bin/Rscript
# Encoding: utf-8
# Creator: mengjx
# Created Date: 2021-09-14
# Modified Date: 2022-7-30

####option####
options(stringsAsFactors = F)
pacman::p_load(optparse,ggplot2,stringr,dplyr,tibble)
option_list <- list(
  make_option(c("-i", "--in"), type = "character", action = "store", help = "CAZy Profile (Geneset)"),
  make_option(c("-p", "--prefix"), type = "character", action = "store", help = "Output Prefix"))
opt_parse <- OptionParser(option_list = option_list, usage = "Usage: %prog [options] \nDescription: This Script Used to Process CAZy Profile about Geneset!")
opt <- parse_args(opt_parse)

if (length(opt) != 2) {
  print_help(opt_parse)
  stop("At least one argument must be applied.", call. = FALSE)
}

input = opt$i
prefix = opt$p

####Initial Process####
cat("=== Enzyme Proportion & Count ===")
profile <- read.table(input, sep = "\t", header = F, row.names = NULL)

CAZy <- profile %>% 
  select(V2) %>% 
  table(.) %>% 
  data.frame(.) %>% 
  rename(., CAZy = .)
  
list0 <- strsplit(as.character(CAZy$CAZy), split = "\\|")
proportion <- lapply(list0, function(x) 1/length(x)) %>% unlist()

dt <- rbind()
for (i in 1:length(list0)) {
  dt_i <- rbind()
  for (j in 1:length(list0[[i]])) {
    dt_j <- data.frame(CAZy = list0[[i]][j], freq = CAZy$Freq[i], proportion = proportion[i])
    dt_i <- rbind(dt_i, dt_j)
  }
  dt <- rbind(dt,dt_i)
}

dt$CAZy <- gsub("_\\d+", "", dt$CAZy)
dt <- dt %>% 
  group_by(CAZy) %>% 
  mutate(res = freq * proportion) %>% 
  summarise(res = sum(res)) %>% 
  mutate(res = round(res, digits = 2))

write.table(dt, paste0(prefix,".profile"), sep = "\t", quote = F, row.names = F)

####Module stat####
cat("=== Module Stat ===")
dt$Module <- str_extract_all(dt$CAZy,"[A-Z]+") %>% 
  unlist()
dt <- dt %>% 
  column_to_rownames(var = "CAZy") %>% 
  group_by(Module) %>% 
  summarise_all(sum)
write.table(dt, paste0(prefix, "_Module.profile"), sep = "\t", quote = F, row.names = F)

####Module plot####
cat("=== Module Plot ===\n")
plot <- dt %>% 
  column_to_rownames(var = "Module") %>% 
  t() %>% 
  colSums() %>% 
  data.frame() %>% 
  rownames_to_column(var = "Module") %>% 
  ggplot(aes(x = Module, y = ., fill = Module, color = Module)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = .), color = "black", vjust = 1.1) +
  scale_fill_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462")) +
  scale_color_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462")) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Count", fill = "Module", color = "Module") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))
ggsave(plot, filename = paste0(prefix, "_All_Module.pdf"), width = 6, height = 4.5)
