#!/home/mengjx/.conda/envs/mengjx/bin/Rscript
# Encoding: utf-8
# Creator: mengjx
# Created Date: 2022-7-25
# Modified Date: 2022-7-26

####option####
options(stringsAsFactors = F)
pacman::p_load(optparse,ggplot2,stringr,dplyr,tibble)
option_list <- list(make_option(c("-i", "--in"), type = "character", action = "store", help = "CAZy Profile (MAGs)"),
                    make_option(c("-p", "--prefix"), type = "character", action = "store", help = "Output Prefix"))
opt_parse <- OptionParser(option_list = option_list, usage = "Usage: %prog [options] \nDescription: This Script Used to Process CAZy Profile about MAGs!\nExample: %prog -i genome_CAZy.profile -p genome")
opt <- parse_args(opt_parse)
if (length(opt) != 3) {
   print_help(opt_parse)
   stop("At least one argument must be applied.", call. = FALSE)
}

input = opt$i
prefix = opt$p

####Initial Process####
cat(":: Enzyme Process ::\n")
profile <- read.table(input, sep = "\t", header = F, stringsAsFactors = F)
genome <- unique(profile$V1)

dt_profile <- rbind()
flag = 0
for (i in genome) {

  profile_i <- profile %>%
    subset(V1%in%i) %>% 
    select(V2) %>% 
    table(.) %>% 
    data.frame(.) %>% 
    rename(., CAZy = .)
  
  list0 <- strsplit(as.character(profile_i$CAZy), split = "\\|")
  proportion <- lapply(list0, function(x) 1/length(x)) %>% unlist()
  
  dt <- rbind()
  for (j in 1:length(list0)) {
    dt_j <- rbind()
    for (k in 1:length(list0[[j]])) {
      dt_k <- data.frame(CAZy = list0[[j]][k], freq = profile_i$Freq[j], proportion = proportion[j])
      dt_j <- rbind(dt_j, dt_k)
    }
    dt <- rbind(dt,dt_j)
  }
  
  dt$CAZy <- gsub("_\\d+", "", dt$CAZy)
  dt <- dt %>% 
    group_by(CAZy) %>% 
    mutate(res = freq * proportion) %>%
    summarise(res = sum(res)) %>%
  dt[,2] <- round(dt[,2], digits = 2) 
  colnames(dt)[2] <- i 
  
  if (flag == 0) {
    dt_profile <- rbind(dt_profile, dt)
    flag = 1
  } else {
    dt_profile <- merge(dt_profile, dt, by = "CAZy", all = T)
  }
  
}
dt_profile[is.na(dt_profile)] <- 0
write.table(dt_profile, paste0(prefix,".profile"), sep = "\t", quote = F, row.names = F)

####Module stat####
cat("::   Module Stat  ::\n")
dt_profile$Module <- str_extract_all(dt_profile$CAZy,"[A-Z]+") %>% unlist()
dt_profile <- dt_profile %>% 
  column_to_rownames(var = "CAZy") %>% 
  group_by(Module) %>% 
  summarise_all(sum)
write.table(dt_profile, paste0(prefix, "_Module.profile"), sep = "\t", quote = F, row.names = F)

####Module plot####
cat("::   Module Plot  ::\n")
plot <- dt_profile %>% 
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
ggsave(plot, filename = paste0(prefix, "_Module_All.pdf"), width = 6, height = 4.5)
cat("::     Save!      ::\n")
