#!/home/mengjx/.conda/envs/mengjx/bin/Rscript
# Encoding: utf-8 -*-
# mengjx 09/08/2021

pacman::p_load(argparse,dplyr,tibble)

parser <- ArgumentParser()
parser$add_argument("-i", help = "input cazy profile")
args <- parser$parse_args()
input = args$i

# 整理profile
print("modify cazy profile!")
print("load cazy profile!")

cazy_profile <- read.table(input, sep = "\t", header = T, stringsAsFactors = F)
sample_name <- gsub(".cazy","",colnames(cazy_profile)[-1])
cazy_name <- cazy_profile[,1]
red_data <- cazy_name[grep("_",cazy_name)]
cazy_name_dered <- gsub("_[0-9]*","",cazy_name)
cazy_profile$name <- cazy_name_dered
colnames(cazy_profile)[-1] <- sample_name

cazy_profile_val <- cazy_profile[,-1]
cazy_profile_val_float <- matrix(nrow = nrow(cazy_profile_val), ncol = ncol(cazy_profile_val))
for (i in 1:nrow(cazy_profile_val)) {
      for (j in 1:ncol(cazy_profile_val)) {
              var1 <- eval(parse(text = cazy_profile_val[i,j]))
    cazy_profile_val_float[i,j] <- var1
      }
}
cazy_profile_val_float <- round(cazy_profile_val_float, digits = 4) 
cazy_profile_val_float <- rownames_to_column(as.data.frame(cazy_profile_val_float))
cazy_profile_val_float$rowname <- cazy_name_dered
colnames(cazy_profile_val_float)[-1] <- sample_name

cazy_profile_val_float_group <- cazy_profile_val_float %>% group_by(rowname) %>% summarise_all(sum)
print("save result profile")
write.table(cazy_profile_val_float_group, file = paste0(input,".fix"), sep = "\t", quote = F, row.names = F, col.names = T)
