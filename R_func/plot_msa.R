#### Jinxin Meng, 20241029, 20241029 v0.1 ####
pacman::p_load(dplyr, tidyr, tibble, purrr, stringr, ggplot2, ggpubr)
library(Biostrings)
library(rlang)

#### plot_msa ####
# input a msa with fa-formated
# muscle/mafft and trimal
plot_msa <- function(path, width = 80, font_size = 2) {
  msa <- Biostrings::readAAMultipleAlignment(path) %>% 
    as.character() %>% 
    as.list() %>% 
    map2_dfc(., names(.), \(x, y) data.frame(aa = unlist(str_split(x, ""))) %>% dplyr::rename(!!y := aa)) %>% 
    rownames_to_column("site")
  
  seq_lens <- nrow(msa)
  num_seqs <- ceiling(seq_lens / width)
  total_lens <- num_seqs * width
  
  data <- map2(seq(0, width * (num_seqs - 1), width), 
               seq(width, width * num_seqs, width), \(x, y)
               filter(msa, site %in% seq(x + 1, y, 1)) %>% 
                 gather(key = "name", value = "aa", -site) %>% 
                 mutate(site = factor(site, seq(x + 1, y, 1)))) %>% 
    setNames(paste0("seq_", 1:num_seqs))
  
  if (seq_lens < total_lens) {
    data[[num_seqs]] <- rbind(data[[num_seqs]], 
                              data.frame(site = rep(seq(seq_lens + 1, total_lens, 1), each = 3), 
                                         name = rep(unique(data[[num_seqs]]$name), times = (total_lens - seq_lens))) %>%
                                add_column(aa = " ")) %>% 
      mutate(site = factor(site, seq(((num_seqs - 1) * width + 1), total_lens, 1)))
  }
  
  aa_col <- c("E" = "#ff6d6d","D" = "#ff6d6d","P" = "#f2be3c","A" = "#f2be3c","V" = "#f2be3c",
              "M" = "#f2be3c","L" = "#f2be3c","I" = "#f2be3c","G" = "#f2be3c","K" = "#769dcc",
              "R" = "#769dcc","H" = "#769dcc","N" = "#74ce98","T" = "#74ce98","C" = "#74ce98",
              "Q" = "#74ce98","S" = "#74ce98","F" = "#ffff66","Y" = "#ffff66","W" = "#ffff66",
              "-" = "#ffffff"," " = "grey95")
  
  map(data, \(x) 
      ggplot(x, aes(site, name, fill = aa)) +
        geom_tile(color = "grey", show.legend = F) +
        geom_text(aes(label = aa), size = font_size) +
        labs(x = "", y = "") +
        scale_x_discrete(breaks = c(as.character(min(as.numeric(as.character(x$site)))),
                                    as.character(max(as.numeric(as.character(x$site))))),
                         labels = c(as.character(min(as.numeric(as.character(x$site)))),
                                    as.character(max(as.numeric(as.character(x$site)))))) +
        scale_fill_manual(values = aa_col) +
        coord_fixed() +
        theme_pubr() +
        theme(axis.line = element_blank(),
              axis.text = element_text(size = 5),
              axis.ticks = element_blank())
  ) %>% cowplot::plot_grid(plotlist = ., ncol = 1)
}