##### Jinxin Meng, 20241121, 20241121, v0.1 ####

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(xcms)
library(MsExperiment)

#### parse_MChromatograms ####
cat("readSRMData() to build MChromatograms\n")
parse_MChromatograms <- function(data, label = "MChromatograms object") {
  lens <- length(data)
  parser <- map(seq_len(lens), \(x) data.frame(rtime = rtime(data[x]), intensity = intensity(data[x]))) %>% 
    setNames(seq_len(lens))
  
  metadata <- map_df(parser, \(x) data.frame(rt_min = min(x$rtime), 
                                             rt_max = max(x$rtime), 
                                             rt_intensity_max = max(pull(slice_max(x, intensity), 1)),
                                             min_intensity = quantile(x$intensity)[1],
                                             q25_intensity = quantile(x$intensity)[2],
                                             median_intensity = quantile(x$intensity)[3],
                                             q75_intensity = quantile(x$intensity)[4],
                                             max_intensity = quantile(x$intensity)[5], row.names = NULL)) %>% 
    rownames_to_column("name") %>% 
    arrange(rt_intensity_max)
  return(metadata)
  }


