# Jinxin Meng, 20241023, 20241023 -------------

# 由于reactome没有对应的层级结构，获取某个ID的上级关系需要用此函数进行确定
reactome_longest_path <- function(path_ID, organism = NULL) {
  path_file <- "/database/Reactome/ReactomePathwaysRelation.txt"
  path_info_file <- "/database/Reactome/ReactomePathways.txt"
  org <- list("BTA" = "Bos taurus", "GGA" = "Gallus gallus", "HSA" = "Homo sapiens", 
              "MMU" = "Mus musculus", "SSC" = "Sus scrofa")
  
  organism <- ifelse(is.null(organism), unlist(strsplit(path_ID, "-"))[2], organism)
  path <- read.delim(path_file, header = F, col.names = c("from", "to")) %>% 
    dplyr::filter(grepl(organism, to)) %>% 
    igraph::graph_from_data_frame(directed = T)
  path_info <- read.delim(path_info_file, header = F, col.names = c("id", "name", "tax")) %>% 
    filter(tax == org[[organism]])
  
  # 查找从所有节点出发的所有路径, 遍历每个节点作为起点
  all_paths <- purrr::map(igraph::V(path), \(x) igraph::all_simple_paths(path, from = x)) %>% 
    purrr::list_flatten()

  paths_with_target <- base::Filter(\(x) path_ID %in% names(igraph::V(path)[x]), all_paths)
  longest_path <- paths_with_target[[which.max(map_vec(paths_with_target, \(x) length(x)))]]
  out <- data.frame(id = (names(igraph::V(path))[longest_path])) %>% 
    mutate(name = path_info$name[match(id, path_info$id)]) %>% 
    add_column(seq_id = seq_len(nrow(.)), .before = 1)
  
  return(out)
}
