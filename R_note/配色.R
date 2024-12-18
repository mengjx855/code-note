# https://zhuanlan.zhihu.com/p/593320758
colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(19)
# 12-color set3
colors <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")
colors <- c("#80b1d3","#b3de69","#fdb462","#8dd3c7","#bc80bd","#fb8072","#ffed6f","#fccde5","#bebada","#ccebc5","#ffffb3","#d9d9d9")
# 12-color paired（浅深）
colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")
# 12-color paired（深浅）
colors <- c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#ff7f00","#fdbf6f","#6a3d9a","#cab2d6","#b15928","#ffff99")
# 20-color paired (深浅)
color <- paletteer::paletteer_d("ggthemes::Classic_20")
# 20-scale_color_hue()
color <- c("#f8766d","#ea8331","#d89000","#c09b00","#a3a500","#7cae00","#39b600","#00bb4e","#00bf7d","#00c1a3","#00bfc4","#00bae0","#00b0f6","#35a2ff","#9590ff","#c77cff","#e76bf3","#fa62db","#ff62bc","#ff6a98")
# 14-color
colors = c("#373735", "#f4e4cb", "#bec3e3", "#e6d2d3", "#f5e9bf", "#faece1", "#8994b4", "#5b5b5d", "#f9bc7d", "#cdc9e0", "#d58b8a", "#bbac85", "#764701", "#5e9bd2")

# 2case-1control
breaks <- c("ctrl", "trt1", "trt2"),
labels <- c("Control", "Treatment 1", "Treatment 2")
colors <- c("#999999", "#E69F00", "#56B4E9")

# case-control
levels <- c("Disease", "Control")
colors <- c("#f46d43","#74add1")  # 多种免疫病病毒组研究中病例和对照组用的颜色
colors <- c("#f77d4d", "#1f78b4")  # MDD,IBS,ACVD病例对照中的颜色方案
colors <- c("#3e4faa", "#ea3726")
colors <- c("#f37f73". "#80b1d3")

# control baseline 2wk 4wk 24wk
# Gout病毒组病例研究的配色方案
colors <- c("#1f78b4","#ff7f00","#33a02c","#e31a1c","#6a3d9a")
labels <- c("Control", "baseline", "2wk", "4wk", "24wk")



# 病毒组分析 物种和宿主信息填充颜色
colors <- c("#4E9F50FF", # Actinobacteriota
            "#EF8A0CFF", # Bacteroidota
            "#FCC66DFF", # Desulfobacterota_A
            "#C3CE3DFF", # Eremiobacterota
            "#98D9E4FF", # Firmicutes
            "#94A323FF", # Fusobacteriota
            "#F7D42AFF", # Multiple
            "#3CA8BCFF", # Proteobacteria
            "#8DBFA8FF", # Unknown 
            "#26897EFF")  # Verrucomicrobiota

colors <- c("#B15928", # Autographiviridae
            "#9c964a", # Circoviridae
            "#8DD3C7FF", # Flandersviridae
            "#FFFFB3FF", # Gratiaviridae
            "#BEBADAFF", # Inoviridae
            "#FB8072FF", # Microviridae
            "#80B1D3FF", # Myoviridae
            "#CCEBC5FF", # Podoviridae
            "#FDB462FF", # Podoviridae_crAss-like
            "#B3DE69FF", # Quimbyviridae
            "#FFED6FFF", # Retroviridae
            "#BCBD22FF", # Salasmaviridae
            "#17BECFFF", # Schitoviridae
            "#BC80BDFF", # Siphoviridae
            "grey85") # Unclassified


