#### GSEA by Jinxin Meng

all_data <- read.csv("F:/Rproject/Rproj_OpDatabase/OpDatabase/result.csv")
#gene.sybmol <- rownames(all_data)
#drugs <- colnames(all_data)
#data.exp[is.na(data.exp)] <- 0
#data.exp <- data.exp + 0.00000001

immune_genes <- read.csv("F:/MineR/immune_rps.csv")
immune_genes <-as.character(immune_genes[,2])
matching_genes <- intersect(immune_genes, all_data[,1])
all_data <- all_data[all_data[,1] %in% as.character(matching_genes),]
gene.sybmol <- all_data[,1]
sampe_name <- colnames(all_data)[-1]
all_data <- all_data[,2:32]

all_data <- data.frame(all_data)

A <- apply(all_data, 2, function(x){order(x, decreasing = TRUE)})

immune_sets <- read.csv("F:/MineR/immune_geneset.gmt", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
set_ids <- immune_sets[,1]
set_names <- immune_sets[,2]
immune_sets <- immune_sets[,-(1:2)]


upES <- matrix(nrow = length(immune_sets[,1])/2, ncol = length(A[1,]))
downES <- matrix(nrow = length(immune_sets[,1])/2, ncol = length(A[1,]))
TES <- matrix(nrow = length(immune_sets[,1])/2, ncol = length(A[1,]))
score_pvalue <- matrix(nrow = length(immune_sets[,1])/2, ncol = length(A[1,]))

for(i in 1:length(A[1,])) {
  gene.list2 <- A[, i]
  #gene.value <- sort(all.data.exp[, i], decreasing = TRUE)
  # calculate background score distribution between gene profile of the drug
  # and a random selection of ranked genes -> (null dist for p-values)
  nTrials = 1e3
  deg_set_len = 250
  # print (Sys.time())
  
  matching_order <- match(matching_genes, gene.sybmol)
  rSc = rand_cmap_scores(drug_sig=gene.list2, m_genes=matching_order,
                         deg_set_len=deg_set_len, nTrials=nTrials)
  
  m=0
  for(j in seq(1, length(immune_sets[,1])-1, 2)){
    imup_gene_order <- unique(match(as.character(immune_sets[j,]), gene.sybmol))
    imup_gene_order <- imup_gene_order[!is.na(imup_gene_order)]
    imdw_gene_order <- unique(match(as.character(immune_sets[j+1,]), gene.sybmol))
    imdw_gene_order <- imdw_gene_order[!is.na(imdw_gene_order)]
    m=m+1
    
    upES[m, i] <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = imup_gene_order, weighted.score.type = 0))
    downES[m, i] <- unlist(GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = imdw_gene_order, weighted.score.type = 0))
    
    if (sign(upES[m, i]) == sign(downES[m, i])) {
      TES[m, i] <- 0
    } else {
      TES[m, i] <- upES[m, i] - downES[m, i]
    }
    
    if (TES[m, i] != 0) {
      score_pvalue[m, i] = length(which(abs(rSc) >= abs(TES[m, i])))  / nTrials
    }
  }
}

colnames(TES) <- sampe_name
rownames(TES) <- set_names[seq(1,607,2)]
write.csv(TES,"immunologicalfeaturesfordrug_TES.csv")

colnames(upES) <- sampe_name
rownames(upES) <- set_names[seq(1,607,2)]
write.csv(upES,"immunologicalfeaturesfordrug_upES.csv")

colnames(downES) <- sampe_name
rownames(downES) <- set_names[seq(2,608,2)]
write.csv(downES,"immunologicalfeaturesfordrug_downES.csv")

colnames(score_pvalue) <- sampe_name
rownames(score_pvalue) <- set_names[seq(1,607,2)]
write.csv(score_pvalue,"immunologicalfeaturesfordrug_score_pvalue.csv")
# gsea caculation, infact the caculation from the ref, Discovery of drug mode of action and drug repositioning from transcriptional responses

# gsea caculation, infact the caculation from the ref, Discovery of drug mode of action and drug repositioning from transcriptional responses
