# Getting base GO set for all methylated genes and all genes identified in the RNA-seq dataset

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/GO_analysis_NEW")

library(readr)

GO_annotations <- read.table("Bumble_bee_ensemble_GO_terms.txt")
colnames(GO_annotations) <- c("geneID","goID")

# ----------------------------------------------------------------  
abd_genes <- read_csv("ABD_all_genes_results.csv") 
abd_genes <- as.data.frame(unique(abd_genes$X1)) #9791
colnames(abd_genes) <- "geneID"

head_genes <- read_csv("HEAD_all_genes_results.csv") 
head_genes <- as.data.frame(unique(head_genes$X1)) #9791
colnames(head_genes) <- "geneID"

exp_gos <- merge(GO_annotations, abd_genes, by="geneID")
length(exp_gos[is.na(exp_gos$geneID)]) # 0
write.table(exp_gos, file = "GOs_all_expression_genes_as_background.txt", sep="\t",
            col.names = T, row.names = F, quote = F)

# ----------------------------------------------------------------  

abd_diff_exp_genes <- read_csv("gene_lists/LOC_abd_diff_exp_genes.txt", 
                                   col_names = FALSE) #8097
colnames(abd_diff_exp_genes) <- "geneID"

head_diff_exp_genes <- read_csv("gene_lists/LOC_head_diff_exp_genes.txt", 
                               col_names = FALSE) #279
colnames(head_diff_exp_genes) <- "geneID"


abd_gos <- merge(GO_annotations, abd_diff_exp_genes, by="geneID")
length(abd_gos[is.na(abd_gos$geneID)]) # 0
write.table(abd_gos, file = "GOs_all_adb_diffexp_genes_as_background.txt", sep="\t",
            col.names = T, row.names = F, quote = F)

head_gos <- merge(GO_annotations, head_diff_exp_genes, by="geneID")
length(head_gos[is.na(head_gos$geneID)]) # 0
write.table(head_gos, file = "GOs_all_head_diffexp_genes_as_background.txt", sep="\t",
            col.names = T, row.names = F, quote = F)

# ---------------------------------------------------------------- 

PoO_genes <- read_csv("PoO_genes_all.csv", 
                        col_names = FALSE) #184
colnames(PoO_genes) <- "geneID"

PoO_gos <- merge(GO_annotations, PoO_genes, by="geneID")
length(PoO_gos[is.na(PoO_gos$geneID)]) # 0
write.table(PoO_gos, file = "GOs_all_PoO_genes_as_background.txt", sep="\t",
            col.names = T, row.names = F, quote = F)
