###########################################
# ~2/3 imprinted genes diff exp in abdomen tissue 
# Not significant overlap but there are a huge num of genes diff exp in abd
###########################################

# From the script overlapping_gene_lists.R

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Gene_lists_and_GO")

library(VennDiagram)
library(gplots)
library(rowr)
library(readr)
library(doBy)
library(tidyr)

all <- read_csv("gene_lists_BB_imprinted.csv")

imprinted_repro_sig <- na.omit(all$sig_repro_imprinted)
imprinted_nonrepro_sig <- na.omit(all$sig_sterile_imprinted)
diff_exp_abd <- na.omit(all$diff_exp_abd)

list2 <- list(imprinted_nonrepro_sig, imprinted_repro_sig, diff_exp_abd)
overlap <- calculate.overlap(list2)

venn.plot2 <- venn.diagram(
  x = list2,
  filename = NULL,
  category = c("Sterile", "Reproductive", "Diff Exp Abd"),
  fill = c("dodgerblue", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot2)

# Of these genes found to be diff exp, which direction does the imprinting go?
# Need to pull out the lists of 80, 18 and 11

# 89 genes
diff_exp_both_sterile_repro_imp <- as.data.frame(overlap[[1]])
colnames(diff_exp_both_sterile_repro_imp) <- "Gene"
# 18 genes
diff_exp_sterile_imp <- as.data.frame(overlap[[3]])
colnames(diff_exp_sterile_imp) <- "Gene"
# 11 genes
diff_exp_repro_imp <- as.data.frame(overlap[[4]])
colnames(diff_exp_repro_imp) <- "Gene"


# Get the file of the imprinted genes which gives if they're paternally or maternally
input <- read_csv("stats_imprinting_model2_total counts per gene reproductive nonreproductive workers.csv")

signbias_repr = (input$signinbothdirs_repr==1)
signbias_nonrepr = (input$signinbothdirs_nonrepr==1)

# 93 genes 
sterile_mat_imp_genes = as.data.frame(unique(input$geneID[signbias_nonrepr&
                                                     (input$avgpropmatexpr>0.6)]))
colnames(sterile_mat_imp_genes) <- "Gene"

# 77 genes 
sterile_pat_imp_genes = as.data.frame(unique(input$geneID[signbias_nonrepr&
                                                      (input$avgpropmatexpr<0.4)]))
colnames(sterile_pat_imp_genes) <- "Gene"

# 89 genes 
repro_mat_imp_genes = as.data.frame(unique(input$geneID[signbias_repr&
                                                  (input$avgpropmatexpr>0.6)]))
colnames(repro_mat_imp_genes) <- "Gene"

# 74 genes 
repro_pat_imp_genes = as.data.frame(unique(input$geneID[signbias_repr&
                                                  (input$avgpropmatexpr<0.4)]))
colnames(repro_pat_imp_genes) <- "Gene"

# Find out if the diff exp genes are upreg in repro or sterile
ABD_upreg_in_nonrepro <- read_csv("ABD_upreg_in_nonrepro.csv")
colnames(ABD_upreg_in_nonrepro)[1]<- "Gene"
ABD_upreg_in_repro <- read_csv("ABD_upreg_in_repro.csv")
colnames(ABD_upreg_in_repro)[1]<- "Gene"



# Get the imprint status and exp status of the 18 genes
sterile_imp_upreg_sterile <- merge(diff_exp_sterile_imp, ABD_upreg_in_nonrepro, by = "Gene") #15/18
sterile_imp_upreg_sterile <- as.data.frame(sterile_imp_upreg_sterile[,1])
sterile_imp_upreg_sterile <- na.omit(sterile_imp_upreg_sterile[,1])
sterile_imp_upreg_repro <- merge(diff_exp_sterile_imp, ABD_upreg_in_repro, by = "Gene") #3/18
sterile_imp_upreg_repro <- as.data.frame(sterile_imp_upreg_repro[,1])
sterile_imp_upreg_repro <- na.omit(sterile_imp_upreg_repro[,1])

sterile_imp_mat_exp <- merge(diff_exp_sterile_imp, sterile_mat_imp_genes, by = "Gene") #13/18
sterile_imp_mat_exp <- na.omit(sterile_imp_mat_exp[,1])
sterile_imp_pat_exp <- merge(diff_exp_sterile_imp, sterile_pat_imp_genes, by = "Gene") #5/18
sterile_imp_pat_exp <- na.omit(sterile_imp_pat_exp[,1])

list1 <- list(sterile_imp_upreg_sterile, sterile_imp_upreg_repro, sterile_imp_mat_exp,sterile_imp_pat_exp)
overlap <- calculate.overlap(list1)

venn.plot1 <- venn.diagram(
  x = list1,
  filename = NULL,
  category = c("UpReg Sterile", "UpReg Repro", "Maternal Exp", "Paternal Exp"),
  fill = c("dodgerblue", "seagreen3", "orchid3","lightgoldenrod2"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3","lightgoldenrod2"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3,
  main = "Diff Exp and Sterile Imprinted Genes"
)
grid.newpage()
grid.draw(venn.plot1)

# Get the imprint status and exp status of the 11 genes
repro_imp_upreg_repro <- merge(diff_exp_repro_imp, ABD_upreg_in_repro, by = "Gene") #6/11
repro_imp_upreg_repro <- as.data.frame(repro_imp_upreg_repro[,1])
repro_imp_upreg_repro <- na.omit(repro_imp_upreg_repro[,1])
repro_imp_upreg_sterile <- merge(diff_exp_repro_imp, ABD_upreg_in_nonrepro, by = "Gene") #5/11
repro_imp_upreg_sterile <- as.data.frame(repro_imp_upreg_sterile[,1])
repro_imp_upreg_sterile <- na.omit(repro_imp_upreg_sterile[,1])

repro_imp_mat_exp <- merge(diff_exp_repro_imp, repro_mat_imp_genes, by = "Gene") #7/11
repro_imp_mat_exp <- na.omit(repro_imp_mat_exp[,1])
repro_imp_pat_exp <- merge(diff_exp_repro_imp, repro_pat_imp_genes, by = "Gene") #4/11
repro_imp_pat_exp <- na.omit(repro_imp_pat_exp[,1])


list3 <- list(repro_imp_upreg_sterile, repro_imp_upreg_repro, repro_imp_mat_exp,repro_imp_pat_exp)
overlap <- calculate.overlap(list3)

venn.plot3 <- venn.diagram(
  x = list3,
  filename = NULL,
  category = c("UpReg Sterile", "UpReg Repro", "Maternal Exp", "Paternal Exp"),
  fill = c("dodgerblue", "seagreen3", "orchid3","lightgoldenrod2"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3","lightgoldenrod2"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3,
  main = "Diff Exp and Repro Imprinted Genes"
)
grid.newpage()
grid.draw(venn.plot3)





# Get the imprint status and exp status of the 89 genes
both_imp_upreg_repro <- merge(diff_exp_both_sterile_repro_imp,ABD_upreg_in_repro, by="Gene") #29 / 89
both_imp_upreg_repro <- as.data.frame(both_imp_upreg_repro[,1])
both_imp_upreg_repro <- na.omit(both_imp_upreg_repro[,1])
both_imp_upreg_sterile <- merge(ABD_upreg_in_nonrepro,diff_exp_both_sterile_repro_imp, by="Gene") #60 /89
both_imp_upreg_sterile <- as.data.frame(both_imp_upreg_sterile[,1])
both_imp_upreg_sterile <- na.omit(both_imp_upreg_sterile[,1])


both_imp_mat_exp_repro <- merge(repro_mat_imp_genes,diff_exp_both_sterile_repro_imp, by = "Gene") #55 /89
both_imp_mat_exp_repro <- na.omit(both_imp_mat_exp_repro[,1])
both_imp_pat_exp_repro <- merge(repro_pat_imp_genes,diff_exp_both_sterile_repro_imp, by = "Gene") #34/89
both_imp_pat_exp_repro <- na.omit(both_imp_pat_exp_repro[,1])

both_imp_mat_exp_sterile <- merge(sterile_mat_imp_genes,diff_exp_both_sterile_repro_imp, by = "Gene") #55/ 89
both_imp_mat_exp_sterile <- na.omit(both_imp_mat_exp_sterile[,1])
both_imp_pat_exp_sterile <- merge(sterile_pat_imp_genes,diff_exp_both_sterile_repro_imp, by = "Gene") #34/89
both_imp_pat_exp_sterile <- na.omit(both_imp_pat_exp_sterile[,1])


check_1 <- merge(both_imp_mat_exp_repro, both_imp_mat_exp_sterile, by="Gene") #55
check_1 <- na.omit(check_1[,1])
check_2 <- merge(both_imp_pat_exp_repro, both_imp_pat_exp_sterile, by="Gene") #34
check_2 <- na.omit(check_2[,1])


list4 <- list(both_imp_upreg_repro, both_imp_upreg_sterile, check_1,check_2)
overlap <- calculate.overlap(list4)

venn.plot4 <- venn.diagram(
  x = list4,
  filename = NULL,
  category = c("UpReg Repro", "UpReg Sterile", "Maternal Exp", "Paternal Exp"),
  fill = c("dodgerblue", "seagreen3", "orchid3","lightgoldenrod2"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3","lightgoldenrod2"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3,
  main = "Diff Exp and Both Imprinted Genes"
)
grid.newpage()
grid.draw(venn.plot4)









list4 <- list(both_imp_upreg_repro, both_imp_upreg_sterile, both_imp_mat_exp_repro,both_imp_pat_exp_repro,both_imp_mat_exp_sterile,both_imp_pat_exp_sterile)
overlap <- calculate.overlap(list4)

venn.plot4 <- venn.diagram(
  x = list3,
  filename = NULL,
  category = c("UpReg Repro", "UpReg Sterile", "Maternal Exp Repro", "Paternal Exp Repro", "Maternal Exp Sterile", "Paternal Exp Sterile"),
  fill = c("dodgerblue", "seagreen3", "orchid3","lightgoldenrod2", "snow4","wheat3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3","lightgoldenrod2", "snow4","wheat3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3,
  main = "Diff Exp and Both Imprinted Genes"
)
grid.newpage()
grid.draw(venn.plot4)



