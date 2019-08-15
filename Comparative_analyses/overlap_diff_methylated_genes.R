###########################################
# Overlapping gene lists: bumblebee imprinting
# Prelim with diff methylated genes from MERN
###########################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Gene_lists_and_GO")

library(VennDiagram)
library(gplots)
library(rowr)
library(readr)

all <- read_csv("gene_lists_BB_imprinted.csv")

# -----------------------------------------
# Make lists for comparison
# -----------------------------------------

imprinted_mat <- na.omit(all$sig_maternal_bias)
imprinted_pat <- na.omit(all$sig_paternal_bias)
imprinted_repro_sig <- na.omit(all$sig_repro_imprinted)
imprinted_nonrepro_sig <- na.omit(all$sig_sterile_imprinted)
diff_exp_abd <- na.omit(all$diff_exp_abd)
diff_exp_head <- na.omit(all$diff_exp_head)
diff_meth_mern <- na.omit(all$MERN_diff_meth)


# Mat/Pat and Meth
list1 <- list(imprinted_mat, imprinted_pat, diff_meth_mern)
overlap <- calculate.overlap(list1)

venn.plot1 <- venn.diagram(
  x = list1,
  filename = NULL,
  category = c("Maternal", "Paternal", "Diff Methylated"),
  fill = c("dodgerblue", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot1)

# Repro/Sterile and Meth
list2 <- list(imprinted_repro_sig, imprinted_nonrepro_sig, diff_meth_mern)
overlap <- calculate.overlap(list2)

venn.plot2 <- venn.diagram(
  x = list2,
  filename = NULL,
  category = c("Repro", "Sterile", "Diff Methylated"),
  fill = c("dodgerblue", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot2)


# Abd/Head and Meth
list3 <- list(diff_exp_abd, diff_exp_head, diff_meth_mern)
overlap <- calculate.overlap(list3)

venn.plot3 <- venn.diagram(
  x = list3,
  filename = NULL,
  category = c("Abdomen", "Head", "Diff Methylated"),
  fill = c("dodgerblue", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot3)



# -----------------------------------------
# Hypergeomatric test for comparisons
# -----------------------------------------

# ?phyper
# phyper(overlap-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE)
# population size is the number of genes in the B.terrestris genome (11030)

Maternal_model1 <- 90
Paternal_model1 <- 64
Repro_model2 <- 170
Sterile_model2 <- 163
Diff_exp_abd <- 8097
Diff_exp_head <- 279
Diff_meth <- 203
Genome <- 11030

# Maternal vs meth
p01 <- phyper(1-1, Diff_meth, Genome-Diff_meth, Maternal_model1, lower.tail = FALSE, log.p = FALSE)
# Repro vs meth
p02 <- phyper(2-1, Diff_meth, Genome-Diff_meth, Repro_model2, lower.tail = FALSE, log.p = FALSE)
# Sterile vs meth
p03 <- phyper(1-1, Diff_meth, Genome-Diff_meth, Sterile_model2, lower.tail = FALSE, log.p = FALSE)
# Abd vs meth
p04 <- phyper(170-1, Diff_meth, Genome-Diff_meth, Diff_exp_abd, lower.tail = FALSE, log.p = FALSE)
# Head vs meth
p05 <- phyper(1-1, Diff_meth, Genome-Diff_meth, Diff_exp_head, lower.tail = FALSE, log.p = FALSE)

all_ps <- c(p01,p02,p03,p04,p05)

# Adjust a vector of p-vales using benjimini-hochberg
adjusted_ps <- p.adjust(all_ps, method="BH")
adjusted_ps

# Looking at these results in conclusion:
# Only significant overlap is with the diff exp genes in abd 170/203, interesting
# as the tissue sample used for the methylation was head tissue



# -----------------------------------------
# Overlapping lists
# -----------------------------------------

ItemsList <- venn(x, show.plot=FALSE)
list_all <- attr(ItemsList, "intersections")
list_all






