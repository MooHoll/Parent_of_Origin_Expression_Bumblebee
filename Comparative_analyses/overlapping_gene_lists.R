###########################################
# Overlapping gene lists: bumblebee imprinting
###########################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Gene_lists_and_GO")

library(VennDiagram)
library(gplots)
library(rowr)
library(readr)

all <- read_csv("gene_lists_BB_imprinted.csv")
head(all)

# -----------------------------------------
# Two models: mat/pat and repro/sterile
# -----------------------------------------

imprinted_mat <- na.omit(all$sig_maternal_bias)
imprinted_pat <- na.omit(all$sig_paternal_bias)
imprinted_repro_sig <- na.omit(all$sig_repro_imprinted)
imprinted_nonrepro_sig <- na.omit(all$sig_sterile_imprinted)

list1 <- list(imprinted_mat, imprinted_nonrepro_sig, imprinted_pat, imprinted_repro_sig)
overlap <- calculate.overlap(list1)

venn.plot <- venn.diagram(
  x = list1,
  filename = NULL,
  category = c("Maternal", "Sterile", "Paternal", "Reproductive"),
  fill = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 4,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot)


# -----------------------------------------
# Imprinted and diff exp
# -----------------------------------------

diff_exp_abd <- na.omit(all$diff_exp_abd)
diff_exp_head <- na.omit(all$diff_exp_head)

# Repro and abd
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

# Repro and head
list3 <- list(imprinted_nonrepro_sig, imprinted_repro_sig, diff_exp_head)
overlap <- calculate.overlap(list3)

venn.plot3 <- venn.diagram(
  x = list3,
  filename = NULL,
  category = c("Sterile", "Reproductive", "Diff Exp Head"),
  fill = c("dodgerblue", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot3)

# Mat/Pat and Abd
list4 <- list(imprinted_mat, imprinted_pat, diff_exp_abd)
overlap <- calculate.overlap(list4)

venn.plot4 <- venn.diagram(
  x = list4,
  filename = NULL,
  category = c("Maternal", "Paternal", "Diff Exp Abd"),
  fill = c("dodgerblue", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot4)

# Mat/Pat and Head
list5 <- list(imprinted_mat, imprinted_pat, diff_exp_head)
overlap <- calculate.overlap(list5)

venn.plot5 <- venn.diagram(
  x = list5,
  filename = NULL,
  category = c("Maternal", "Paternal", "Diff Exp Head"),
  fill = c("dodgerblue", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot5)


# Mat/Pat and just repro
list6 <- list(imprinted_mat, imprinted_pat, imprinted_repro_sig)
overlap <- calculate.overlap(list6)

venn.plot6 <- venn.diagram(
  x = list6,
  filename = NULL,
  category = c("Maternal", "Paternal", "Repro"),
  fill = c("dodgerblue", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot6)


# Mat/Pat and just sterile
list7 <- list(imprinted_mat, imprinted_pat, imprinted_nonrepro_sig)
overlap <- calculate.overlap(list7)

venn.plot7 <- venn.diagram(
  x = list7,
  filename = NULL,
  category = c("Maternal", "Paternal", "Sterile"),
  fill = c("dodgerblue", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot7)


# Sterile and Repro alone
list8 <- list(imprinted_repro_sig, imprinted_nonrepro_sig)
overlap <- calculate.overlap(list8)

venn.plot8 <- venn.diagram(
  x = list8,
  filename = NULL,
  category = c("Reproductive", "Sterile"),
  fill = c("dodgerblue", "seagreen3"),
  cat.col = c("dodgerblue", "seagreen3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot8)

# Diff exp abd and head alone
list9 <- list(diff_exp_abd, diff_exp_head)
overlap <- calculate.overlap(list9)

venn.plot9 <- venn.diagram(
  x = list9,
  filename = NULL,
  category = c("Abdomen", "Head"),
  fill = c("dodgerblue", "seagreen3"),
  cat.col = c("dodgerblue", "seagreen3"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot9)


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
Genome <- 11030
  
# Maternal vs repro
p01 <- phyper(76-1, Maternal_model1, Genome-Maternal_model1, Repro_model2, lower.tail = FALSE, log.p = FALSE)
# Paternal vs repro
p02 <- phyper(57-1, Paternal_model1, Genome-Paternal_model1, Repro_model2, lower.tail = FALSE, log.p = FALSE)
# Maternal vs sterile
p03 <- phyper(71-1, Maternal_model1, Genome-Maternal_model1, Sterile_model2, lower.tail = FALSE, log.p = FALSE)
# Paternal vs sterile
p04 <- phyper(56-1, Paternal_model1, Genome-Paternal_model1, Sterile_model2, lower.tail = FALSE, log.p = FALSE)
# Maternal vs head
p05 <- phyper(6-1, Maternal_model1, Genome-Maternal_model1, Diff_exp_head, lower.tail = FALSE, log.p = FALSE)
# Paternal vs head
p06 <- phyper(3-1, Paternal_model1, Genome-Paternal_model1, Diff_exp_head, lower.tail = FALSE, log.p = FALSE)
# Maternal vs abd
p07 <- phyper(65-1, Maternal_model1, Genome-Maternal_model1, Diff_exp_abd, lower.tail = FALSE, log.p = FALSE)
# Paternal vs abd
p08 <- phyper(36-1, Paternal_model1, Genome-Paternal_model1, Diff_exp_abd, lower.tail = FALSE, log.p = FALSE)
# Repro vs abd
p09 <- phyper(100-1, Repro_model2, Genome-Repro_model2, Diff_exp_abd, lower.tail = FALSE, log.p = FALSE)
# Repro vs head
p10 <- phyper(9-1, Repro_model2, Genome-Repro_model2, Diff_exp_head, lower.tail = FALSE, log.p = FALSE)
# Repro vs sterile
p11 <- phyper(149-1, Repro_model2, Genome-Repro_model2, Sterile_model2, lower.tail = FALSE, log.p = FALSE)
# Sterile vs abd
p12 <- phyper(107-1, Sterile_model2, Genome-Sterile_model2, Diff_exp_abd, lower.tail = FALSE, log.p = FALSE)
# Sterile vs head
p13 <- phyper(11-1, Sterile_model2, Genome-Sterile_model2, Diff_exp_head, lower.tail = FALSE, log.p = FALSE)
# Abdomen vs head
p14 <- phyper(230-1, Diff_exp_abd, Genome-Diff_exp_abd, Diff_exp_head, lower.tail = FALSE, log.p = FALSE)

all_ps <- c(p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14)

# Adjust a vector of p-vales using benjimini-hochberg
adjusted_ps <- p.adjust(all_ps, method="BH")
adjusted_ps

# Looking at these results in conclusion:
# Repro and Sterile imprinted gene significantly overlap with mat/pat
# No overlap with abd diff exp sig (prob beccause so many)
# Just about sig overlap with head and maternal (not paternal), only 6 genes though so this doesn't mean much
# Sig overlap with head and repro/sterile, 9 and 11 genes though
# Sig overlap of head/abd diff exp


# -----------------------------------------
# Overlapping lists
# -----------------------------------------

ItemsList <- venn(list3, show.plot=FALSE)
list_all <- attr(ItemsList, "intersections")
list_all

# 9 genes found to be imprinted in both repro and sterile overlapping with diff exp
# 2 in just sterile overlapping with diff exp

both_overlap <- list_all$`A:B:C`
just_sterile <- list_all$`A:C`

overlap_diff_exp_head <- cbind.fill(both_overlap, just_sterile, fill =NA)
colnames(overlap_diff_exp_head) <- c("diff_exp_imprint_sterile_repro", "diff_exp_imprint_just_sterile")

write.csv(overlap_diff_exp_head, file="overlapping_imprint_diff_exp_head.csv")
