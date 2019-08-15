##########################################
# Looking for common imprinted genes HB and BB
##########################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Honeybee")

library(readr)
library(tidyverse)
library(rowr)
library(VennDiagram)
library(stringr)

# Dataframe with BB LOC and HB GB corresponding gene codes
orthologs <- read_csv("reciprocal_blast_BB_HB.csv")
head(orthologs)
# Filter so only working with 6953 genes but they are reliable matches
orthologs <- subset(orthologs, orthologs$is_unique_hit == "yes")
orthologs <- subset(orthologs, orthologs$match_type == "matches_in_both_blasts")

# Gene lists from Galbriath et al. (2016) NEED TO GET MATERNAL/PATERNAL IMPRINTED LISTS
honeybee_gene_lists <- read_csv("honeybee_genes_galbraith.csv")
head(honeybee_gene_lists)
length(na.omit(honeybee_gene_lists$Diff_exp_repro_sterile)) #2842 genes diff exp
length(na.omit(honeybee_gene_lists$Repro_imprinted)) #201 repro imrprinted
length(na.omit(honeybee_gene_lists$Sterile_imprinted)) #164 sterile imprinted
length(na.omit(honeybee_gene_lists$Diff_exp_repro_imprint_overlap))#49 repro imprinted and diff exp
length(na.omit(honeybee_gene_lists$Diff_exp_sterile_imprint_overlap)) #40 sterile imrprint and diff exp

# Relevant gene lists from our BB analysis
bumblebee_gene_lists <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Gene_lists_and_GO/gene_lists_BB_imprinted.csv")
head(bumblebee_gene_lists)
length(na.omit(bumblebee_gene_lists$all_imprinted)) #850 all imprinted
length(na.omit(bumblebee_gene_lists$sig_maternal_bias)) #90 maternal bias
length(na.omit(bumblebee_gene_lists$sig_paternal_bias)) #64 paternal bias
length(na.omit(bumblebee_gene_lists$sig_sterile_imprinted)) #170 sterile imprinted
length(na.omit(bumblebee_gene_lists$sig_repro_imprinted)) #163 repro imprinted
length(na.omit(bumblebee_gene_lists$diff_exp_abd)) #8097 diff exp abdomen
length(na.omit(bumblebee_gene_lists$diff_exp_head)) #279 diff exp head

# ------------------------------------------
# Looking for repro/sterile imprinted overlap
# ------------------------------------------
repro_look <- cbind.fill(honeybee_gene_lists$Repro_imprinted, honeybee_gene_lists$Sterile_imprinted,
                         bumblebee_gene_lists$sig_repro_imprinted, bumblebee_gene_lists$sig_sterile_imprinted,
                         fill = NA)
repro_look <- repro_look[rowSums(is.na(repro_look)) != ncol(repro_look), ]
colnames(repro_look) <- c("HB_repro", "HB_sterile", "BB_repro", "BB_sterile")
head(repro_look)

# Convert the bumblbee gene IDs to the matching HB ortholog so can compare
sterile_merge <- as.data.frame(repro_look$BB_sterile)
colnames(sterile_merge) <- "BB"
sterile_merge_both <- merge(sterile_merge, orthologs, by = "BB") #77/170 with orthologs
sterile_BB_orthologes <- as.data.frame(sterile_merge_both$HB)
colnames(sterile_BB_orthologes) <- "sterile_BB_orthologes"

repro_merge <- as.data.frame(repro_look$BB_repro)
colnames(repro_merge) <- "BB"
repro_merge_both <- merge(repro_merge, orthologs, by = "BB") #71/163 with orthologs
repro_BB_orthologes <- as.data.frame(repro_merge_both$HB)
colnames(repro_BB_orthologes) <- "repro_BB_orthologes"

honeybee_repro_imprint <- cbind.fill(sterile_BB_orthologes, repro_BB_orthologes,
                                     repro_look$HB_repro, repro_look$HB_sterile,
                                     fill = NA)
colnames(honeybee_repro_imprint)[3:4] <- c("HB_repro", "HB_sterile")
honeybee_repro_imprint <- as.data.frame(apply(honeybee_repro_imprint,
                                              2,function(x)gsub('\\s+', '',x)))

# For venn diagram need each overlapping dataset as a vector
sterile_BB_orthologes <- na.omit(honeybee_repro_imprint$sterile_BB_orthologes)
repro_BB_orthologes <- na.omit(honeybee_repro_imprint$repro_BB_orthologes)
HB_sterile <- na.omit(honeybee_repro_imprint$HB_sterile)
HB_repro <- na.omit(honeybee_repro_imprint$HB_repro)

venn.plot <- venn.diagram(
  x = list(sterile_BB_orthologes,repro_BB_orthologes,HB_sterile,HB_repro),
  filename = NULL,
  category = c("A.mel Orthologs of B.ter Sterile", 
               "A.mel Orthologs of B.ter Reproductive", 
               "A.mel Steile", "A.mel Reproductive"),
  fill = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 1,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot)

# Genes overlapping:
list1234 <- list(sterile_BB_orthologes,repro_BB_orthologes,HB_sterile,HB_repro)
cal1 <- calculate.overlap(list1234)
cal1 # "GB43942" "GB46687"


# ------------------------------------------
# Look for overlap in the opposite direction of ortholog search incase
# ------------------------------------------

# Try the other way, HB IDs to matching BB ortholog to see if makes a diff
sterile_merge <- as.data.frame(repro_look$HB_sterile)
colnames(sterile_merge) <- "HB"
sterile_merge_both <- merge(sterile_merge, orthologs, by = "HB") #89/164 with orthologs
sterile_HB_orthologes <- as.data.frame(sterile_merge_both$BB)
colnames(sterile_HB_orthologes) <- "sterile_HB_orthologes"

repro_merge <- as.data.frame(repro_look$HB_repro)
colnames(repro_merge) <- "HB"
repro_merge_both <- merge(repro_merge, orthologs, by = "HB") #106/201 with orthologs
repro_HB_orthologes <- as.data.frame(repro_merge_both$BB)
colnames(repro_HB_orthologes) <- "repro_HB_orthologes"

bumblebee_repro_imprint <- cbind.fill(sterile_HB_orthologes, repro_HB_orthologes,
                                     repro_look$BB_repro, repro_look$BB_sterile,
                                     fill = NA)
colnames(bumblebee_repro_imprint)[3:4] <- c("BB_repro", "BB_sterile")
bumblebee_repro_imprint <- as.data.frame(apply(bumblebee_repro_imprint,
                                              2,function(x)gsub('\\s+', '',x)))


# For venn diagram need each overlapping dataset as a vector
sterile_HB_orthologes <- na.omit(bumblebee_repro_imprint$sterile_HB_orthologes)
repro_HB_orthologes <- na.omit(bumblebee_repro_imprint$repro_HB_orthologes)
BB_repro <- na.omit(bumblebee_repro_imprint$BB_repro)
BB_sterile <- na.omit(bumblebee_repro_imprint$BB_sterile)

venn.plot1 <- venn.diagram(
  x = list(sterile_HB_orthologes,repro_HB_orthologes,BB_repro,BB_sterile),
  filename = NULL,
  category = c("B.ter Orthologs of A.mel Sterile", 
               "B.ter Orthologs of A.mel Reproductive", 
               "B.ter Steile", "B.ter Reproductive"),
  fill = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 1,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot1)

# Genes overlapping:
list1234 <- list(sterile_HB_orthologes,repro_HB_orthologes,BB_repro,BB_sterile)
cal1 <- calculate.overlap(list1234)
cal1 # "LOC100644680" "LOC100648162"


# Checked the ortholog file and these genes from both overlaps are orthologs
# "LOC100644680" putative serine protease K12H4.7
# "GB43942" putative serine protease K12H4.7

#"LOC100648162" uncharacterized
#"GB46687" uncharacterized
