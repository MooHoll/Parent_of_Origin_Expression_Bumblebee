###----------------------------------------------
# UpSetR plots for honeybee overlap
###----------------------------------------------

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Honeybee_overlap")

library(readr)
library(UpSetR)
library(reshape2)

###----------------------------------------------
# Get data
###----------------------------------------------

# Dataframe with BB LOC and HB GB corresponding gene codes (homolog lists)
orthologs <- read_csv("reciprocal_blast_BB_HB.csv")
head(orthologs)
# Filter so only working with 6953 genes but they are reliable matches
orthologs <- subset(orthologs, orthologs$is_unique_hit == "yes")
orthologs <- subset(orthologs, orthologs$match_type == "matches_in_both_blasts")
orthologs <- orthologs[,c(1,2)]

# Gene lists from Galbriath et al. (2016) 
hb_repro_imp_direction<- read_csv("honeybee_repro_imp_direction.csv")
hb_sterile_imp_direction<- read_csv("honeybee_sterile_imp_direction.csv")
hb_diff_exp <- read_csv("honeybee_diff_exp_direction.csv")

# Relevant gene lists from our BB analysis
bumblebee_imp_lists <- read.csv("imprinted_bb_directional.csv", row.names = 1)
bb_diff_exp_lists <- read.csv("diff_exp_directional.csv", row.names = 1)

# Make one bumblebee dataframe
all_bumble <- merge(bb_diff_exp_lists, bumblebee_imp_lists, by = "gene", all=T)
all_bumble[is.na(all_bumble)] <- 0

###----------------------------------------------
# Convert bumblbee geneids to honeybee homologous gene id
###----------------------------------------------

colnames(orthologs) <- c("gene", "hb_geneid")
bumble_with_hb_id <- merge(all_bumble, orthologs, by="gene") #5544/8210 genes matched
sum(bumble_with_hb_id$Abdomen_Diff_Exp) #5498/8097
sum(bumble_with_hb_id$Abdomen_Upreg_Repro) #2753/4365
sum(bumble_with_hb_id$Abdomen_Upreg_Sterile) #2745/3732
sum(bumble_with_hb_id$Head_Diff_Exp) #119/279
sum(bumble_with_hb_id$Head_Upreg_Sterile) #45/105
sum(bumble_with_hb_id$Head_Upreg_Repro) #74/174
sum(bumble_with_hb_id$Sterile_Maternal_Bias) #39/93
sum(bumble_with_hb_id$Sterile_Paternal_Bias) #38/77
sum(bumble_with_hb_id$Reproductive_Maternal_Bias) #35/89
sum(bumble_with_hb_id$Reproductive_Paternal_Bias) #36/74
sum(bumble_with_hb_id$all_imp) #84/184
#Percentage with honeybee homolog: 47%(+/-12%)


###----------------------------------------------
# Put honeybee data in Upset format
###----------------------------------------------

hb_diff_exp <- dcast(hb_diff_exp, ID ~ `Up-regulated in`)
colnames(hb_diff_exp) <- c("hb_geneid", "hb_repro_diff_exp", "hb_sterile_diff_exp")
hb_diff_exp$hb_repro_diff_exp <- (hb_diff_exp$hb_repro_diff_exp == "Reproductive")*1
hb_diff_exp$hb_sterile_diff_exp <- (hb_diff_exp$hb_sterile_diff_exp == "Sterile")*1
hb_diff_exp$hb_diffexp <- 1

hb_repro_imp_direction <- dcast(hb_repro_imp_direction, ID ~ bias)
colnames(hb_repro_imp_direction) <- c("hb_geneid", "hb_maternal_repro", "hb_paternal_repro")
hb_repro_imp_direction$hb_maternal_repro <- (hb_repro_imp_direction$hb_maternal_repro == "maternal")*1
hb_repro_imp_direction$hb_paternal_repro <- (hb_repro_imp_direction$hb_paternal_repro == "paternal")*1

hb_sterile_imp_direction <- dcast(hb_sterile_imp_direction, ID ~ bias)
colnames(hb_sterile_imp_direction) <- c("hb_geneid", "hb_maternal_sterile", "hb_paternal_sterile")
hb_sterile_imp_direction$hb_maternal_sterile <- (hb_sterile_imp_direction$hb_maternal_sterile == "maternal")*1
hb_sterile_imp_direction$hb_paternal_sterile <- (hb_sterile_imp_direction$hb_paternal_sterile == "paternal")*1

all_honeybee <- merge(hb_sterile_imp_direction, hb_repro_imp_direction, by="hb_geneid", all=T)
all_honeybee$hb_all_imp <- 1
all_honeybee <- merge(all_honeybee, hb_diff_exp, by="hb_geneid", all=T)
all_honeybee[is.na(all_honeybee)] <- 0

###----------------------------------------------
# Identify number of Galbraith's genes in database
###----------------------------------------------

honeybee_in_db <- merge(all_honeybee, orthologs, by="hb_geneid") #1585/3049 genes matched
sum(honeybee_in_db$hb_maternal_sterile) #6/12
sum(honeybee_in_db$hb_paternal_sterile) #77/143
sum(honeybee_in_db$hb_maternal_repro) #9/13
sum(honeybee_in_db$hb_paternal_repro) #95/181
sum(honeybee_in_db$hb_all_imp) #146/273
sum(honeybee_in_db$hb_repro_diff_exp) #988/1638
sum(honeybee_in_db$hb_sterile_diff_exp) #487/1204
sum(honeybee_in_db$hb_diffexp) #1475/2842
#Percentage with honeybee in database: 54% (+/-9%)

###----------------------------------------------
# Combine HB data with bumlbe HB homologs from this study
###----------------------------------------------

overlap <- merge(all_honeybee, bumble_with_hb_id, by = "hb_geneid",all=T) #1305 genes in common (from 5544 BB and 3049 HB)
overlap[is.na(overlap)] <- 0

# Look at imprinting overlap in general (only 2 genes)
test1 <- overlap[,c(6,21)]
colnames(test1) <- c("Honeybee", "Bumblebee")

upset(test1, nsets=2,order.by = "freq",
      text.scale = 2,
      point.size = 4)
look <- subset(overlap, hb_all_imp=='1' & all_imp=='1')
# HB: LOC411889 putative serine protease K12H4.7, BB: LOC100644680 putative serine protease K12H4.7
# HB: LOC552195 uncharacterized, BB: LOC100648162 uncharacterized


# Look at diff exp overlap in general
test2 <- overlap[,c(9,11,14)]
colnames(test2) <- c("Honeybee_Ov_FB", "Bumblebee_Abd","Bumblebee_Head")

upset(test2, nsets=3,order.by = "freq",
      text.scale = 2,
      point.size = 4)


### --------------------------------------
## Hypergeometric Test for Overlap Signifiance
### --------------------------------------

# Hypergeometric test:
# Population size: 6953 homologous genes - the number of successes
# Number of successes in population: number hb imprinted genes
# Sample size: number of bb imprintied genes
# Number of successes in sample: overlapping number 

# lower-tail = T checks if below 95%
# lower-tail = F checks if under-represented compared to chance


# Overlap imprinted
sum(overlap$hb_all_imp) #273
sum(overlap$all_imp) #84
phyper(2, 273, (6953-273), 84, lower.tail = F) # p-val: 0.6469935

# Overlap diff exp with head
sum(overlap$hb_diffexp) #2842
sum(overlap$Head_Diff_Exp) #119
phyper(18, 2842, (6953-2842), 119, lower.tail = F) # p-val: 1

# Overlap diff exp with abd
sum(overlap$hb_diffexp) #2842
sum(overlap$Abdomen_Diff_Exp) #5498
phyper(1213, 2842, (6953-2842), 5498, lower.tail = F) # p-val: 1
