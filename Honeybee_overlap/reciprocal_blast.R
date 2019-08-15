##########################################
# Making bumblbee and honeybee homologs
##########################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Honeybee")


library(readr)
library(tidyverse)


honeybee_genome_to_bb <- read_delim("honeybee_genome_to_bb.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE,
                                    col_names = F)
colnames(honeybee_genome_to_bb) <- c("hb_geneID", "bb_geneID", "percent_identical_matches",
                                     "alignment_length", "num_mismatches", "num_gaps", 
                                     "alignment_start_hb", "alignment_end_hb", 
                                     "alignment_start_bb", "alignment_end_bb", "evalue",
                                     "bitscore")

bumblebee_genome_to_hb <- read_delim("bumblebee_genome_to_hb.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE,
                                     col_names = F)
colnames(bumblebee_genome_to_hb) <- c("bb_geneID", "hb_geneID", "percent_identical_matches",
                                     "alignment_length", "num_mismatches", "num_gaps", 
                                     "alignment_start_bb", "alignment_end_bb", 
                                     "alignment_start_hb", "alignment_end_hb", "evalue",
                                     "bitscore")

# Genes appear more than once in the blast results as multiple parts of the gene
# align, so checked that it's only the best match, i.e. BB gene1 has only one matching 
# HB gene, not more than one, then remove duplicate rows so left with list of BB/HB genes 
# and their matches only once
n_occur <- data.frame(table(bumblebee_genome_to_hb$bb_geneID))
n_occur[n_occur$Freq > 1,]
bumblebee_genome_to_hb[bumblebee_genome_to_hb$bb_geneID %in% n_occur$Var1[n_occur$Freq > 1],]

have_a_look<-subset(bumblebee_genome_to_hb, bumblebee_genome_to_hb$bb_geneID=="LOC100631054")

bb_to_hb_no_duplicates<-bumblebee_genome_to_hb [!duplicated(bumblebee_genome_to_hb[c(1,2)]),]
n_occur_no_dups <- data.frame(table(bb_to_hb_no_duplicates$bb_geneID))
n_occur_no_dups[n_occur_no_dups$Freq > 1,]

# and for honeybee ...
n_occur <- data.frame(table(honeybee_genome_to_bb$hb_geneID))
n_occur[n_occur$Freq > 1,]
honeybee_genome_to_bb[honeybee_genome_to_bb$hb_geneID %in% n_occur$Var1[n_occur$Freq > 1],]

have_a_look<-subset(honeybee_genome_to_bb, honeybee_genome_to_bb$hb_geneID=="GB40007")

hb_to_bb_no_duplicates<-honeybee_genome_to_bb [!duplicated(honeybee_genome_to_bb[c(1,2)]),]
n_occur_no_dups <- data.frame(table(bb_to_hb_no_duplicates$hb_geneID))
n_occur_no_dups[n_occur_no_dups$Freq > 1,]


# Now check the BB-HB best match is the same as the HB-BB best match 
# NOTE: not all genes had a match between the species so there are a handful of 
# missing genes from both lists

bumblbee_matches <- as.data.frame(paste(bb_to_hb_no_duplicates$bb_geneID,
                                        bb_to_hb_no_duplicates$hb_geneID))
colnames(bumblbee_matches) <- c("bb_list")
bumblbee_matches$bb_list<-gsub(" ", "", bumblbee_matches$bb_list)


honeybee_matches <- as.data.frame(paste(hb_to_bb_no_duplicates$bb_geneID,
                                        hb_to_bb_no_duplicates$hb_geneID))
colnames(honeybee_matches) <- c("hb_list")
honeybee_matches$hb_list<-gsub(" ", "", honeybee_matches$hb_list)

# Most of the BB-HB best matches are the same as the HB-BB best matches
honeybee_matches$bb_matches <- bumblbee_matches$bb_list[match(honeybee_matches$hb_list, bumblbee_matches$bb_list)]
length(honeybee_matches$bb_matches[!is.na(honeybee_matches$bb_matches)]) #7345/7973


# Start to put a file together: matches in both blasts, total genes:7345
final_data <- honeybee_matches
final_data <- as.data.frame(final_data[,-1])
final_data <- as.data.frame(final_data[!is.na(final_data)])
colnames(final_data) <- c("match_in_both_blasts")

final_data$BB <- gsub("GB.+$", "", final_data$match_in_both_blasts)
final_data$HB <- gsub("^.+GB", "GB", final_data$match_in_both_blasts)
final_data <- final_data[,-1]
final_data$match_type <- "matches_in_both_blasts"


# Start to put a file together: matches in only honeybee blast, total genes:1899
honeybee_data <- honeybee_matches
honeybee_data <- honeybee_data[is.na(honeybee_data$bb_matches),]

honeybee_data$BB <- gsub("GB.+$", "", honeybee_data$hb_list)
honeybee_data$HB <- gsub("^.+GB", "GB", honeybee_data$hb_list)
honeybee_data <- honeybee_data[,-c(1,2)]
honeybee_data$match_type <- "matches_in_only_honeybee"
nrow(honeybee_data) # 1899
final_data <- rbind(final_data, honeybee_data)


# Start to put a file together: matches in only bumblebee blast, total genes: 628
bumblbee_matches$hb_matches <- honeybee_matches$hb_list[match(bumblbee_matches$bb_list, honeybee_matches$hb_list)]
length(bumblbee_matches$hb_matches[is.na(bumblbee_matches$hb_matches)]) #628

bumblebee_data <- bumblbee_matches
bumblebee_data <- bumblebee_data[is.na(bumblebee_data$hb_matches),]

bumblebee_data$BB <- gsub("GB.+$", "", bumblebee_data$bb_list)
bumblebee_data$HB <- gsub("^.+GB", "GB", bumblebee_data$bb_list)
bumblebee_data <- bumblebee_data[,-c(1,2)]
bumblebee_data$match_type <- "matches_in_only_bumblebee"
final_data <- rbind(final_data, bumblebee_data)


# Now check if there were matches made differently between the unique HB and BB matches
one_match <- subset(final_data, !final_data$match_type=="matches_in_both_blasts")

# More than one HB gene for one BB gene
n_occur <- data.frame(table(one_match$BB))
homo_to_more_than_one_HB_gene<-n_occur[n_occur$Freq > 1,]
nrow(homo_to_more_than_one_HB_gene) #409
subset(one_match, one_match$BB=="LOC100652167")

# More than one BB gene for one HB gene
n_occur <- data.frame(table(one_match$HB))
homo_to_more_than_one_BB_gene<-n_occur[n_occur$Freq > 1,]
nrow(homo_to_more_than_one_BB_gene) #203
subset(one_match, one_match$HB=="GB55943")

# Label up these offending genes in the final dataframe to make sure informed decisions
# are made in future analysis (then keep only the one's that occur in both matches if
# more are present in individual blasts)
test <- final_data
nrow(test) #9872

colnames(homo_to_more_than_one_HB_gene) <- c("BB", "Num_gene_matches")
colnames(homo_to_more_than_one_BB_gene) <- c("HB", "Num_gene_matches")

test_HB <- merge(test, homo_to_more_than_one_BB_gene, by ="HB")
test_HB$is_unique_hit <- "no"
nrow(test_HB) #442

test_BB <- merge(test, homo_to_more_than_one_HB_gene, by ="BB")
test_BB$is_unique_hit <- "no"
nrow(test_BB) #1469

new_test_bb <- test_BB[,c(1,2)]
new_test_hb <- test_HB[,c(1,2)]
common <- intersect(new_test_bb, new_test_hb) 
nrow(common) #134 in common

# Make a dataframe with the rows not accounted for above
# Should be 9872-((442+1469)-134)=8095
test_unique <- test[!(test$BB %in% test_BB$BB),] 
test_unique_2 <- test_unique[!(test_unique$HB %in% test_HB$HB),] 
nrow(test_unique_2)#8095
test_unique_2$Num_gene_matches <- 1
test_unique_2$is_unique_hit <- "yes"

very_final_data <- rbind(test_unique_2, test_BB, test_HB)
very_final_data <- very_final_data[!duplicated(very_final_data),]
nrow(very_final_data) #9958 (86 more than final_data above, weird, but doesn't change the
# core set to be used which is the 8095 above, unlikly to use the others, and if they
# are to be used they will be assessed per gene)

write.csv(very_final_data, file = "reciprocal_blast_BB_HB.csv")



