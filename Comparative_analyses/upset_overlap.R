###----------------------------------------------
# UpSetR plots (better than venn diagrams)
###----------------------------------------------

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Final_imprinting_model/Imprinting_graphs")
library(readr)
library(UpSetR)
library(rowr)

###----------------------------------------------
# Get Imprinting data
###----------------------------------------------


# Data from model 2
model2 <- read_csv("stats_imprinting_model2_total counts per gene reproductive nonreproductive workers.csv")

signbias_repr <- (model2$signinbothdirs_repr==1)
signbias_nonrepr <- (model2$signinbothdirs_nonrepr==1)

#170 overall significantly nonrepro imprinted
verysig_nonrepro <- as.data.frame(unique(model2$geneID[signbias_nonrepr&
                                                         (model2$avgpropmatexpr>0.6|model2$avgpropmatexpr<0.4)]))
#163 overall significantly repro imprinted
very_sig_repro <- as.data.frame(unique(model2$geneID[signbias_repr&
                                                       (model2$avgpropmatexpr>0.6|model2$avgpropmatexpr<0.4)]))

#repro_imp <- model2[model2$geneID %in% very_sig_repro$`unique(model2$geneID[signbias_repr & (model2$avgpropmatexpr > 0.6 | model2$avgpropmatexpr < 0.4)])`,]
#sterile_imp <- model2[model2$geneID %in% verysig_nonrepro$`unique(model2$geneID[signbias_nonrepr & (model2$avgpropmatexpr > 0.6 | model2$avgpropmatexpr < 0.4)])`,]
#write.csv(repro_imp, file="repro_imprinted_genes.csv")
#write.csv(sterile_imp, file="sterile_imprinted_genes.csv")

# 93 nonrepro maternal
maternal_nonrepro <- as.data.frame(unique(model2$geneID[signbias_nonrepr&(model2$avgpropmatexpr>0.6)]))
colnames(maternal_nonrepro)<-"maternal_nonrepro"
# 77 nonrepro paternal
paternal_nonrepro <- as.data.frame(unique(model2$geneID[signbias_nonrepr&(model2$avgpropmatexpr<0.4)]))
colnames(paternal_nonrepro)<-"paternal_nonrepro"

# 89 repro maternal
maternal_repro <- as.data.frame(unique(model2$geneID[signbias_repr&(model2$avgpropmatexpr>0.6)]))
colnames(maternal_repro)<-"maternal_repro"
# 74 repro paternal
paternal_repro <- as.data.frame(unique(model2$geneID[signbias_repr&(model2$avgpropmatexpr<0.4)]))
colnames(paternal_repro)<-"paternal_repro"


haigs_test_genes <- cbind.fill(maternal_nonrepro,paternal_nonrepro,maternal_repro,paternal_repro,
                               fill = NA)

maternal_nonrepro <- na.omit(haigs_test_genes$maternal_nonrepro)
maternal_repro <- na.omit(haigs_test_genes$maternal_repro)
paternal_nonrepro <- na.omit(haigs_test_genes$paternal_nonrepro)
paternal_repro <- na.omit(haigs_test_genes$paternal_repro)

###----------------------------------------------
# Upset: make matrix and plot: overlap within imprinted gens
###----------------------------------------------


# Get full list of genes for making the right input for upsetR
colnames(very_sig_repro) <- "gene"
colnames(verysig_nonrepro) <- "gene"

full_list_imp <- rbind(very_sig_repro, verysig_nonrepro)
full_list_imp <- as.data.frame(full_list_imp[!duplicated(full_list_imp),])
colnames(full_list_imp) <- "gene"

# Add new columns with 1 to show presence of that gene in the corresponding category
full_list_imp$Sterile_Maternal_Bias = (full_list_imp$gene %in% maternal_nonrepro)*1
full_list_imp$Reproductive_Maternal_Bias = (full_list_imp$gene %in% maternal_repro)*1
full_list_imp$Sterile_Paternal_Bias = (full_list_imp$gene %in% paternal_nonrepro)*1
full_list_imp$Reproductive_Paternal_Bias = (full_list_imp$gene %in% paternal_repro)*1
full_list_imp$all_imp = 1
#write.csv(full_list_imp, file="imprinted_bb_directional.csv")

full_list_imp1 <- full_list_imp[,-6]
upset(full_list_imp1, order.by = "freq",
      text.scale = 2,
      point.size = 4)


###----------------------------------------------
# Get Differential Exp Data
###----------------------------------------------
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Differential_expression")
exp_files <- list.files(path="~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Differential_expression/", pattern="LOC")


# Ascdataframes to can do the rbind to get a concensus 
for (i in 1:length(exp_files)) assign(exp_files[i], read_csv(exp_files[i], col_names = 'gene'))
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/")

###----------------------------------------------
# Upset: make matrix and plot: imprinted and diff exp
###----------------------------------------------


full_list_diff_exp <- rbind(LOC_abd_diff_exp_genes.txt, LOC_abd_upreg_nonrepro.txt,
                            LOC_abd_upreg_repro.txt, LOC_head_diff_exp_genes.txt,
                            LOC_head_upreg_nonrepro.txt, LOC_head_upreg_repro.txt)
full_list_diff_exp <- as.data.frame(full_list_diff_exp[!duplicated(full_list_diff_exp),])
colnames(full_list_diff_exp) <- "gene"


# Read in same data as vectors now to get the logical output for upset (messy to read in twice I know)
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Differential_expression")
for (i in 1:length(exp_files)) assign(exp_files[i], read.csv(exp_files[i], header = F)$V1)
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/")


# Add new columns with 1 to show presence of that gene in the corresponding category
full_list_diff_exp$Abdomen_Diff_Exp = (full_list_diff_exp$gene %in% LOC_abd_diff_exp_genes.txt)*1
full_list_diff_exp$Abdomen_Upreg_Repro = (full_list_diff_exp$gene %in% LOC_abd_upreg_nonrepro.txt)*1
full_list_diff_exp$Abdomen_Upreg_Sterile = (full_list_diff_exp$gene %in% LOC_abd_upreg_repro.txt)*1
full_list_diff_exp$Head_Diff_Exp = (full_list_diff_exp$gene %in% LOC_head_diff_exp_genes.txt)*1
full_list_diff_exp$Head_Upreg_Sterile = (full_list_diff_exp$gene %in% LOC_head_upreg_nonrepro.txt)*1
full_list_diff_exp$Head_Upreg_Repro = (full_list_diff_exp$gene %in% LOC_head_upreg_repro.txt)*1
write.csv(full_list_diff_exp, file="diff_exp_directional.csv")

upset(full_list_diff_exp, nsets=6,order.by = "freq",
      text.scale = 2,
      point.size = 4)

diff_exp_and_imp <- merge(full_list_imp, full_list_diff_exp, by="gene", all=T)
diff_exp_and_imp[is.na(diff_exp_and_imp)] <- 0

# Look at diff exp in abdomen and all imprinting
test <- diff_exp_and_imp[,c(6,7)]

upset(test, nsets=2,order.by = "freq",
      text.scale = 2,
      point.size = 4)


# Look at diff exp in abdomen and all imprinting
test0 <- diff_exp_and_imp[,c(6,10)]

upset(test0, nsets=2,order.by = "freq",
      text.scale = 2,
      point.size = 4)


# Look at diff exp in abdomen and directional imprinting
test1 <- diff_exp_and_imp[,c(1:5,7)]

upset(test1, nsets=6,order.by = "freq",
      text.scale = 2,
      point.size = 4)


# Look at diff exp in head and directional imprinting
test2 <- diff_exp_and_imp[,c(1:5,10)]

upset(test2, nsets=6,order.by = "freq",
      text.scale = 2,
      point.size = 4)


# Upreg sterile/repro and directional imprinting: head
test4 <- diff_exp_and_imp[,c(1:5,11,12)]

upset(test4, nsets=7,order.by = "freq",
      text.scale = 2,
      point.size = 4)

# Upreg sterile/repro and directional imprinting: abd
test3 <- diff_exp_and_imp[,c(1:5,8,9)]

upset(test3, nsets=7,order.by = "freq",
      text.scale = 2,
      point.size = 4)





# NOTE ABOUT UPSET OVERLAPS:
# Overlapping repro/sterile maternal imprinted gives 79, but when include 
# diff exp in the same upset plot the overlap between just repro/sterile maternal imp
# decreases to 24, this is because the 55 not included here also overlap with 
# the diff exp genes so they aren't included in the new overlap 


### --------------------------------------
## Hypergeometric Test for Overlap Signifiance
### --------------------------------------

# Hypergeometric test:
# Population size: 11030 genes in genome (grep gff file) - the number of successes
# Number of successes in population: number imprinted genes
# Sample size: number of diff expressed genes
# Number of successes in sample: overlapping number 

# lower-tail = T checks if below 95%
# lower-tail = F checks if under-represented compared to chance

### --------------------------------------
## Overlap Imprinted Genes
### --------------------------------------

# Overlap of maternal imprinted genes in repro/sterile
phyper(79, 163, 10867, 170, lower.tail = F) # p-val: 9.200035e-108***

# Overlap of paternal imprinted genes in repro/sterile
phyper(70, 163, 10867, 170, lower.tail = F) # p-val: 7.665742e-90***

### --------------------------------------
## Overlap Imprinted Genes with Head Diff Exp
### --------------------------------------

# Overlap of all imprinted genes in with head diff exp (184 = unique imprinted genes)
phyper(11, 184, 10846, 279, lower.tail = F) # p-val: 0.002483512***

# Unique maternal imprinted genes
union(maternal_nonrepro, maternal_repro) #103
# Overlap of head differential expression with maternal bias 
# 7 total overlap, 6 overlap both repro/sterile
phyper(7, 103, 10927, 279, lower.tail = F) # p-val: 0.004528954***
# Look at genes
look <- subset(diff_exp_and_imp, Head_Diff_Exp=='1' & all_imp=='1')
# Papilin: a serine protease inhibitor comes up: sterile/repro paternal bias, upreg in repro head and abdomen
# it's highly expressed in embryonic stages of drosophilia 
# Serine Protease Inhibitor 3/4: maternal bias in both, only diff exp in head, upreg in repro, 
# AKA sepins, associated with protease inhibition but also immune response and as a storage in eggs

# Unique paternal imprinted genes
union(paternal_nonrepro, paternal_repro) #81
# Overlap of head differential expression with paternal bias 
phyper(3, 81, 10949, 279, lower.tail = F) # p-val: 0.1488569

# Overlap of head differential expression with sterile maternal bias (0 are only repro maternal bias)
phyper(1, 93, 10937, 279, lower.tail = F) # p-val: 0.6861941


### --------------------------------------
## Overlap Imprinted Genes with Abd Diff Exp
### --------------------------------------

# Overlap of all imprinted genes in with abd diff exp (184 = unique imprinted genes)
phyper(118, 184, 10846, 8097, lower.tail = F) # p-val: 0.996705

# Overlap of abd differential expression with maternal bias 
phyper(75, 103, 10927, 8097, lower.tail = F) # p-val: 0.5167988

# Overlap of abd differential expression with paternal bias 
phyper(43, 81, 10949, 8097, lower.tail = F) # p-val: 0.9999344


### --------------------------------------
## Pull out the lists of diff exp and imp overlap
### --------------------------------------

imp_diff_exp_head <- subset(diff_exp_and_imp, Head_Diff_Exp=='1' & all_imp=='1')
imp_diff_exp_head_genes <- imp_diff_exp_head$gene
write.csv(imp_diff_exp_head_genes, file="genes_overlapping_imp_diffexp_head.csv")

imp_diff_exp_abd <- subset(diff_exp_and_imp, Abdomen_Diff_Exp=='1' & all_imp=='1')
imp_diff_exp_abd_genes <- imp_diff_exp_abd$gene
write.csv(imp_diff_exp_abd_genes, file="genes_overlapping_imp_diffexp_abd.csv")

all_gos_with_description <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Final_imprinting_model/GO_analysis/all_gos_with_description.csv")

Bumble_bee_ensemble_GO_terms <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Final_imprinting_model/GO_analysis/Bumble_bee_ensemble_GO_terms.txt", 
                                           "\t", escape_double = FALSE, col_names = FALSE, 
                                           trim_ws = TRUE)
colnames(Bumble_bee_ensemble_GO_terms) <- c("gene", "go_id")

imp_diff_exp_head_GO <- merge(imp_diff_exp_head, Bumble_bee_ensemble_GO_terms, by="gene")
imp_diff_exp_abd_GO <- merge(imp_diff_exp_abd, Bumble_bee_ensemble_GO_terms, by="gene")

imp_diff_exp_head_GO <- merge(imp_diff_exp_head_GO, all_gos_with_description, by = "go_id")
imp_diff_exp_head_GO <- unique(imp_diff_exp_head_GO)
write.csv(imp_diff_exp_head_GO, file = "imprinted_diff_exp_head_GOterms.csv")

imp_diff_exp_abd_GO <- merge(imp_diff_exp_abd_GO, all_gos_with_description, by = "go_id")
imp_diff_exp_abd_GO <- unique(imp_diff_exp_abd_GO)
write.csv(imp_diff_exp_abd_GO, file = "imprinted_diff_exp_abd_GOterms.csv")
