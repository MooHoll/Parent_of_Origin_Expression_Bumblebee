###########################################
# Looking for maternal and paternal bias in repro/sterile imprinted
###########################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Gene_lists_and_GO")

library(readr)
library(stringr)
library(VennDiagram)
library(rowr)
library(gplots)

# Data from model 1
#model1 <- read_csv("stats_imprinting_model1_total_counts_per_gene.csv")
#model1 <- subset(model1, model1$signimprinbothcrosses==1)
#signmatbias <- (model1$signimprinbothcrosses==1&model1$avgpropmatexpr>0.6)
#signpatbias <- (model1$signimprinbothcrosses==1&model1$avgpropmatexpr<0.4)

#maternal_bias <- as.data.frame(unique(input$geneID[signmatbias]))
#paternal_bias <- as.data.frame(unique(input$geneID[signpatbias]))


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


#### These gene lists may hold the answer to Haig's kinship theory test

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

write.csv(haigs_test_genes, file="./Haigs_test/haigs_test_genes.csv")

maternal_nonrepro <- na.omit(haigs_test_genes$maternal_nonrepro)
maternal_repro <- na.omit(haigs_test_genes$maternal_repro)
paternal_nonrepro <- na.omit(haigs_test_genes$paternal_nonrepro)
paternal_repro <- na.omit(haigs_test_genes$paternal_repro)

list2 <- list(maternal_nonrepro, maternal_repro, 
              paternal_nonrepro, paternal_repro)
overlap <- calculate.overlap(list2)

venn.plot2 <- venn.diagram(
  x = list2,
  filename = NULL,
  category = c("Maternal Repro", "Maternal Sterile", "Paternal Repro", "Paternal Sterile"),
  fill = c("dodgerblue", "seagreen3", "orchid3", "darkorange1"),
  cat.col = c("dodgerblue", "seagreen3", "orchid3","darkorange1"),
  cat.cex = 2,
  margin = 0.05,
  cex = 3,
  fontface ='bold',
  cat.fontface ="bold"
)
grid.newpage()
grid.draw(venn.plot2)


ItemsList <- venn(list2, show.plot=FALSE)
list_all <- attr(ItemsList, "intersections")
list_all

maternal_both <- list_all$`A:B`
paternal_both <- list_all$`C:D`
maternal_just_repro <- list_all$A
paternal_just_repro <- list_all$C
maternal_just_sterile <- list_all$B
paternal_just_sterile <- list_all$D

overlapping_unique <- cbind.fill(maternal_both, paternal_both, maternal_just_repro, paternal_just_repro,
                                 maternal_just_sterile, paternal_just_sterile, fill =NA)
colnames(overlapping_unique) <- c("maternal_both", "paternal_both", "maternal_just_repro", 
                                  "paternal_just_repro",
                                  "maternal_just_sterile", "paternal_just_sterile")

write.csv(overlapping_unique, file="overlapping_and_unique_matpat_repro_sterile.csv")




