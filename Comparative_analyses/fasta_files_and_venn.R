###########################################
# Getting fasta seq for genes of interest (Eamonn's script)
# Also making venn diagram of overlapping genes
###########################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Imprinting_graphs")

#Get seqinr package. This lets you read the fasta files
#install.packages("seqinr") 
library(seqinr)
library(readr)
library(VennDiagram)

#Get the fasta file in and make it useable
genome<-read.fasta (file = "final_GFF_bter_1.0.fasta", as.string = TRUE, strip.desc=TRUE)
seq.data<-as.data.frame(do.call(rbind, genome))
seq.data<-data.frame(as.character(names(genome)),seq.data)
colnames(seq.data)=c("Genes","Sequences")

# -----------------------------------------
# Model 1 all genes 
# -----------------------------------------

# Subset the full gene list from the imprinting model to get the significant genes (850 total)
input <- read_csv("stats_imprinting_model1_total_counts_per_gene.csv")
input <- subset(input, input$signimprinbothcrosses==1)
query <- as.data.frame(unique(input$geneID))
colnames(query)<-"Genes"

signmatbias = (input$signimprinbothcrosses==1&input$avgpropmatexpr>0.6)
signpatbias = (input$signimprinbothcrosses==1&input$avgpropmatexpr<0.4)

# 90 genes
maternal_bias <- as.data.frame(unique(input$geneID[signmatbias]))
colnames(maternal_bias)<-"Genes"

# 64 genes
paternal_bias <- as.data.frame(unique(input$geneID[signpatbias]))
colnames(paternal_bias)<-"Genes"

# Merging for fasta info
output <- merge(query,seq.data,by="Genes")
output_mat <- merge(maternal_bias,seq.data,by="Genes")
output_pat <- merge(paternal_bias,seq.data,by="Genes")

#Getting it out of R
write.csv(output, file="imprinted_genes_model1_850total.csv")
write.csv(output_mat, file="imprinted_genes_model1_maternal_90total.csv")
write.csv(output_pat, file="imprinted_genes_model1_paternal_64total.csv")

write.fasta(as.list(output$Sequences), output$Genes, open = "w", nbchar = 100,file.out="imprinted_genes_model1_850total.fasta")
write.fasta(as.list(output_mat$Sequences), output_mat$Genes, open = "w", nbchar = 100,file.out="imprinted_genes_model1_maternal_90total.fasta")
write.fasta(as.list(output_pat$Sequences), output_pat$Genes, open = "w", nbchar = 100,file.out="imprinted_genes_model1_paternal_64total.fasta")


# -----------------------------------------
# Model 2 all genes 
# -----------------------------------------

# Subset the full gene list from the imprinting model to get the significant genes
input <- read_csv("stats_imprinting_model2_total counts per gene reproductive nonreproductive workers.csv")

signbias_repr = (input$signinbothdirs_repr==1)
signbias_nonrepr = (input$signinbothdirs_nonrepr==1)

# 170 genes 
query_nonrepro = as.data.frame(unique(input$geneID[signbias_nonrepr&
                                (input$avgpropmatexpr>0.6|input$avgpropmatexpr<0.4)]))
colnames(query_nonrepro)<-"Genes"

# 163 genes 
query_repro = as.data.frame(unique(input$geneID[signbias_repr&
                                   (input$avgpropmatexpr>0.6|input$avgpropmatexpr<0.4)]))
colnames(query_repro)<-"Genes"

# 699 genes 
query_repro_all <- as.data.frame(unique(input$geneID[signbias_repr]))
colnames(query_repro_all)<-"Genes"

# 746 genes 
query_nonrepro_all <- as.data.frame(unique(input$geneID[signbias_nonrepr]))
colnames(query_nonrepro_all)<-"Genes"



# Getting fasta seq for those genes
output_repro <- merge(query_repro,seq.data,by="Genes")
output_nonrepro <- merge(query_nonrepro,seq.data,by="Genes")
output_repro_all <- merge(query_repro_all,seq.data,by="Genes")
output_nonrepro_all <- merge(query_nonrepro_all,seq.data,by="Genes")

#Getting it out of R
write.csv(output_repro, file="imprinted_genes_model2_repro_163total.csv")
write.csv(output_nonrepro, file="imprinted_genes_model2_nonrepro_170total.csv")
write.csv(output_repro_all, file="imprinted_genes_model2_repro_699total.csv")
write.csv(output_nonrepro_all, file="imprinted_genes_model2_nonrepro_746total.csv")


write.fasta(as.list(output_repro$Sequences), output_repro$Genes, open = "w", nbchar = 100,file.out="imprinted_genes_model2_repro_163total.fasta")
write.fasta(as.list(output_nonrepro$Sequences), output_nonrepro$Genes, open = "w", nbchar = 100,file.out="imprinted_genes_model2_nonrepro_170total.fasta")
write.fasta(as.list(output_repro_all$Sequences), output_repro_all$Genes, open = "w", nbchar = 100,file.out="imprinted_genes_model2_repro_699total.fasta")
write.fasta(as.list(output_nonrepro_all$Sequences), output_nonrepro_all$Genes, open = "w", nbchar = 100,file.out="imprinted_genes_model2_nonrepro_746total.fasta")


# -----------------------------------------
# Making a Venn 
# -----------------------------------------

imprinted_mat <- output_mat$Genes
imprinted_pat <- output_pat$Genes
imprinted_repro_sig <- output_repro$Genes
imprinted_nonrepro_sig <- output_nonrepro$Genes

list1 <- list(imprinted_mat, imprinted_nonrepro_sig, imprinted_pat, imprinted_repro_sig)
overlap <- calculate.overlap(list1)

# Had to run this for all comparisons to manually put in numbers below
length(Reduce(intersect, list(imprinted_nonrepro_sig,imprinted_pat, imprinted_repro_sig)))

grid.newpage()
venn.plot.new <- draw.quad.venn(
  area1 = 90,
  area2 = 170,
  area3 = 64,
  area4 = 163,
  n12 = 71,
  n13 = 0,
  n14 = 76,
  n23 = 56,
  n24 = 149,
  n34 = 57,
  n123 = 0,
  n124 = 69,
  n134 = 0,
  n234 = 55,
  n1234 = 0,
  category = c("Maternal", "Sterile", "Paternal", "Reproductive"),
  fill = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 4,
  margin = 0.05,
  ind = TRUE,
  cex = 3
)

