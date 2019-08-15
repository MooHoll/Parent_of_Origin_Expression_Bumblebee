###########################################
# Getting fasta seq for genes of interest (Eamonn's script)
# Also making venn diagram of overlapping genes
###########################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Honeybee")

#Get seqinr package. This lets you read the fasta files
#install.packages("seqinr") 
library(seqinr)
library(readr)


#Get the fasta file in and make it useable

genome<-read.fasta(file = "a.mellifera_genes_fasta.fa", as.string = TRUE, strip.desc=TRUE)
seq.data<-as.data.frame(do.call(rbind, genome))
seq.data<-data.frame(as.character(names(genome)),seq.data)
colnames(seq.data)=c("Genes","Sequences")

# -----------------------------------------
# Imprinted Gene Lists 
# -----------------------------------------

# Subset the full gene list from the imprinting model to get the significant genes (850 total)
all_gene_lists <- read_csv("honeybee_genes_galbraith.csv")

diff_exp <- as.data.frame(unique(all_gene_lists$Diff_exp_repro_sterile))
colnames(diff_exp)<-"Genes"

repro_imprinted <- as.data.frame(unique(all_gene_lists$Repro_imprinted))
colnames(repro_imprinted)<-"Genes"

sterile_imprinted <- as.data.frame(unique(all_gene_lists$Sterile_imprinted))
colnames(sterile_imprinted)<-"Genes"

repro_diff_exp_overlap <- as.data.frame(unique(all_gene_lists$Diff_exp_repro_imprint_overlap))
colnames(repro_diff_exp_overlap)<-"Genes"

sterile_diff_exp_overlap <- as.data.frame(unique(all_gene_lists$Diff_exp_sterile_imprint_overlap))
colnames(sterile_diff_exp_overlap)<-"Genes"


# Merging for fasta info
output_diff_exp <- merge(diff_exp,seq.data,by="Genes")
output_repro_imprinted <- merge(repro_imprinted,seq.data,by="Genes")
output_sterile_imprinted <- merge(sterile_imprinted,seq.data,by="Genes")
output_repro_overlap_diffexp <- merge(repro_diff_exp_overlap,seq.data,by="Genes")
output_sterile_overlap_diffexp <- merge(sterile_diff_exp_overlap,seq.data,by="Genes")


#Getting it out of R
write.csv(output_diff_exp, file="diffexp_genes_fasta.csv")
write.csv(output_repro_imprinted, file="repro_imrptined_genes_fasta.csv")
write.csv(output_sterile_imprinted, file="sterile_imrptined_genes_fasta.csv")
write.csv(output_repro_overlap_diffexp, file="repro_imrptined_genes_overlap_diffexp_fasta.csv")
write.csv(output_sterile_overlap_diffexp, file="sterile_imrptined_genes_overlap_diffexp_fasta.csv")


write.fasta(as.list(output_diff_exp$Sequences), output_diff_exp$Genes, open = "w", nbchar = 100,file.out="diffexp_genes.fasta")
write.fasta(as.list(output_repro_imprinted$Sequences), output_repro_imprinted$Genes, open = "w", nbchar = 100,file.out="repro_imrptined_genes.fasta")
write.fasta(as.list(output_sterile_imprinted$Sequences), output_sterile_imprinted$Genes, open = "w", nbchar = 100,file.out="sterile_imrptined_genes.fasta")
write.fasta(as.list(output_repro_overlap_diffexp$Sequences), output_repro_overlap_diffexp$Genes, open = "w", nbchar = 100,file.out="repro_imrptined_genes_overlap_diffexp.fasta")
write.fasta(as.list(output_sterile_overlap_diffexp$Sequences), output_sterile_overlap_diffexp$Genes, open = "w", nbchar = 100,file.out="sterile_imrptined_genes_overlap_diffexp.fasta")
