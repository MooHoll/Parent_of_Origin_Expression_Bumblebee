###-------------------------------------------------------------------
# Sanity check the variation within the pooled data
###-------------------------------------------------------------------

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Final_imprinting_model")

library(readr)
library(ggplot2)
library(reshape)
library(doBy)

###-------------------------------------------------------------------
all_info <- read_delim("all_information_for_proportion_sanity_checks.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(all_info)
all_info$propmat <- 1-all_info$proppaternal

all_info$reprstatus1 <- "NONREPR" 
all_info$reprstatus1[all_info$worker_type=="DOM"|all_info$worker_type=="SUB"] <-"REPR"

###-------------------------------------------------------------------
# Distribution checks

# All
ggplot(all_info, aes(x=propmat))+
  geom_histogram(colour="black",bins = 50)+
  geom_vline(xintercept =0.4, size=1, colour="red")+
  geom_vline(xintercept = 0.6, size=1, colour="red")+
  theme_bw()+
  xlab("Proportion of Maternal Expression")+
  ylab("Number of Genes")+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20))

# by tissue
ggplot(all_info, aes(x=propmat, fill=tissue))+
  geom_histogram(colour="black",bins = 50)+
  geom_vline(xintercept =0.4, size=1, colour="red")+
  geom_vline(xintercept = 0.6, size=1, colour="red")+
  theme_bw()+
  xlab("Proportion of Maternal Expression")+
  ylab("Number of Genes")+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20))

# by work type
ggplot(all_info, aes(x=propmat, fill=worker_type))+
  geom_histogram(colour="black",bins = 50)+
  geom_vline(xintercept =0.4, size=1, colour="red")+
  geom_vline(xintercept = 0.6, size=1, colour="red")+
  theme_bw()+
  xlab("Proportion of Maternal Expression")+
  ylab("Number of Genes")+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20))

# by repro status
ggplot(all_info, aes(x=propmat, fill=reprstatus1))+
  geom_histogram(colour="black",bins = 50)+
  geom_vline(xintercept =0.4, size=1, colour="red")+
  geom_vline(xintercept = 0.6, size=1, colour="red")+
  theme_bw()+
  xlab("Proportion of Maternal Expression")+
  ylab("Number of Genes")+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20))

###-------------------------------------------------------------------
# Check what the varience in exp is when we expect no variation, i.e. between 
# libraries of the same tissue/worker type, how much does it deviate from 0.5 when not showing ASE?

# Get list of genes which were sig parent of origin exp
overallmeanpreds_with_sig <- read_delim("overallmeanpreds_with_sig.txt", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
head(overallmeanpreds_with_sig)
sig_genes_repro <- as.data.frame(unique(overallmeanpreds_with_sig$geneID[overallmeanpreds_with_sig$signinbothdirs_repr==1]))
sig_genes_sterile <- as.data.frame(unique(overallmeanpreds_with_sig$geneID[overallmeanpreds_with_sig$signinbothdirs_nonrepr==1]))

colnames(sig_genes_repro) <- "sig_genes"
colnames(sig_genes_sterile) <- "sig_genes"

all_sig <- merge(sig_genes_repro, sig_genes_sterile, all=T)

all_info_nonsig <- all_info[!(all_info$geneID %in% all_sig$sig_genes),]
head(all_info_nonsig)

ggplot(all_info_nonsig, aes(x=propmat))+
  geom_histogram(colour="black",bins = 50)+
  geom_vline(xintercept =0.4, size=1, colour="red")+
  geom_vline(xintercept = 0.6, size=1, colour="red")+
  theme_bw()+
  xlab("Proportion of Maternal Expression")+
  ylab("Number of Genes")+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20))

# Check what is the deviation from 0.5 for those which were deemed not to show PoO exp
all_info$sig <- "no"
all_info$sig[(all_info$geneID %in% sig_genes_repro$sig_genes) &
               all_info$reprstatus1 == "REPR"] <- "yes"
all_info$sig[(all_info$geneID %in% sig_genes_sterile$sig_genes) &
               all_info$reprstatus1 == "NONREPR"] <- "yes"

ggplot(all_info, aes(y=propmat, x = seq(1, length(all_info$propmat)), colour=sig))+
  geom_point()

# what I need to do is plot the deviation from 0.5 between replicates of sig and non-sig in like a boxplot
# step 1: calculate deviation between replicates...
gene1 <- all_info[all_info$geneID=="Def",]
gene1 <- as.data.frame(gene1)
look <- summaryBy(propmat ~ geneID+reprstatus1+tissue, data=gene1, FUN=sd)

all_info <- as.data.frame(all_info)
summarised_all_info <- summaryBy(propmat ~ geneID+reprstatus1, data=all_info, FUN=sd)
summarised_all_info$sig <- "no"
summarised_all_info$sig[summarised_all_info$geneID %in% all_sig$sig_genes] <- "yes"

ggplot(summarised_all_info, aes(x = sig, y= propmat.sd))+
  geom_boxplot()

ggplot(summarised_all_info, aes(y= propmat.sd))+
  geom_boxplot()



mean_all_info <- summaryBy(propmat ~ geneID+reprstatus1, data=all_info, FUN=mean)
head(mean_all_info)
head(summarised_all_info)

all_data_new <- merge(mean_all_info, summarised_all_info, by=c("geneID","reprstatus1"))
head(all_data_new)

all_data_new$sig[(all_data_new$propmat.mean > 0.6 | all_data_new$propmat.mean < 0.4) & 
                   all_data_new$geneID %in% all_sig$sig_genes ] <- "true_sig"

ggplot(all_data_new, aes(y= propmat.mean, x=propmat.sd, colour=sig))+
  geom_point()+
  theme_bw()+
  xlab("Standard Deviation")+
  ylab("Mean Proportion of Maternal Expression")+
  geom_hline(yintercept = 0.4, colour ="red")+
  geom_hline(yintercept = 0.6, colour ="red")+
  scale_colour_manual(limits = c("true_sig","yes","no"),
                    labels = c("Imprinted","Significant","Non-Significant"),
                    values = c("red","orange","grey48"))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        legend.title = element_blank(),
        legend.text = element_text(size=18))
table(all_data_new$sig)

  
gene_2 <- all_data_new[all_data_new$geneID=="Def",]

# Make the graph as Eamonn suggest
head(all_data_new)

data1 <- all_data_new
data1$sig <- "all_genes"

data2 <- all_data_new
#data2$sig[data2$sig=="true_sig"] <- "yes"
data2 <- data2[!data2$sig == "true_sig",]

data2$sig[data2$sig == "yes" & data2$propmat.mean < 0.5] <- "sig_pat"
data2$sig[data2$sig == "yes" & data2$propmat.mean > 0.5] <- "sig_mat"

data3 <- all_data_new
data3 <- data3[data3$sig =="true_sig",]

data3$sig[data3$propmat.mean < 0.5] <- "true_sig_pat"
data3$sig[data3$propmat.mean > 0.5] <- "true_sig_mat"

all_for_plot <- rbind(data1,data2,data3)

ggplot(all_for_plot, aes(x= sig, y=propmat.mean, fill=sig))+
  geom_boxplot()+
  theme_bw()+
  xlab("Gene Category")+
  ylab("Mean Proportion of Materal Expression")+
  scale_x_discrete(limits = c("all_genes","no","sig_mat","sig_pat","true_sig_mat","true_sig_pat"),
                    labels = c("All Genes", "None\nSignificant","Significant\nMaternal\nExpression",
                               "Significant\nPaternal\nExpression","Imprinted\nMaternal","Imprinted\nPaternal"))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        legend.position = "none")+
  geom_hline(yintercept = 0.4, colour = "red")+
  geom_hline(yintercept = 0.6, colour = "red")

head(all_for_plot)
table(all_for_plot$sig)

