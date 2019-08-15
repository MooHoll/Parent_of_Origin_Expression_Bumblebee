# Graph for SNPs called by freebayes 

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/1. Bombus epigenetics/mapping_stats_etc/snp_info")

library(ggplot2)
library(reshape2)

dat<-read.csv("snp_data.csv", header=T)
dat_subset<-dat[-c(2,4,5)]
dat_subset<-reshape2::melt(dat_subset)
head(dat_subset)

ggplot(data=dat_subset, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat="identity", position="fill")+
  ylab("SNP Count")+
  ggtitle("Number of SNPs per Sample")+
  scale_fill_manual(breaks=c("Total_SNPs","Homo_alt","Unique_to_male_or_queen"),
                    labels = c("Total SNPs","Homozygous Alternative","Unique"),
                    values= c("dodgerblue","skyblue2","turquoise"))

  
ggplot(data=dat, aes(x=Sample))+
  geom_bar(aes(y=Total_SNPs), stat="identity",fill="dodgerblue",show.legend=T)+
  geom_bar(aes(y=Homo_alt), stat="identity",fill="skyblue2",show.legend=T)+
  geom_bar(aes(y=Unique_to_male_or_queen), stat="identity",fill="turquoise",show.legend=T)+
  ylab("SNP Count")+
  geom_text(data = dat, aes(label = Total_SNPs, y = Total_SNPs -20000),
            size=5)+
  geom_text(data=dat, aes(label=Unique_to_male_or_queen, y=Unique_to_male_or_queen -20000),
            size=5)+
  geom_text(data = dat, aes(label = Homo_alt, y = Homo_alt - 20000),
            size=5)+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20))
  
## Made both graphs and put in powerpoint to put the legend of the first graph onto the second



males <- dat[1:4,]
females <- dat[5:8,]

sd(males$Total_SNPs) #144951
sd(females$Total_SNPs) #167392
