########################################################
####### Making Nice Imprinting Graphs for Poster 
########################################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Final_imprinting_model/Imprinting_graphs")
library(readr)
library(pbapply)
library(plyr)
library(doBy)
library(devtools)
# install_github("dgrtwo/broom")
library(broom)
library(dplyr)
library(afex)
library(lsmeans)
library(car)
library(ggplot2)
library(ggthemes)
library(tidyr)

# -----------------------------------------
###### Model 2: Imprinting by Repro Status
# -----------------------------------------

# Read in data
model2 <- read_csv("stats_imprinting_model2_total counts per gene reproductive nonreproductive workers.csv")
model2<- model2[,-1]

# Make lists of significantly bias genes
model2<-as.data.frame(model2)

df = summaryBy(sign + matexprbias ~ geneID+reprstatus1, data=model2, FUN=sum)
signimprinbothcrosses_repr = df[df$reprstatus1=="REPR",]$geneID[(df[df$reprstatus1=="REPR",]$sign.sum==2&df[df$reprstatus1=="REPR",]$matexprbias.sum==2)|
                                                                  (df[df$reprstatus1=="REPR",]$sign.sum==2&df[df$reprstatus1=="REPR",]$matexprbias.sum==0)]
signimprinbothcrosses_nonrepr = df[df$reprstatus1=="NONREPR",]$geneID[(df[df$reprstatus1=="NONREPR",]$sign.sum==2&df[df$reprstatus1=="NONREPR",]$matexprbias.sum==2)|
                                                                        (df[df$reprstatus1=="NONREPR",]$sign.sum==2&df[df$reprstatus1=="NONREPR",]$matexprbias.sum==0)]

# convert from long to wide format to plot
model2$condition = interaction(model2$direction_cross,model2$reprstatus1)
model2$signinbothdirs = ((model2$signinbothdirs_repr==1)|(model2$signinbothdirs_nonrepr==1))*1
model2_wide = spread(model2[,c(1,2,3,4,10,17)], direction_cross, propmatexpr)
model2_wide$signinbothdirs[model2_wide$reprstatus1=="REPR"] = 0
model2_wide$signinbothdirs[model2_wide$reprstatus1=="REPR"&(model2_wide$geneID %in% signimprinbothcrosses_repr)] = 1
model2_wide$signinbothdirs[model2_wide$reprstatus1=="NONREPR"] = 0
model2_wide$signinbothdirs[model2_wide$reprstatus1=="NONREPR"&(model2_wide$geneID %in% signimprinbothcrosses_nonrepr)] = 1
model2_wide$reprstatus1=factor(model2_wide$reprstatus1,levels=c("NONREPR","REPR"),labels=c("Sterile","Reproductive"))
head(model2_wide)


qplot(x=initial,
      y=reciprocal,
      data=model2_wide,
      #colour=ifelse((signinbothdirs==1)&(avgpropmatexpr>0.6|avgpropmatexpr<0.4),I("red"),I("grey")),
      col=ifelse((signinbothdirs==1)&(avgpropmatexpr>0.9|avgpropmatexpr<0.1), "red", 
                 ifelse((signinbothdirs==1)&(avgpropmatexpr>0.6|avgpropmatexpr<0.4), "blue", "grey")),
      size=I(3),
      pch = I(16),
      xlab = expression(paste("Maternal expression: Initial Cross ")),
      ylab = expression(paste("Maternal expression: Reciprocal Cross")),
      geom="point") +
  facet_grid(~reprstatus1) +
  scale_color_identity() +
  coord_equal() +
  geom_hline(yintercept=0.5, colour=I("grey")) +
  geom_vline(xintercept=0.5, colour=I("grey")) +
  theme_bw() +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22),
        strip.text.x = element_text(size = 22))+
  annotate("text",x=0.75,y=0.75, label="paste(Maternal)",size=10,parse=T)+
  annotate("text",x=0.25,y=0.25, label="paste(Paternal)",size=10,parse=T)+
  annotate("text",x=0.25,y=0.75, label="paste(Audax)",size=10,parse=T)+
  annotate("text",x=0.8,y=0.25, label="paste(Dalmatinus)",size=10,parse=T)
 # geom_text(data=model2_wide, 
  #          aes(label=ifelse((signinbothdirs==1)&(avgpropmatexpr>0.6|avgpropmatexpr<0.4),
   #                          as.character(geneID),""), colour="black"), hjust=0.5, vjust=2, size=2) 

signbias_repr <- (model2$signinbothdirs_repr==1)
signbias_nonrepr <- (model2$signinbothdirs_nonrepr==1)

#17 extreme
ext_sterile <- as.data.frame(unique(model2$geneID[signbias_nonrepr&
                                                         (model2$avgpropmatexpr>0.9|model2$avgpropmatexpr<0.1)]))
colnames(ext_sterile) <- "gene"
#17 extreme
ext_repro <- as.data.frame(unique(model2$geneID[signbias_repr&
                                                        (model2$avgpropmatexpr>0.9|model2$avgpropmatexpr<0.1)]))
colnames(ext_repro) <- "gene"
both_ext <- merge(ext_repro, ext_sterile, by="gene", all=T) #17 so top are in common in both sterile and repro


maternal_both_all_data <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Final_imprinting_model/GO_analysis/maternal_both_all_data.csv")


gene_names <- subset(maternal_both_all_data, maternal_both_all_data$Symbol %in% both_ext$gene)
# Majority uncharacterised 

## Testing if more maternal or paternal in general
observed = c(103, 81)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

chisq.test(x = observed,
           p = expected)


## Testing if more maternal in repro/sterile
observed = c(89, 93)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

chisq.test(x = observed,
           p = expected)
