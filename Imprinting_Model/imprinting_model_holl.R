###-------------------------------------------------------------------
# BINOMIAL GLMs TO DETERMINE IMPRINTING STATUS OF EACH GENE
###-------------------------------------------------------------------

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Final_imprinting_model")

library(pbapply)
library(plyr)
library(doBy)
library(export) 
library(effects)
library(devtools)
library(broom)
library(dplyr)
library(afex)ccc
library(lsmeans)
library(car)
library(ggplot2)
library(ggthemes)
library(tidyr)

###------------------------------------------------------------------
# Read in count data
###------------------------------------------------------------------

df = readRDS("maternal_filtered_count_table.rds") 
# Convert family id's into correct format
df$familyID=factor(df$familyID,levels=levels(df$familyID),labels=c("01_02","03_04"))
# Add a total count row
df$total = df$maternal+df$paternal

# Remove snps with less total coverage than 5
df = df[df$total>5,] 

# add two additional factors to represent reproductive / nonreproductive workers
df$reprstatus1="NONREPR" 
df$reprstatus1[df$worker_type=="DOM"|df$worker_type=="SUB"]="REPR"
df$reprstatus2="NONREPR"
df$reprstatus2[df$worker_type=="DOM"]="DOMREPR"
df$reprstatus2[df$worker_type=="SUB"]="SUBREPR"
df$reprstatus1=factor(df$reprstatus1,levels=c("NONREPR","REPR")) 
df$reprstatus2=factor(df$reprstatus2,levels=c("NONREPR","SUBREPR","DOMREPR")) 
df$worker_type=factor(df$worker_type,levels=c("NUR","FOR","SUB","DOM"))


###------------------------------------------------------------------
# Look at data
###------------------------------------------------------------------

# Look at the summed counts over all the SNP positions for each gene
 dftotal=summaryBy(maternal+paternal~geneID, data=df, FUN=sum)
 dftotal$total<-dftotal$maternal.sum +dftotal$paternal.sum
 dftotal$proppaternal=dftotal$paternal.sum/(dftotal$paternal.sum+dftotal$maternal.sum)
 range(dftotal$proppaternal) # 0.0000000-0.9688553
 hist(dftotal$proppaternal[dftotal$total>100],col="steelblue",
      main="Proportion of paternal expression",
      xlab="Proportion of paternal expression")

# Look at genes with extreme maternal/paternal bias
unique(dftotal$geneID[dftotal$proppaternal>0.9]) # paternal expr bias 
length(unique(dftotal$geneID[dftotal$proppaternal>0.9])) # 15 genes
unique(dftotal$geneID[dftotal$proppaternal<0.1]) # maternal expr bias 
length(unique(dftotal$geneID[dftotal$proppaternal<0.1])) # 78 genes

# Make a column taking into account all experimental information 
df$treatm = factor(paste(df$direction_cross, df$familyID, df$worker_type, df$tissue, sep="."))

# Look how many genes are in every condition within our experiment 
dftreatmspergene = summaryBy(maternal+paternal~geneID+treatm, data=df, FUN=length )
dftreatmspergene2 = summaryBy(maternal.length+paternal.length~geneID, data=dftreatmspergene, FUN=length )

hist(dftreatmspergene2$maternal.length.length, main="Number of individuals with data per gene", xlab="Nr of individuals")
genesinallcrossesandtissues = dftreatmspergene2$geneID[(dftreatmspergene2$maternal.length.length==32)]
length(genesinallcrossesandtissues) # 7412 genes from all conditions/crosses left - nice!

# Look at the number of SNPs per gene
dfnrSNPspergene = summaryBy(maternal+paternal~geneID+snp, data=df, FUN=length)
dfnrSNPspergene2 = summaryBy(maternal.length+paternal.length~geneID, data=dfnrSNPspergene, FUN=length)
head(dfnrSNPspergene2)
hist(dfnrSNPspergene2$maternal.length.length, main="Histogram of nr of SNPs per gene", xlab="Count")
geneswithatleast2SNPs = dfnrSNPspergene2$geneID[(dfnrSNPspergene2$maternal.length.length>=2)]
length(geneswithatleast2SNPs) # 10211 genes with at least two SNPs per gene

# Count total number of SNPs to be used in analysis
totalsnps <- dfnrSNPspergene2[dfnrSNPspergene2$geneID %in% genesinallcrossesandtissues,]
head(totalsnps)
sum(totalsnps$maternal.length.length) #1,143,218

###------------------------------------------------------------------
# Filter data to keep only genes in all crosses/conditions
###------------------------------------------------------------------


# Make a total maternal and paternal counts per gene across all tissues, 
# to be able to test overall imprinting per gene (split by reproductive status)
dftotcountspergene2 = summaryBy(maternal+paternal ~ geneID+direction_cross+familyID+reprstatus1, data=df, FUN=sum)

# See which genes occur in all experimental conditions/crosses
dfcrossfamreprcount = summaryBy(maternal.sum+paternal.sum ~ geneID, data=dftotcountspergene2, FUN=length)

# Keep only genes for which we have counts in all families & crosses
selectedgenes = dfcrossfamreprcount$geneID[dfcrossfamreprcount$maternal.sum.length==8]
length(selectedgenes) # this would retain 7508 genes
dftotcountspergene2 = dftotcountspergene2[dftotcountspergene2$geneID %in% selectedgenes,]


# Check the proportion spead and look at the most meternal/paternal proportions
dftotcountspergene2$total.sum = dftotcountspergene2$maternal.sum + dftotcountspergene2$paternal.sum
dftotcountspergene2$proppaternal = dftotcountspergene2$paternal.sum/dftotcountspergene2$total.sum
dftotcountspergene2$geneID = droplevels(dftotcountspergene2$geneID)
head(dftotcountspergene2)
hist(dftotcountspergene2$proppaternal,col="steelblue", 
     main="Proportion of paternal expression",
     xlab="Proportion of paternal expression", breaks=40)
range(dftotcountspergene2$proppaternal) # 0.0000000-0.9956617
unique(dftotcountspergene2$geneID[dftotcountspergene2$proppaternal>0.9]) # paternal expr bias: 35
length(unique(dftotcountspergene2$geneID[dftotcountspergene2$proppaternal>0.9]))
unique(dftotcountspergene2$geneID[dftotcountspergene2$proppaternal<0.1]) # maternal expr bias: 70 
length(unique(dftotcountspergene2$geneID[dftotcountspergene2$proppaternal<0.1]))


###------------------------------------------------------------------
# Run model
###------------------------------------------------------------------


# MODEL 2: ANALYSIS OF OVERALL IMPRINTING PER GENE ACROSS ALL TISSUES BASED ON DATAFRAME dftotcountspergene2
# (split up by workers' reproductive status)
# (Tom also tried model 1. not split by repro status)
set_sum_contrasts()
head(dftotcountspergene2)
length(unique(dftotcountspergene2$geneID)) # 7508

# Test with 1 gene
gene = "LOC100650436" # the real vitellogenin? (more highly expressed in repr workers)
dfsubs = dftotcountspergene2[dftotcountspergene2$geneID==gene,]
dfsubs$maternal.sum[dfsubs$maternal.sum==0] = 1 # to avoid problems of complete separation (not really relevent for maternal as we removed all 0 anyway earlier)
dfsubs$paternal.sum[dfsubs$paternal.sum==0] = 1 # to avoid problems of complete separation

# Take into account overdispersion using quasibinomial
fit = glm( cbind(maternal.sum, paternal.sum) ~ direction_cross+familyID+reprstatus1, 
           family = quasibinomial(link=logit), data=dfsubs)
summary(fit) # Wald z tests
Anova(fit,type="III") # type III Anova tests

# Proportion of maternal expression in different crosses & families according to repr status
meanpreds = data.frame(summary(lsmeans(fit,~ direction_cross+familyID+reprstatus1)))
meanpreds$z = meanpreds$lsmean/meanpreds$SE # z scores for deviation from 50:50
meanpreds$p = 2*pnorm(abs(meanpreds$z),lower.tail = FALSE) # 2 sided p values for deviation from 50:50 (not corrected for mult testing)
meanpreds$lsmean = plogis(meanpreds$lsmean) # backtransform to original scale
meanpreds$asymp.LCL = plogis(meanpreds$asymp.LCL)
meanpreds$asymp.UCL = plogis(meanpreds$asymp.UCL)
meanpreds
# or more elegantly:
meanpreds=data.frame(lsmeans::test(lsmeans(fit,~direction_cross*familyID*reprstatus1,type="response"), adjust = "none"))[,-c(5,6)] # without multiplicity correction, otherwise use adjust = "mvt"
cints=data.frame(confint(lsmeans(fit,~direction_cross*familyID*reprstatus1,type="response")))[,c(7:8)]
meanpreds=data.frame(meanpreds[,c(1:4)],cints,meanpreds[,c(5:6)])
meanpreds

# Overall mean proportion of maternal expression in 2 reciprocal crosses across both families for repr & nonrepr workers (pooling two familiies together)
overallmeanpreds = data.frame(summary(lsmeans(fit,~ direction_cross*reprstatus1)))
overallmeanpreds$z = overallmeanpreds$lsmean/overallmeanpreds$SE # z scores for deviation from 50:50
overallmeanpreds$p = 2*pnorm(abs(overallmeanpreds$z),lower.tail = FALSE) # 2 sided p values for deviation from 50:50 (not corrected for mult testing)
overallmeanpreds$lsmean = plogis(overallmeanpreds$lsmean) # backtransform to original scale
overallmeanpreds$asymp.LCL = plogis(overallmeanpreds$asymp.LCL)
overallmeanpreds$asymp.UCL = plogis(overallmeanpreds$asymp.UCL)
overallmeanpreds
# or more elegantly:
overallmeanpreds=data.frame(lsmeans::test(lsmeans(fit,~direction_cross*reprstatus1,type="response"), adjust = "none"))[,-c(4,5)] # without multiplicity correction, otherwise use adjust = "mvt"
cints=data.frame(confint(lsmeans(fit,~direction_cross*reprstatus1,type="response")))[,c(6:7)]
overallmeanpreds=data.frame(overallmeanpreds[,c(1:3)],cints,overallmeanpreds[,c(4:5)])
overallmeanpreds

# Look how proportion of expression changes given direction cross, family, repro status
plot(allEffects(fit), type="response", ylab="Proportion of paternal expression")

# Model for all genes
overallmeanpreds = do.call(rbind,
                           pblapply(levels(dftotcountspergene2$geneID),
                                    function (gene) { 
                                      dfsubs = dftotcountspergene2[dftotcountspergene2$geneID==gene,]
                                      dfsubs$maternal.sum[dfsubs$maternal.sum==0] = 1 # to avoid problems of complete separation
                                      dfsubs$paternal.sum[dfsubs$paternal.sum==0] = 1 # to avoid problems of complete separation
                                      fit = glm( cbind(maternal.sum, paternal.sum) ~ 
                                                   direction_cross+familyID+reprstatus1, 
                                                 family = quasibinomial(link=logit), data=dfsubs)
                                      # overall mean proportion of maternal expression in 2 reciprocal crosses across both families
                                      overallmeanpreds = data.frame(summary(suppressMessages(lsmeans(fit,~ direction_cross*reprstatus1))))
                                      overallmeanpreds$z = overallmeanpreds$lsmean/overallmeanpreds$SE # z scores for deviation from 50:50
                                      overallmeanpreds$p = 2*pnorm(abs(overallmeanpreds$z),lower.tail = FALSE) # 2 sided p values for deviation from 50:50 (not corrected for mult testing)
                                      overallmeanpreds$lsmean = plogis(overallmeanpreds$lsmean) # backtransform to original scale
                                      overallmeanpreds$asymp.LCL = plogis(overallmeanpreds$asymp.LCL)
                                      overallmeanpreds$asymp.UCL = plogis(overallmeanpreds$asymp.UCL)
                                      df=data.frame(geneID=gene,overallmeanpreds)
                                      df=df[,-6]
                                      colnames(df)=c("geneID","direction_cross","reprstatus1","propmatexpr","SE","propmatexpr.LCL","propmatexpr.UCL","z","p")
                                      df$avgpropmatexpr=mean(df$propmatexpr)
                                      df
                                    }
                           ) )

# Benjamini-Hochberg correction for multiple testing
overallmeanpreds$padj = p.adjust(overallmeanpreds$p, method="BH") 

# Add a column showing the significant rows
overallmeanpreds$sign = (overallmeanpreds$padj<0.05)*1

# Add a column showing if the expression bias is materal (1) or paternal (0)
overallmeanpreds$matexprbias = (overallmeanpreds$propmatexpr>0.5)*1

# Pool data across the direction of the cross
df1 = summaryBy(sign + matexprbias ~ geneID+reprstatus1, data=overallmeanpreds, FUN=sum)

# Pull out data that is signifiacnt in both crosses in the same direction (both maternal/paternal bias) in repro/sterile
signimprinbothcrosses_repr = df1[df1$reprstatus1=="REPR",]$geneID[(df1[df1$reprstatus1=="REPR",]$sign.sum==2&df1[df1$reprstatus1=="REPR",]$matexprbias.sum==2)|
                                                                  (df1[df1$reprstatus1=="REPR",]$sign.sum==2&df1[df1$reprstatus1=="REPR",]$matexprbias.sum==0)]
signimprinbothcrosses_nonrepr = df1[df1$reprstatus1=="NONREPR",]$geneID[(df1[df1$reprstatus1=="NONREPR",]$sign.sum==2&df1[df1$reprstatus1=="NONREPR",]$matexprbias.sum==2)|
                                                                        (df1[df1$reprstatus1=="NONREPR",]$sign.sum==2&df1[df1$reprstatus1=="NONREPR",]$matexprbias.sum==0)]
length(signimprinbothcrosses_repr ) # 700
length(signimprinbothcrosses_nonrepr ) # 747

# Have a look at the overlap
library(venn)
venn(list(signimprinbothcrosses_repr,signimprinbothcrosses_nonrepr), snames=c("reproductive W","nonreproductive W"), 
     ilab=TRUE, zcolor = "style", cexil=2, cexsn=2)


# Add new columns to show if the row is significant in both directions for repro/sterile
overallmeanpreds$signinbothdirs_repr = 0
overallmeanpreds$signinbothdirs_repr[overallmeanpreds$geneID %in% signimprinbothcrosses_repr] = 1
overallmeanpreds$signinbothdirs_nonrepr = 0
overallmeanpreds$signinbothdirs_nonrepr[overallmeanpreds$geneID %in% signimprinbothcrosses_nonrepr] = 1


# Pull out the repro/sterile rows which are signigicant in both direction
signbias_repr = (overallmeanpreds$signinbothdirs_repr==1)
signbias_nonrepr = (overallmeanpreds$signinbothdirs_nonrepr==1)

# Filter these rows with more stringent proportions of expression
signbiasedgenes_repr = unique(overallmeanpreds$geneID[signbias_repr&(overallmeanpreds$avgpropmatexpr>0.6|overallmeanpreds$avgpropmatexpr<0.4)])
length(signbiasedgenes_repr) # 163
signbiasedgenes_nonrepr = unique(overallmeanpreds$geneID[signbias_nonrepr&(overallmeanpreds$avgpropmatexpr>0.6|overallmeanpreds$avgpropmatexpr<0.4)])
length(signbiasedgenes_nonrepr) # 170

#write.csv(overallmeanpreds,"imprinting_model2_total counts per gene reproductive nonreproductive workers.csv")

# Filter these rows with more stringent proportions of expression (have a look at top ones)
top_repr = as.data.frame(unique(overallmeanpreds$geneID[signbias_repr&(overallmeanpreds$avgpropmatexpr>0.9|overallmeanpreds$avgpropmatexpr<0.1)]))
length(top_repr) # 17
colnames(top_repr) <- "gene"
#[1] LOC100642787 LOC100649867 LOC100651112 LOC100651450 LOC100651509
#[6] LOC100652265 LOC105666923 LOC105666969 LOC105666994 LOC105666995
#[11] LOC105667016 LOC105667042 LOC105667069 LOC105667070 LOC105667072
#[16] LOC105667073 LOC105667082
top_nonrepr = unique(overallmeanpreds$geneID[signbias_nonrepr&(overallmeanpreds$avgpropmatexpr>0.9|overallmeanpreds$avgpropmatexpr<0.1)])
length(top_nonrepr) # 17 (same as above)

look <- overallmeanpreds[overallmeanpreds$geneID %in% top_repr$gene, ]# all maternally exp bias

###------------------------------------------------------------------
# Make plot
###------------------------------------------------------------------

# Convert from long to wide format dataframe to plot
overallmeanpreds$condition = interaction(overallmeanpreds$direction_cross,overallmeanpreds$reprstatus1)
overallmeanpreds$signinbothdirs = ((overallmeanpreds$signinbothdirs_repr==1)|(overallmeanpreds$signinbothdirs_nonrepr==1))*1
overallmeanpreds_wide = spread(overallmeanpreds[,c(1,2,3,4,10,17)], direction_cross, propmatexpr)

overallmeanpreds_wide$signinbothdirs[overallmeanpreds_wide$reprstatus1=="REPR"] = 0
overallmeanpreds_wide$signinbothdirs[overallmeanpreds_wide$reprstatus1=="REPR"&(overallmeanpreds_wide$geneID %in% signimprinbothcrosses_repr)] = 1
overallmeanpreds_wide$signinbothdirs[overallmeanpreds_wide$reprstatus1=="NONREPR"] = 0
overallmeanpreds_wide$signinbothdirs[overallmeanpreds_wide$reprstatus1=="NONREPR"&(overallmeanpreds_wide$geneID %in% signimprinbothcrosses_nonrepr)] = 1
overallmeanpreds_wide$reprstatus1=factor(overallmeanpreds_wide$reprstatus1,levels=c("NONREPR","REPR"),labels=c("Nonreproductive workers","Reproductive workers"))
head(overallmeanpreds_wide)

# Make a nice plot
qplot(x=initial,
      y=reciprocal,
      data=overallmeanpreds_wide,
      colour=ifelse((signinbothdirs==1)&(avgpropmatexpr>0.6|avgpropmatexpr<0.4),I("red"),I("black")),
      size=I(3),
      pch = I(16),
      xlab = expression(paste("Prop. maternal expression  ", italic("B. t. dalmatinus"), " queen x ", italic("B. t. audax")," male")),
      ylab = expression(paste("Prop. maternal expression  ", italic("B. t. audax"), " queen x ", italic("B. t. dalmatinus")," male")),
      geom="point") +
  facet_grid(~reprstatus1) +
  scale_color_identity() +
  coord_equal() +
  geom_hline(yintercept=0.5, colour=I("grey")) +
  geom_vline(xintercept=0.5, colour=I("grey")) +
  geom_text(data=overallmeanpreds_wide, 
            aes(label=ifelse((signinbothdirs==1)&(avgpropmatexpr>0.6|avgpropmatexpr<0.4),
                             as.character(geneID),""), colour="black"), hjust=0.5, vjust=2, size=2) 

