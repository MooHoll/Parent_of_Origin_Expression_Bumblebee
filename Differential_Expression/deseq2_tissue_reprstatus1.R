##### Differential Expression Analysis, Repro and NonRepro Workers: DESeq2
# Edited from Anneleenes script Leuven, also see https://www.bioconductor.org/help/workflows/rnaseqGene/#aligning-reads-to-a-reference-genome

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/1. Bombus epigenetics/differential_expression_reproductive_status")

library(DESeq2)
library(reshape2)
library(pheatmap)
library(genefilter)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(PoiClaClu)
library(ReporteRsjars)
library(ReporteRs)#needed for export
library(rtable)#needed for export
library(rgl)#needed for export
library(export)
library(ggbeeswarm)
library(readr)

#################################################################################
# RUNNING FOR JUST ABD or KOP SAMPLES for REPRO/NONREPRO

# Count table made in .R script "count_table_diff_exp.R"
countdata <- read.csv("~/Dropbox/PhD Folder/Virtual Lab Book/1. Bombus epigenetics/differential_expression_reproductive_status/counts_diff_exp.csv")
# exclude these genes as they are massivly overexpressed in one/two samples compared to evrything else
countdata<-countdata[!(countdata$geneID== "LOC100645800"),]
#countdata<-countdata[!(countdata$geneID== "LOC100644839"),]
head(countdata)

# subset out count table to get metadata (rows of metadata must == columns of count data, so number of samples must match)
coldata<-countdata[,c(2,3,4,5,6,7,8,9,10,13,14)]
head(coldata)
coldata1<-coldata[!duplicated(coldata), ]
nrow(coldata1)

# Make count table matrix
countmatr<-countdata[,c(2,11,12)]
countmatr<-reshape2::melt(countmatr)
countmatr1<-dcast(countmatr, geneID ~ sample)
head(countmatr1)
countmatr1<-countmatr1[-c(1:5),]# Remove ambiguous genes etc
row.names(countmatr1)<-(countmatr1$geneID)
countmatr1<-countmatr1[,-1]
countmatr<-as.matrix(countmatr1)

# Change nest_worker to colony
coldata1$colony <- sub("_.*","",coldata1$nest_worker)
head(coldata1)
# Add new column for tissue with reprostatus, can then check interactions later separately 
coldata1$tissue_reprstatus1 = interaction(coldata1$tissue,coldata1$reprstatus1)
head(coldata1)
coldata1$weight<-as.numeric(coldata1$weight)
coldata1$age<-as.numeric(coldata1$age)
coldata1$colony<-as.factor(coldata1$colony)
# make a DESeq object taking into account colony and spliting by tissue/reprostatus (couldn't use colony as it's a combination of direction_cross and familyID)
dds= DESeqDataSetFromMatrix(countData = countmatr, colData = coldata1, design = ~ familyID+age+weight+direction_cross+tissue_reprstatus1)

#################################################################################

# remove features with low counts
dds = dds[ rowMeans(counts(dds)) > 10, ] 
nrow(dds) #9790 genes left

# rlog transform counts
rld = rlog(dds, blind=FALSE)

# PCA plot
data = plotPCA(rld, intgroup = c("tissue_reprstatus1"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=tissue_reprstatus1)) + geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+
  scale_color_manual(breaks = c("ABD.NONREPR","HEAD.NONREPR","ABD.REPR","HEAD.REPR"),
                     values=c("#00CC66","#66FFCC","#0066CC","#66CCFF"),
                     labels=c("Abdomen Sterile","Head Sterile", 
                              "Abdomen Repro", "Head Repro"))+
 #                    name=c("Tissue and Reproductive Status"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20),legend.title = element_blank())
#element_text(size=18))
#graph2ppt(file="pca_repro_nonrepro.pptx")


# two 1st samples plotted against each other to check consistency (for rlog and log2) for kop
par( mfrow = c( 1, 2 ) )
dds = estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,c(2,4)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
plot(assay(rld)[,c(2,4)],
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
par( mfrow = c( 1, 1 ) )
#graph2ppt(file="check_log_tranform_2_samples.pptx")

#################################################################################

# estimate size factors = normalize for library size
dds = DESeq2::estimateSizeFactors(dds)
dds = estimateDispersions(dds)
plotDispEsts(dds, xlab= "Mean of Normalised Counts",
             ylab= "Dispersion", cex=1.0, cex.lab=1.45, cex.axis=1.45)
#graph2ppt(file="check_DESeq2_dispersion.pptx")

# check sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )

new_colnames<-c("02_10_ABD_REPRO" , "02_10_HEAD_REPRO" ,"02_12_ABD_REPRO"  ,"02_12_HEAD_REPRO" ,"02_18_ABD_REPRO" ,
                "02_18_HEAD_REPRO", "02_22_ABD_REPRO"  ,"02_22_HEAD_REPRO" ,"02_35_ABD_NONREPRO"  ,"02_35_HEAD_NONREPRO",
                "02_37_ABD_NONREPRO" , "02_37_HEAD_NONREPRO" ,"02_44_ABD_NONREPRO"  ,"02_44_HEAD_NONREPRO" ,"02_54_ABD_NONREPRO", 
                "02_54_HEAD_NONREPRO", "12_74_ABD_NONREPRO"  ,"12_74_HEAD_NONREPRO" ,"12_75_ABD_NONREPRO"  ,"12_75_HEAD_NONREPRO",
                "12_76_ABD_NONREPRO" , "12_76_HEAD_NONREPRO" ,"12_78_ABD_NONREPRO"  ,"12_78_HEAD_NONREPRO" ,"12_82_ABD_REPRO", 
                "12_82_HEAD_REPRO", "12_84_ABD_REPRO"  ,"12_84_HEAD_REPRO" ,"12_85_ABD_REPRO"  ,"12_85_HEAD_REPRO",
                "12_89_ABD_REPRO" , "12_89_HEAD_REPRO" ,"22_09_ABD_REPRO"  ,"22_09_HEAD_REPRO" ,"22_18_ABD_REPRO", 
                "22_18_HEAD_REPRO", "22_21_ABD_REPRO"  ,"22_21_HEAD_REPRO" ,"22_23_ABD_REPRO"  ,"22_23_HEAD_REPRO",
                "22_26_ABD_NONREPRO" , "22_26_HEAD_NONREPRO" ,"22_30_ABD_NONREPRO"  ,"22_30_HEAD_NONREPRO" ,"22_35_ABD_NONREPRO", 
                "22_35_HEAD_NONREPRO", "22_49_ABD_NONREPRO"  ,"22_49_HEAD_NONREPRO" ,"31_44_ABD_NONREPRO"  ,"31_44_HEAD_NONREPRO",
                "31_45_ABD_NONREPRO" , "31_45_HEAD_NONREPRO" ,"31_64_ABD_NONREPRO"  ,"31_64_HEAD_NONREPRO" ,"31_65_ABD_NONREPRO", 
                "31_65_HEAD_NONREPRO", "31_79_ABD_REPRO"  ,"31_79_HEAD_REPRO" ,"31_84_ABD_REPRO"  ,"31_84_HEAD_REPRO",
                "31_86_ABD_REPRO" , "31_86_HEAD_REPRO" ,"31_89_ABD_REPRO"  ,"31_89_HEAD_REPRO")
new_colnames<-gsub("NONREPRO","Sterile", new_colnames)
new_colnames<-gsub("REPRO","Repro", new_colnames)
new_colnames<-gsub("HEAD","Head", new_colnames)
new_colnames<-gsub("ABD","Abd", new_colnames)
colnames(sampleDistMatrix) <- new_colnames
colnames(sampleDistMatrix)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_colnames = T,
         fontsize = 14)
#graph2ppt(file="sample_distances.pptx")

# check sample distances using the poisson method
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- new_colnames
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         show_colnames = T,
         fontsize = 14)
#graph2ppt(file="sample_distances_poisson.pptx")

#################################################################################

# differential expression (used Benjamini-Hochberg adjustment and a p-val of <0.1)
dds = DESeq2::DESeq(dds, parallel=TRUE)
resultsNames(dds)

# combinations of reprostatus and tissue 
res=results(dds, contrast=c("tissue_reprstatus1","ABD.REPR","ABD.NONREPR"))
summary(res) # 3854 (39%) upreg and 4501 (46%) downreg
#res_less_0.01<-subset(res, res$padj<0.01) 

res_plot<-lfcShrink(dds, contrast=c("tissue_reprstatus1","ABD.REPR","ABD.NONREPR"), res=res) #shrink log2 fold change to make comparison easier
plotMA(res_plot, ylim=c(-5,5), main="Repro vs Sterile (Abdomen Tissue)", cex=1.0, cex.lab=1.45, cex.axis=1.45, xlab='Mean of Normalised Counts',
       ylab='Log Fold Change') #distribution of coefficents of the model
#graph2ppt(file="MA_plot_ABD.pptx")
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
     cex.axis=1.45, cex.lab=1.45, cex.main=1.45, main="Repro vs Sterile (Abdomen Tissue)") # plot of p-vals excluding genes with very small counts
#graph2ppt(file="ABD_p-val_distribution_exclusing_small_counts.pptx")

resOrd=res[order(res$log2FoldChange),]
nrow(resOrd)
write.csv(as.data.frame(resOrd), file="ABD_all_genes_results.csv")

resOrd_significant<-subset(resOrd, padj<0.05)
nrow(resOrd_significant)#8097
write.csv(as.data.frame(resOrd_significant), file="ABD_diff_exp_genes.csv")

upreg_repro_sig<-subset(resOrd_significant, log2FoldChange>0 )
nrow(upreg_repro_sig)#3732
write.csv(upreg_repro_sig, file="ABD_upreg_in_repro.csv")

upreg_nonrepro_sig<-subset(resOrd_significant, log2FoldChange<0 )
nrow(upreg_nonrepro_sig)#4365
write.csv(upreg_nonrepro_sig, file="ABD_upreg_in_nonrepro.csv")

###

res1=results(dds, contrast=c("tissue_reprstatus1","HEAD.REPR","HEAD.NONREPR"))
summary(res1) # huge difference 235 (2.4% ) up and 150 (1.6%) down
#res1_less_0.01<-subset(res1, res1$padj<0.01)
res_plot1<-lfcShrink(dds, contrast=c("tissue_reprstatus1", "HEAD.REPR", "HEAD.NONREPR"), res1=res1) #shrink log2 fold change to make comparison easier
plotMA(res_plot1, ylim=c(-5,5), main="Repro vs Sterile (Head Tissue)", cex=1.0, cex.lab=1.45, cex.axis=1.45, xlab='Mean of Normalised Counts',
       ylab='Log Fold Change') #distribution of coefficents of the model
#graph2ppt(file="MA_plot_KOP.pptx")
hist(res1$pvalue[res1$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
     cex.axis=1.45, cex.lab=1.45, cex.main=1.45, main="Repro vs Sterile (Head Tissue)") # plot of p-vals excluding genes with very small counts
#graph2ppt(file="KOP_p-val_distribution_exclusing_small_counts.pptx")

resOrd1=res1[order(res1$log2FoldChange),]
nrow(resOrd1)
write.csv(as.data.frame(resOrd1), file="HEAD_all_genes_results.csv")

resOrd_significant<-subset(resOrd1, padj<0.05)
nrow(resOrd_significant)#279
write.csv(as.data.frame(resOrd_significant), file="HEAD_diff_exp_genes.csv")

upreg_repro_sig<-subset(resOrd_significant, log2FoldChange>0 )
nrow(upreg_repro_sig)#174
write.csv(upreg_repro_sig, file="HEAD_upreg_in_repro.csv")

upreg_nonrepro_sig<-subset(resOrd_significant, log2FoldChange<0 )
nrow(upreg_nonrepro_sig)#105
write.csv(upreg_nonrepro_sig, file="HEAD_upreg_in_nonrepro.csv")


###

# list of genes that are differential between repr and nonrepr W in either head or abdomen 
n=50
topdiff_down_abdreprw = head(c(1:nrow(res))[order(res$log2FoldChange)],n) # rownames(res)
topdiff_up_abdreprw = tail(c(1:nrow(res))[order(res$log2FoldChange)],n)
topdiff_abdreprw = union(topdiff_down_abdreprw,topdiff_up_abdreprw )
topdiff_down_headreprw = head(c(1:nrow(res1))[order(res1$log2FoldChange)],n) # rownames(res)
topdiff_up_headreprw = tail(c(1:nrow(res1))[order(res1$log2FoldChange)],n)
topdiff_headreprw = union(topdiff_down_headreprw,topdiff_up_headreprw )

################################################################################# 

# topVarGenes = head(order(rowVars(assay(rld)),decreasing=TRUE),50)




# heatmap of top differential genes in abdomen of repr vs nonrepr workers
my_colors = list(
  Status = c(Sterile = "#00CC66", Reproductive ="#0066CC"))

mat = assay(rld)[ topdiff_abdreprw, grepl("ABD",colnames(assay(rld)))]
mat = mat - rowMeans(mat)
df = as.data.frame(colData(rld)[,c("reprstatus1"),drop=FALSE])
colnames(df)<-c("Status")
df$Status<-gsub("NONREPR", "Sterile", df$Status)
df$Status<-gsub("REPR", "Reproductive", df$Status)

pheatmap(mat, annotation_col=df,
         show_rownames = F,
         fontsize = 16,
         annotation_colors = my_colors)
#graph2ppt(file="heatmap_100topdiffgenes_reprnonreprW_abd.pptx")

# heatmap of top differential genes in head of repr vs nonrepr workers
mat = assay(rld)[ topdiff_headreprw, grepl("KOP",colnames(assay(rld)))]
mat = mat - rowMeans(mat)
df = as.data.frame(colData(rld)[,c("reprstatus1"),drop=FALSE])
colnames(df)<-c("Status")
df$Status<-gsub("NONREPR", "Sterile", df$Status)
df$Status<-gsub("REPR", "Reproductive", df$Status)

pheatmap(mat, annotation_col=df,
         show_rownames = F,
         fontsize = 16,
         annotation_colors = my_colors)
#graph2ppt(file="heatmap_100topdiffgenes_reprnonreprW_head.pptx")

#################################################################################
#################################################################################

# plot expression levels of some genes

# ABD
# most downregulated gene in Repro vs Nonrepro
topGene = rownames(res)[which.min(res$log2FoldChange)] 
topGene #LOC100650745 hexamerin
data = plotCounts(dds, gene=topGene, intgroup=c("tissue_reprstatus1"), returnData=TRUE)
ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) + 
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + 
  xlab("Repro Status and Tissue") + ylab("Normalized read count") + 
  ggtitle(topGene) + scale_y_log10() +
  scale_fill_manual(name ="",
                    values = c("blue", "lightblue","violetred","violet"),
                    labels = c("ABD.NONREPR" = "Abdomen Non-Repro", 
                               "ABD.REPR" = "Abdomen Repro", 
                               "HEAD.NONREPR" = "Head Non-Repro", 
                               "HEAD.REPR"="Head Repro"),
                    limits=c("ABD.NONREPR","HEAD.NONREPR","ABD.REPR","HEAD.REPR"))+
  theme(axis.text.x = element_blank())

#graph2ppt(file="boxplot_most_downreg_gene.pptx")

# most upregulated in Repro vs nonrepro
topGene = rownames(res)[which.max(res$log2FoldChange)] 
topGene #LOC105665921 uncharacterized
data = plotCounts(dds, gene=topGene, intgroup=c("tissue_reprstatus1"), returnData=TRUE)
ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) +  #geom_jitter(alpha = I(1/2), position=position_dodge(width=0.5), aes(color = reprstatus1)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + 
  xlab("Repro Status and Tissue") + ylab("Normalized read count") + 
  ggtitle(topGene) + scale_y_log10() +
  scale_fill_manual(name ="",
                    values = c("blue", "lightblue","violetred","violet"),
                    labels = c("ABD.NONREPR" = "Abdomen Sterile", 
                               "ABD.REPR" = "Abdomen Repro", 
                               "HEAD.NONREPR" = "Head Sterile", 
                               "HEAD.REPR"="Head Repro"),
                    limits=c("ABD.NONREPR","HEAD.NONREPR","ABD.REPR","HEAD.REPR"))+
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30),
        legend.text=element_text(size=26),
        axis.text.x = element_blank(),
        title=element_text(size=30))







#ABD

#graph2ppt(file="boxplot_most_upreg_gene.pptx")

# 12 most downregulated genes in repro vs nonrepro
n=12
selGenes = head(rownames(res)[order(res$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("tissue_reprstatus1"), returnData=TRUE))))
ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) + # geom_jitter(alpha = I(1/2), position=position_dodge(width=0.5), aes(color = tissue)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status and Tissue") + ylab("Normalized read count") + 
  scale_y_log10() + ggtitle("Top Downregulated") +
  scale_fill_manual(name ="",
                    values = c("#00CC66", "#66FFCC","#0066CC","#66CCFF"),
                    labels = c("ABD.NONREPR" = "Abdomen Sterile", 
                               "ABD.REPR" = "Abdomen Repro", 
                               "HEAD.NONREPR" = "Head Sterile", 
                               "HEAD.REPR"="Head Repro"),
                    limits=c("ABD.NONREPR","HEAD.NONREPR","ABD.REPR","HEAD.REPR"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))

#graph2ppt(file="boxplot_12most_downreg_gene.pptx")


# 12 most upregulated genes in repro vs nonrepro
n=12
selGenes = head(rownames(res)[order(-res$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("tissue_reprstatus1"), returnData=TRUE))))
ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) + # geom_jitter(alpha = I(1/2), position=position_dodge(width=0.5), aes(color = tissue)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status and Tissue") + ylab("Normalized read count") + 
  scale_y_log10() + ggtitle("Top Upregulated") +
  scale_fill_manual(name ="",
                    values = c("#00CC66", "#66FFCC","#0066CC","#66CCFF"),
                    labels = c("ABD.NONREPR" = "Abdomen Sterile", 
                               "ABD.REPR" = "Abdomen Repro", 
                               "HEAD.NONREPR" = "Head Sterile", 
                               "HEAD.REPR"="Head Repro"),
                    limits=c("ABD.NONREPR","HEAD.NONREPR","ABD.REPR","HEAD.REPR"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))

#graph2ppt(file="boxplot_12most_upreg_gene.pptx")




# Head

# 12 most downregulated genes in repro vs nonrepro
n=12
selGenes = head(rownames(res1)[order(res$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("tissue_reprstatus1"), returnData=TRUE))))
ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) + # geom_jitter(alpha = I(1/2), position=position_dodge(width=0.5), aes(color = tissue)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status and Tissue") + ylab("Normalized read count") + 
  scale_y_log10() + ggtitle("Top Downregulated") +
  scale_fill_manual(name ="",
                    values = c("#00CC66", "#66FFCC","#0066CC","#66CCFF"),
                    labels = c("ABD.NONREPR" = "Abdomen Sterile", 
                               "ABD.REPR" = "Abdomen Repro", 
                               "HEAD.NONREPR" = "Head Sterile", 
                               "HEAD.REPR"="Head Repro"),
                    limits=c("ABD.NONREPR","HEAD.NONREPR","ABD.REPR","HEAD.REPR"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))

#graph2ppt(file="boxplot_12most_downreg_gene.pptx")


# 12 most upregulated genes in repro vs nonrepro
n=12
selGenes = head(rownames(res1)[order(-res$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("tissue_reprstatus1"), returnData=TRUE))))
ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) + # geom_jitter(alpha = I(1/2), position=position_dodge(width=0.5), aes(color = tissue)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status and Tissue") + ylab("Normalized read count") + 
  scale_y_log10() + ggtitle("Top Upregulated") +
  scale_fill_manual(name ="",
                    values = c("#00CC66", "#66FFCC","#0066CC","#66CCFF"),
                    labels = c("ABD.NONREPR" = "Abdomen Sterile", 
                               "ABD.REPR" = "Abdomen Repro", 
                               "HEAD.NONREPR" = "Head Sterile", 
                               "HEAD.REPR"="Head Repro"),
                    limits=c("ABD.NONREPR","HEAD.NONREPR","ABD.REPR","HEAD.REPR"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))







####Plotting 12 for head without abd information

n=3
res1sig <- subset(res1, padj > 0.05)
selGenes = head(rownames(res1sig)[order(res1sig$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("tissue_reprstatus1"), returnData=TRUE))))
data <- subset(data, !tissue_reprstatus1 == "ABD.NONREPR" & !tissue_reprstatus1 == "ABD.REPR" )
ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) + # geom_jitter(alpha = I(1/2), position=position_dodge(width=0.5), aes(color = tissue)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status") + ylab("Normalized read count") + 
  scale_y_log10() + ggtitle("Top Upregulated in Sterile Workers: Head") +
  scale_fill_manual(name ="",
                    values = c( "#66FFCC","#66CCFF"),
                    labels = c("HEAD.NONREPR" = "Sterile", 
                               "HEAD.REPR"="Reproductive"),
                    limits=c("HEAD.NONREPR","HEAD.REPR"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=20))


# 12 most upregulated genes in repro vs nonrepro
n=3
res1sig <- subset(res1, padj > 0.05)
selGenes = head(rownames(res1sig)[order(-res1sig$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("tissue_reprstatus1"), returnData=TRUE))))
data <- subset(data, !tissue_reprstatus1 == "ABD.NONREPR" & !tissue_reprstatus1 == "ABD.REPR" )
ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) + # geom_jitter(alpha = I(1/2), position=position_dodge(width=0.5), aes(color = tissue)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status and Tissue") + ylab("Normalized read count") + 
  scale_y_log10() + ggtitle("Top Upregulated in Reproductive Workers: Head") +
  scale_fill_manual(name ="",
                    values = c( "#66FFCC","#66CCFF"),
                    labels = c("HEAD.NONREPR" = "Sterile", 
                               "HEAD.REPR"="Reproductive"),
                    limits=c("HEAD.NONREPR","HEAD.REPR"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))




####Plotting 12 for abd without head information

n=3
res1sig <- subset(res, padj > 0.05)
selGenes = head(rownames(res1sig)[order(res1sig$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("tissue_reprstatus1"), returnData=TRUE))))
data <- subset(data, !tissue_reprstatus1 == "HEAD.NONREPR" & !tissue_reprstatus1 == "HEAD.REPR" )
ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) + # geom_jitter(alpha = I(1/2), position=position_dodge(width=0.5), aes(color = tissue)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status") + ylab("Normalized read count") + 
  scale_y_log10() + ggtitle("Top Upregulated in Sterile Workers: Abdomen") +
  scale_fill_manual(name ="",
                    values = c( "#00CC66","#0066CC"),
                    labels = c("ABD.NONREPR" = "Sterile", 
                               "ABD.REPR"="Reproductive"),
                    limits=c("ABD.NONREPR","ABD.REPR"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=20))


# 12 most upregulated genes in repro vs nonrepro
n=3
res1sig <- subset(res, padj > 0.05)
selGenes = head(rownames(res1sig)[order(-res1sig$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("tissue_reprstatus1"), returnData=TRUE))))
data <- subset(data, !tissue_reprstatus1 == "HEAD.NONREPR" & !tissue_reprstatus1 == "HEAD.REPR" )

ggplot(data, aes(x=tissue_reprstatus1, y=count, fill=tissue_reprstatus1)) + # geom_jitter(alpha = I(1/2), position=position_dodge(width=0.5), aes(color = tissue)) +
  geom_boxplot(outlier.color=NA, position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status and Tissue") + ylab("Normalized read count") + 
  scale_y_log10() + ggtitle("Top Upregulated in Reproductive Workers: Abdomen") +
  scale_fill_manual(name ="",
                    values = c( "#00CC66","#0066CC"),
                    labels = c("ABD.NONREPR" = "Sterile", 
                               "ABD.REPR"="Reproductive"),
                    limits=c("ABD.NONREPR","ABD.REPR"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))








###################################################
### Get fpkm values (normalised FPM's)
###################################################

head(res)
signif_res_abd<-subset(res, padj<0.05)
signif_res_abd<- as.data.frame(signif_res_abd)
signif_res_abd$gene<-rownames(signif_res_abd)
nrow(signif_res_abd) #8097

head(res1)
signif_res_head<-subset(res1, padj<0.05)
signif_res_head<- as.data.frame(signif_res_head)
signif_res_head$gene<-rownames(signif_res_head)
nrow(signif_res_head) #279

# Import gene length of every gene in Bter_1.0 
Bter_gene_length <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Transcription/diff_exp_all_genes/Bter_gene_length.csv")
Bter_gene_length<-Bter_gene_length[,c(2,6)]
Bter_gene_length<-as.data.frame(Bter_gene_length)
rownames(Bter_gene_length)<-Bter_gene_length$gene
Bter_gene_length <- Bter_gene_length[ order(row.names(Bter_gene_length)), ]

# Add gene length to the correct genes
Bter_gene_length <- Bter_gene_length[match(rownames(dds), Bter_gene_length$gene),]
Bter_gene_length <- Bter_gene_length[,2]
mcols(dds)$basepairs<-Bter_gene_length

# Obtain fpkm values (not all genes had a length imported so some lost)
fpkm_values<-fpkm(dds)
fpkm_values_with_na<-apply(fpkm_values, 1, function(x){any(is.na(x))})
fpkm_values <- fpkm_values[!fpkm_values_with_na,]
nrow(fpkm_values) #8232 genes left

fpkm_values<-as.data.frame(fpkm_values)
fpkm_values$gene<-rownames(fpkm_values)
rownames(fpkm_values)<-c()
fpkm_values$abd_repro_fpkm_mean<-apply(fpkm_values[,c(1,3,5,7,25,27,29,31,33,35,37,39,57,59,61,63)],1,mean)
fpkm_values$abd_nonrepro_fpkm_mean<-apply(fpkm_values[,c(9,11,13,15,17,19,21,23,41,43,45,47,49,51,53,55)],1,mean)
fpkm_values$head_repro_fpkm_mean<-apply(fpkm_values[,c(2,4,6,8,26,28,30,32,34,36,38,40,58,60,62,64)],1,mean)
fpkm_values$head_nonrepro_fpkm_mean<-apply(fpkm_values[,c(10,12,14,16,18,20,22,24,42,44,46,48,50,52,54,56)],1,mean)
write.csv(as.data.frame(fpkm_values), file="all_genes_fpkm_values.csv")

fpkm_diff_exp_genes<-merge(signif_res_abd, fpkm_values, by = "gene")
nrow(fpkm_diff_exp_genes)#6858
head(fpkm_diff_exp_genes)
write.csv(as.data.frame(fpkm_diff_exp_genes), file="abd_diff_exp_genes_fpkm_values.csv")

fpkm_diff_exp_genes<-merge(signif_res_head, fpkm_values, by = "gene")
nrow(fpkm_diff_exp_genes)#225
head(fpkm_diff_exp_genes)
write.csv(as.data.frame(fpkm_diff_exp_genes), file="head_diff_exp_genes_fpkm_values.csv")

## Should also write out files with just the gene names only (LOC id) for the diff exp genes
## and also another file for upreg and one for downreg per tissue (6 files total), to make GO analysis
## easier. If not: on command line: cut -f 2 -d ',' <file> > <output> then sed -i 's/"//g' <file>