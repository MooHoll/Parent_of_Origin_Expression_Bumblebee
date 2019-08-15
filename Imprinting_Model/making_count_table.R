### Making Count Table ###
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Count_Table")
library(plyr)

# Read in the male count table, add column names and remove information not needed
male<-read.table("./read_counts/M_trimmed_02_10_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)

# Read in the queen count table for the same sample
queen<-read.table("./read_counts/Q_trimmed_02_10_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)

# Merge the queen and male information with additonal sample specific data
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"

# Remove overlapping SNP/gene merge (would be better to reslove this than remove but only lose ~130SNPs)
abd0210<-both[!grepl("ambig", both$geneID),]


########## Do above for all samples: ABD
male<-read.table("./read_counts/M_trimmed_02_12_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_12_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd0212<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_18_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_18_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd0218<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_22_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_22_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd0222<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_35_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_35_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd0235<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_37_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_37_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd0237<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_44_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_44_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd0244<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_54_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_54_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd0254<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_74_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_74_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd1274<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_75_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_75_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd1275<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_76_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_76_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd1276<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_78_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_78_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"NUR"
abd1278<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_82_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_82_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd1282<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_84_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_84_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd1284<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_85_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_85_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd1285<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_89_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_89_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd1289<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_09_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_09_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd2209<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_18_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_18_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd2218<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_21_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_21_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd2221<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_23_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_23_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd2223<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_26_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_26_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"NUR"
abd2226<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_30_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_30_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd2230<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_35_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_35_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd2235<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_49_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_49_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd2249<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_44_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_44_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd3144<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_45_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_45_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd3145<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_64_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_64_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd3164<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_65_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_65_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"FOR"
abd3165<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_79_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_79_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd3179<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_84_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_84_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd3184<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_86_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_86_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd3186<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_89_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_89_ABD_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"ABD"
both["worker_type"]<-"SUB"
abd3189<-both[!grepl("ambig", both$geneID),]

########## Do above for all samples: KOP

male<-read.table("./read_counts/M_trimmed_02_10_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_10_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP0210<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_12_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_12_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP0212<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_18_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_18_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP0218<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_22_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_22_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP0222<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_35_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_35_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP0235<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_37_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_37_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP0237<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_44_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_44_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP0244<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_02_54_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_02_54_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP0254<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_74_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_74_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP1274<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_75_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_75_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP1275<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_76_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_76_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP1276<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_78_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_78_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"NUR"
KOP1278<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_82_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_82_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP1282<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_84_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_84_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP1284<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_85_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_85_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP1285<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_12_89_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_12_89_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "initial"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP1289<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_09_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_09_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP2209<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_18_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_18_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP2218<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_21_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_21_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP2221<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_23_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_23_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP2223<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_26_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_26_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"NUR"
KOP2226<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_30_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_30_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP2230<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_35_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_35_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP2235<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_22_49_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_22_49_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"1"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP2249<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_44_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_44_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP3144<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_45_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_45_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP3145<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_64_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_64_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP3164<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_65_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_65_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"FOR"
KOP3165<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_79_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_79_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP3179<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_84_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_84_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP3184<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_86_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_86_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP3186<-both[!grepl("ambig", both$geneID),]

male<-read.table("./read_counts/M_trimmed_31_89_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(male)=c("count_male","chrom","SNP","geneID")
male<-male[,-2]
male<-subset(male, male$count>4)
queen<-read.table("./read_counts/Q_trimmed_31_89_KOP_counts_geneOnly.txt", quote="\"", comment.char="")
colnames(queen)=c("count_queen","chrom","SNP","geneID")
queen<-queen[,-2]
queen<-subset(queen, queen$count>4)
both<-merge(male,queen, by=c("SNP","geneID"), all=T)
both["direction_cross"]<- "reciproc"
both["familyID"]<-"2"
both["tissue"]<-"KOP"
both["worker_type"]<-"SUB"
KOP3189<-both[!grepl("ambig", both$geneID),]

rm(male)
rm(queen)

### PUTTING COUNT TABLES TOGETHER
test<-rbind.data.frame(abd0210,abd0212,abd0218)

count_table_02_abd<-rbind.data.frame(abd0210,abd0212,abd0218,abd0222,abd0235,abd0237,abd0244,abd0254)
cotun_table_12_abd<-rbind.data.frame(abd1274,abd1275,abd1276,abd1278,abd1282,abd1284,abd1285,abd1289)
count_table_22_abd<-rbind.data.frame(abd2209,abd2218,abd2221,abd2223,abd2226,abd2230,abd2235,abd2249)
count_table_31_abd<-rbind.data.frame(abd3144,abd3145,abd3164,abd3165,abd3179,abd3184,abd3186,abd3189)
count_table_abd<-rbind.data.frame(count_table_02_abd,cotun_table_12_abd,count_table_22_abd,count_table_31_abd)

count_table_02_kop<-rbind.data.frame(KOP0210,KOP0212,KOP0218,KOP0222,KOP0235,KOP0237,KOP0244,KOP0254)
count_table_12_kop<-rbind.data.frame(KOP1274,KOP1275,KOP1276,KOP1278,KOP1282,KOP1284,KOP1285,KOP1289)
count_table_22_kop<-rbind.data.frame(KOP2209,KOP2218,KOP2221,KOP2223,KOP2226,KOP2230,KOP2235,KOP2249)
count_table_31_kop<-rbind.data.frame(KOP3144,KOP3145,KOP3164,KOP3165,KOP3179,KOP3184,KOP3186,KOP3189)
count_table_kop<-rbind.data.frame(count_table_02_kop,count_table_12_kop,count_table_22_kop,count_table_31_kop)

count_table<-rbind.data.frame(count_table_abd,count_table_kop)
colnames(count_table)<-c("SNP_pos", "geneID","paternal_count","maternal_count","direction_cross","familyID","tissue","worker_type")

###### DID DO BELOW AS REQUIRE MIN COVERAGE OF 4 READS ABOVE
# Removing all positions from all crosses where maternal count == 0 in 1 sample 
#df<-test[!test$geneID %in% test$geneID[test$maternal_count=="0"],]
#filtered_count_table<-count_table[!count_table$geneID %in% count_table$geneID[count_table$count_queen=="0"],]

write.csv(count_table, file="count_table_holl.csv")

