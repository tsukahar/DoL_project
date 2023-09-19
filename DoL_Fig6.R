#Now I want make boxplots for the effect of MM401 on Bic responsive and non-responsive genes.
#Load DEseq2
library(DESeq2)

#Load forestmngr to do rounding of numbers
library(forestmangr)

#import  count table 
MERGED_BRUSEQ_WHOLEGENE_MM10_FEATURECOUNTS_OUTPUT_BT2 <- read.delim("/Users/hioutsukahara/Dropbox_UM/DoL_project/BrU-seq/MERGED_BRUSEQ_WHOLEGENE_MM10_FEATURECOUNTS_OUTPUT_BT2")

#remove unessesary columns 
Count_Table <- MERGED_BRUSEQ_WHOLEGENE_MM10_FEATURECOUNTS_OUTPUT_BT2[,-c(2,3,4,5,6)]
rownames(Count_Table) <- Count_Table$Geneid
Count_Table$Geneid <- NULL
colnames(Count_Table)

#rename columns
colnames(Count_Table) <- c("DMSO_1h_Bic_rep1","DMSO_1h_Bic_rep2","DMSO_1h_Bic_rep3",
                           "DMSO_4h_Bic_rep1","DMSO_4h_Bic_rep2","DMSO_4h_Bic_rep3",
                           "DMSO_0h_Bic_rep1","DMSO_0h_Bic_rep2","DMSO_0h_Bic_rep3",
                           "MM401_1h_Bic_rep1","MM401_1h_Bic_rep2","MM401_1h_Bic_rep3",
                           "MM401_4h_Bic_rep1","MM401_4h_Bic_rep2","MM401_4h_Bic_rep3",
                           "MM401_1h_water_rep1","MM401_1h_water_rep2","MM401_1h_water_rep3",
                           "MM401_4h_water_rep1","MM401_4h_water_rep2","MM401_4h_water_rep3")


#round numbers and make gene_id as rownames
Count_Table <- round_df(Count_Table, digits = 0)
View(Count_Table)

#Reorder the columns
colnames(Count_Table)
Count_Table <- Count_Table[c("DMSO_0h_Bic_rep1", "DMSO_0h_Bic_rep2", "DMSO_0h_Bic_rep3",
                             
                             "DMSO_1h_Bic_rep1", "DMSO_1h_Bic_rep2", "DMSO_1h_Bic_rep3",
                             
                             "DMSO_4h_Bic_rep1", "DMSO_4h_Bic_rep2", "DMSO_4h_Bic_rep3",
                             
                             "MM401_1h_Bic_rep1", "MM401_1h_Bic_rep2", "MM401_1h_Bic_rep3",
                             
                             "MM401_4h_Bic_rep1", "MM401_4h_Bic_rep2", "MM401_4h_Bic_rep3",
                             
                             "MM401_1h_water_rep1","MM401_1h_water_rep2","MM401_1h_water_rep3",
                             
                             "MM401_4h_water_rep1","MM401_4h_water_rep2","MM401_4h_water_rep3")]



#remove MM401_4h_water_rep3. This sample has poor quality.
Count_Table[,"MM401_4h_water_rep3"]<-NULL

#generate experimental design file so called “se”. 
Design_Bru_1216 <- data.frame(row.names = colnames(Count_Table), 
                              
                              Conditions = c("DMSO_0h","DMSO_0h","DMSO_0h","DMSO_1h","DMSO_1h","DMSO_1h","DMSO_4h","DMSO_4h","DMSO_4h",
                                             
                                             "MM401_1h","MM401_1h","MM401_1h","MM401_4h","MM401_4h","MM401_4h",
                                             
                                             "MM401_W1h","MM401_W1h","MM401_W1h","MM401_W4h","MM401_W4h"))

#check the design
Design_Bru_1216

#Run DEseq2 
dds_1216 <- DESeqDataSetFromMatrix(Count_Table, colData=Design_Bru_1216, design= ~ Conditions)
dds_1216_run <- DESeq(dds_1216, betaPrior=TRUE)
resultsNames(dds_1216_run)

#make results table for 4hBic vs 0hbic and 1hBic vs 0h Bic
##Bic-regulated genes in DMSO and MM401
Restable_4hrBic_vs_0hrBic_DMSO <- results(dds_1216_run, contrast=c("Conditions", "DMSO_4h", "DMSO_0h"), alpha=0.05)
Restable_1hrBic_vs_0hrBic_DMSO <- results(dds_1216_run, contrast=c("Conditions", "DMSO_1h", "DMSO_0h"), alpha=0.05)
Restable_4hrMM401Bic_vs_0hrBic_DMSO <- results(dds_1216_run, contrast=c("Conditions", "MM401_4h", "DMSO_0h"), alpha=0.05)
Restable_1hrMM401Bic_vs_0hrBic_DMSO <- results(dds_1216_run, contrast=c("Conditions", "MM401_1h", "DMSO_0h"), alpha=0.05)

#write tables
write.table(Restable_4hrBic_vs_0hrBic_DMSO,file='Restable_4hrBic_vs_0hrBic_DMSO.txt',sep="\t")
write.table(Restable_1hrBic_vs_0hrBic_DMSO,file='Restable_1hrBic_vs_0hrBic_DMSO.txt',sep="\t")
write.table(Restable_4hrMM401Bic_vs_0hrBic_DMSO,file='Restable_4hrMM401Bic_vs_0hrBic_DMSO.txt',sep="\t")
write.table(Restable_1hrMM401Bic_vs_0hrBic_DMSO,file='Restable_1hrMM401Bic_vs_0hrBic_DMSO.txt',sep="\t")


#load restult tables
Restable_4hrBic_vs_0hrBic_DMSO<-read.delim("/Users/hioutsukahara/Restable_4hrBic_vs_0hrBic_DMSO.txt")
Restable_1hrBic_vs_0hrBic_DMSO<-read.delim("/Users/hioutsukahara/Restable_1hrBic_vs_0hrBic_DMSO.txt")
Restable_4hrMM401Bic_vs_0hrBic_DMSO<-read.delim("/Users/hioutsukahara/Restable_4hrMM401Bic_vs_0hrBic_DMSO.txt")
Restable_1hrMM401Bic_vs_0hrBic_DMSO<-read.delim("/Users/hioutsukahara/Restable_1hrMM401Bic_vs_0hrBic_DMSO.txt")



#combine files
Bic_log2fc<-cbind(Count_Table,Restable_1hrBic_vs_0hrBic_DMSO,Restable_1hrMM401Bic_vs_0hrBic_DMSO,Restable_4hrBic_vs_0hrBic_DMSO,Restable_4hrMM401Bic_vs_0hrBic_DMSO)
View(Bic_log2fc)
colnames(Bic_log2fc)

#remove unnecessary columns
Bic_log2fc<-Bic_log2fc[,c(22,28,34,40,26,32,38,44)]
#rename columns
colnames(Bic_log2fc) <- c("log2FoldChange_1hr", "log2FoldChange_MM401_1hr","log2FoldChange_4hr","log2FoldChange_MM401_4hr","padj_1hr","padj_MM401_1hr","padj_4hr","padj_MM401_4hr")
colnames(Bic_log2fc)
#Select All Bic-regulated genes. Bicuculline responsive genes are defined by padj<0.05
Bic_RGenes <- subset(Bic_log2fc, padj_1hr < 0.05 | padj_4hr < 0.05)

write.table(Bic_RGenes, file='Bic_RGenes.txt', sep='\t')
View(Bic_RGenes)

#Divide Bic response genes into up and down regulated genes
up <-subset(Bic_RGenes, log2FoldChange_1hr>0 & log2FoldChange_4hr>0)
down <-subset(Bic_RGenes, log2FoldChange_1hr<0 & log2FoldChange_4hr<0)


#divide these gene into early response genes (ERGs) and late response genes (LRGs)
upERGs <-subset(up, log2FoldChange_1hr>log2FoldChange_4hr)
upLRGs <-subset(up, log2FoldChange_4hr>log2FoldChange_1hr)
downERGs <-subset(down, log2FoldChange_1hr<log2FoldChange_4hr)
downLRGs <-subset(down, log2FoldChange_4hr<log2FoldChange_1hr)

#make box plots for the effect of MM401 on Bic-responsive genes
#ggplot2 ref: https://appsilon.com/ggplot2-boxplots/
library(ggplot2)
require(reshape2)

ggplot(data = melt(upERGlog2fc<-upERGs[,(1:4)]), aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(width=0.8)+ ylim(-0.5,2)
coord_fixed(1.5) +
  ylab('Log2FC') +
  geom_hline(yintercept=0, linetype='dashed', color='#9b9ba3') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

ggplot(data = melt(upLRGlog2fc<-upLRGs[,(1:4)]), aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(width=0.8)+ ylim(-0.5,2)
coord_fixed(1.5) +
  ylab('Log2FC') +
  geom_hline(yintercept=0, linetype='dashed', color='#9b9ba3') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = melt(downERGlog2fc<-downERGs[,(1:4)]), aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(width=0.8)+ ylim(-1,0.5)
coord_fixed(1.5) +
  ylab('Log2FC') +
  geom_hline(yintercept=0, linetype='dashed', color='#9b9ba3') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

ggplot(data = melt(downLRGlog2fc<-downsLRG[,(1:4)]), aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(width=0.8)+ ylim(-1,0.5)
coord_fixed(1.5) +
  ylab('Log2FC') +
  geom_hline(yintercept=0, linetype='dashed', color='#9b9ba3') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
 

#Generate data.frame for Bic non-responsive genes
#remove Bic response genes
Bic_nonRGlog2fc<-Bic_log2fc[! rownames(Bic_log2fc)%in%rownames(Bic_RGenes),]
View(Bic_nonRGlog2fc)
#remove rows with NA
Bic_nonRGlog2fc<-Bic_nonRGlog2fc[complete.cases(Bic_nonRGlog2fc), ]

#make box plots
ggplot(data = melt(Bic_nonRGlog2fc), aes(x=variable, y=value, fill=variable)) + 
geom_boxplot(width=0.8)+ ylim(-1,1)
coord_fixed(1.5) +
ylab('Log2FC') +
geom_hline(yintercept=0, linetype='dashed', color='#9b9ba3') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 