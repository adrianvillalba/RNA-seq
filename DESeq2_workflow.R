## Load DESeq2, stringr and ggplot2 libraries 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
library(stringr)
library(ggplot2)

 
## First set the working directory
setwd("/home/adri/Desktop/RNA/files/raw_data/Total_bam")
getwd()

## Read the output from featurecounts (counts.txt here) and  save the object. 
counts=read.csv("counts.txt", sep="", head=T, skip=1, row.names = "Geneid")
View(counts)

## Create an object for storing names of samples and condition of the samples (Basal/ConditionX) 
colnames(counts)[6:11] <- c("A", "B", "C", "Ax", "Bx", "Cx")
colnames(counts)[6:11]
countData <- counts[6:11]
View(countData)


samples=cbind(colnames(counts)[6:11])
meta_data <- matrix(samples)
metadata <- cbind(meta_data,c("Basal","Basal","Basal","Condition","Condition","Condition"))
View(metadata)

## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metadata, 
                              design=~V2)

# Run the DESeq model
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE)) 
summary(res)
write.csv(res, file = "Gene_list.csv",row.names=TRUE)
## Just in case: sort by adj p-value:
res <- res[order(res$padj),]
head(res)


# Plot a Volcano plot
## Reset par
par(mfrow=c(1,1))
## Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

## Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
