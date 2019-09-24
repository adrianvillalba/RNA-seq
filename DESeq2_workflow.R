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
#Note this is the important file for later enrichment.

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


##Gene enrichment script:
#Install patfindR packages and fix execution issues
install.packages("pathfindR")

detach("package:pathfindR", unload=TRUE)
detach("package:pathview", unload=TRUE)

library(pathfindR) 
library(pathview)

install.packages("tidyverse")
library(tidyverse)

source('https://bioconductor.org/biocLite.R')
biocLite('org.Hs.eg.db')
library('org.Hs.eg.db')

#Prepare your working directory
setwd(choose.dir())
getwd()

#Tidy your data and prepare as input for patfindR
input <- read.csv2(file.choose(),sep=",")
input1<-select(input, X, log2FoldChange, padj)

real_input <- na.omit(input1)
real_input$padj <- as.numeric(as.character(real_input$padj))
real_input$log2FoldChange <- as.numeric(as.character(real_input$log2FoldChange))
real_input$X <- as.character(real_input$X)

symbols <- mapIds(org.Hs.eg.db, real_input$X, 'SYMBOL', 'ENSEMBL')
a<-cbind(symbols, real_input)

test_input<-na.omit(a)
test_input$X <- NULL
View(test_input)

#Run pathfindR function
output1 <- run_pathfindR(test_input)
print("Do not forget to manually export the bubble chart from Plots window!")
View(output1)

#Run clustering
clustered <- cluster_pathways(output1)
print("Do not forget to manually export the connectivity plot!")
enrichment_chart(clustered, plot_by_cluster = TRUE)
print("Do not forget to manually export the hierarchical plot!")

#Pathway scoring
## Select the representative pathways
pws_table <- clustered[clustered$Status == "Representative", ]

## Load the differential expression matrix
diff_matrix <- file.choose()

## Vector of x-axis labe
cases <- c("label1", "label2")

## Calculate pathway scores and plot heatmap
score_matrix <- calculate_pw_scores(pws_table, diff_matrix, cases)

#Finally, do not forget to export the heatmap!
print("Do not forget to manually export the hierarchical plot!")
