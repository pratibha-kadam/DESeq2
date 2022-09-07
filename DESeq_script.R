
#upload DESeq2 and ggplot2 library

library(DESeq2)
library(ggplot2)

#read the count data file

Counts <- read.csv("IHF_34_gene_count.csv", header = TRUE, row.names = 1, sep = ',')
Counts 

#Remove the lower expressed data whose sum in all the samples are less than 30 (you can take any number)
Counts <- Counts[which(rowSums(Counts) > 30),]

#Provide the conditions for each sample
condition <- factor(c("Control","T1","T1","T1","T2","T2" ))

#Assign the conditions to the columns
coldata <- data.frame(row.names = colnames(Counts),condition)
coldata


#Create DESeq dataset from matrix
dds <- DESeqDataSetFromMatrix(countData = Counts,colData = coldata, design = ~condition)

dds

##define the Control sasmple as reference
dds$condition <- relevel(dds$condition, ref = "Control")
dds$condition

# Run DESeq
dds <- DESeq(dds)

#plot despersion plot
plotDispEsts(dds)

#Check differential expression for specific condition from the data
res <- results(dds, contrast = c("condition", "T1", "Control"))
res
