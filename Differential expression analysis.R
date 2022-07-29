setwd("/Users/gsetia1/Desktop")

library(DESeq2)
library(ggplot2)

Counts = read.delim("count_table.csv", header = TRUE, row.names = 1, sep = ",")

Counts

#pre processing-- Filter out those which has very few mapped reads i.e less than 50
Counts = Counts[which(rowSums(Counts) > 50),]

Counts

#Assign the conditions of the different samples in experiment
condition = factor(c("C","C","C","C","C","C","S","S","S","S","S","S"))

#Make dataframe of these conditions
coldata = data.frame(row.names = colnames(Counts), condition)
coldata

dds = DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)

dds = DESeq(dds)
vsdata = vst(dds, blind=FALSE)

#make PCA plot
PlotPCA(vsdata, intgroup = "condition")

#make dispersion plot: The plot shows dispersion of the data-- as the read counts get higher there is less variation-- indicated by the decreasing slope-- the line (red) should be below 1-- the lower it is-- the better
plotDispEsts(dds)

#result table
res = results(dds, contrast= c("condition", "S", "C"))

#if log2foldchain is negative-- it means that the gene is downregulated and positive means that the gene is upregulated
#Also look at the padj value to see if that gene is significant

sigs = na.omit(res)

sigs = sigs[sigs$padj < 0.05,]

sigs 

#look the number of genes that pass the significance values

write.csv(sigs, file = "deseqresults.csv")