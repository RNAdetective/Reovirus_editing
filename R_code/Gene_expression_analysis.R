library(DESeq2)
library(tidyverse)
#read in data
Counts <- as.matrix (read.csv("gene_count_matrix.csv", row.names = "gene_id"))

Coldata <- read.csv("Coldata.csv", row.names = 1)
#run DeSeq
dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = Coldata,
                              design= ~ Condition)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#ReoV vs Control
res1 <- results(dds, contrast=c("Condition","ReoV.Infection","Control"))
# Remove rows with missing values
res1 <- res1[complete.cases(res1), ]  
#Filter for P_adj =< 0.05, logFold2 
res1 <- res1[res1$padj <= 0.05 & res1$log2FoldChange >= 1, ]

write.csv(res1, "ReoVvsControlexpression.csv")
#muReoV VS Control
res2 <- results(dds, contrast=c("Condition","muReoV.Infection", "Control"))
#Remove rows with missing values
res2 <- res2[complete.cases(res2), ]  
#Filter for P_adj =< 0.05, logFold2 
res2 <- res2[res2$padj <= 0.05 & res2$log2FoldChange >= 1, ]

write.csv(res2, "muReoVvsControlexpression.csv")

#For transcript_count_matrix

Trans_Counts <- as.matrix (read.csv("transcript_count_matrix.csv", row.names = "Transcript_id"))
Coldata <- read.csv("Coldata.csv", row.names = 1)

#run DeSeq
dds <- DESeqDataSetFromMatrix(countData = Trans_Counts,
                              colData = Coldata,
                              design= ~ Condition)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#ReoV vs Control
Trans_counts_1 <- results(dds, contrast=c("Condition","ReoV.Infection","Control"))
# Remove rows with missing values
Trans_counts_1<- Trans_counts_1[complete.cases(Trans_counts_1), ]  
#Filter for P_adj =< 0.05, logFold2 
Trans_counts_1 <- Trans_counts_1[Trans_counts_1$padj <= 0.05 & Trans_counts_1$log2FoldChange >= 1, ]

write.csv(Trans_counts_1, "ReoVvsControl_Transcript_Count.csv")
#muReoV VS Control
Trans_counts_2 <- results(dds, contrast=c("Condition","muReoV.Infection", "Control"))
#Remove rows with missing values
Trans_counts_2  <- Trans_counts_2 [complete.cases(Trans_counts_2 ), ]  
#Filter for Significant DEGs P_adj =< 0.05, logFold2 
Trans_counts_2  <- Trans_counts_2 [Trans_counts_2 $padj <= 0.05 & Trans_counts_2 $log2FoldChange >= 1, ]

write.csv(Trans_counts_2 , "muReoVvsControl_Trnascript_Count.csv")

#Annotate Gene Ids
library(annotables)
gene_names <- grcm38 %>% filter(ensgene %in% ReoVvsControlexpression$Gene_id)
names(gene_names)[1] <- "Gene_id"
#join annotable and Deseq2 results
ReoVvsControlexpression <- left_join(ReoVvsControlexpression, gene_names, by= "Gene_id")
write.csv(ReoVvsControlexpression, "ReoVvsControlexpression.csv")

#Repeat for other data frames

gene_names <- grcm38 %>% filter(ensgene %in% muReoVvsControlexpression$Gene_id)
names(gene_names)[1] <- "Gene_id"
muReoVvsControlexpression <- left_join(muReoVvsControlexpression, gene_names, by= "Gene_id")
 write.csv(muReoVvsControlexpression, "muReoVvsControlexpression.csv")




















