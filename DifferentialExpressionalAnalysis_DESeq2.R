# DifferentialExpressionAnalysis: using DESeq2 (count based method) that fits a negative binomial glm (generalized linear model) and estimates for dispersion

##################################
## load files
##################################
load("/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/SRP016059_SummarizedOverlap.rda")
#rowRanges(SRP016059_SummarizedOverlap)
#colData(SRP016059_SummarizedOverlap)
#dim(assay(SRP016059_SummarizedOverlap))
#head(assay(SRP016059_SummarizedOverlap))

##################################
## Differential Gene Expression Analysis
##################################
library(DESeq2)
# The last variable in the design is used by default for building results tables. make sure the "control" or "untreated" level is the first level, such that log fold changes will be treated over control, and not control over treated.
levels(SRP016059_SummarizedOverlap$experiment_title)
# SRP016059_SummarizedOverlap$experiment_title <- as.factor(SRP016059_SummarizedOverlap$experiment_title)
SRP016059_SummarizedOverlap$experiment_title <- relevel(SRP016059_SummarizedOverlap$experiment_title, ref="normal")
levels(SRP016059_SummarizedOverlap$experiment_title)

# Creating a DESeqDataSet object.
# The *DESeqDataSet* object is just an extension of the *SummarizedExperiment* object, with a few changes. The matrix in `assay` is now accessed with `counts` and the elements of this matrix are required to be non-negative integers (0,1,2,...).
SRP016059_dds <- DESeqDataSet(SRP016059_SummarizedOverlap, design= ~ experiment_title)
class(SRP016059_dds)
#run the *DESeq2* model. 
SRP016059_dds <- DESeq(SRP016059_dds)
# build a results table, which by default will compare the levels in the last variable in the design
SRP016059_results_0.1 <- results(SRP016059_dds)

### Examining results tables
table(SRP016059_results_0.1$padj < 0.1)
summary(SRP016059_results_0.1,alpha=0.1)
## use threshold of 0.05 for adjusted Pvalue
SRP016059_results <- results(SRP016059_dds,alpha=0.05)
table(SRP016059_results$padj < 0.05)
summary(SRP016059_results,alpha=0.05)
SRP016059_results
dim(SRP016059_results)
SRP016059_results_final <- subset(SRP016059_results,padj<0.05,)
dim(SRP016059_results_final)
SRP016059_results_final_sorted <- SRP016059_results_final[order(-abs(SRP016059_results_final$log2FoldChange)),]
save(SRP016059_results_final_sorted,file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/SRP016059_results_final_sorted.rda")

# use threshold of 0.05 for adjusted Pvalue AND Foldchange greater than 2
SRP016059_results_2fold <- results(SRP016059_dds,alpha=0.05,lfcThreshold=1)
table(SRP016059_results_2fold$padj < 0.05)
summary(SRP016059_results_2fold,alpha=0.05)
# use threshold of 0.05 for adjusted Pvalue AND Foldchange greater than 1.5
SRP016059_results_1.5fold <- results(SRP016059_dds,alpha=0.05,lfcThreshold=0.585)
table(SRP016059_results_1.5fold$padj < 0.05)
summary(SRP016059_results_1.5fold,alpha=0.05)
SRP016059_results_1.5fold <- subset(SRP016059_results_1.5fold,padj<0.05,abs(log2FoldChange)>0.585)
SRP016059_results_1.5fold_sorted <- SRP016059_results_1.5fold[order(-abs(SRP016059_results_1.5fold$log2FoldChange)),]
save(SRP016059_results_1.5fold_sorted,file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/SRP016059_results_1.5fold_sorted.rda")
##################################
### Visualizing results
##################################
# MA plot
summary(SRP016059_results$log2FoldChange)
par(mfrow=c(2,2))
plotMA(SRP016059_results_0.1,main=c("alpha=0.1","lfcThreshold=0"),ylim=c(-21,15))
plotMA(SRP016059_results,alpha=0.05,main=c("alpha=0.05","lfcThreshold=0"),ylim=c(-21,15))
plotMA(SRP016059_results_1.5fold,alpha=0.05,main=c("alpha=0.05","lfcThreshold=1.5"),ylim=c(-21,15))
plotMA(SRP016059_results_2fold,alpha=0.05,main=c("alpha=0.05","lfcThreshold=2"),ylim=c(-21,15))

# plot dispersion eastimates
plotDispEsts(SRP016059_dds)
#A p-value histogram:
hist(SRP016059_results$pvalue[SRP016059_results$baseMean > 1], 
     col="grey", border="white", xlab="", ylab="", main="")
# Examine the counts for the top gene, sorting by p-value:
plotCounts(SRP016059_dds, gene=which.min(SRP016059_results$padj), intgroup="experiment_title")
# gene with max (up) log fold change
SRP016059_results_final_sorted[which.max(SRP016059_results_final_sorted$log2FoldChange),]
plotCounts(SRP016059_dds, gene=which.max(SRP016059_results_final_sorted$log2FoldChange), intgroup="experiment_title")
# gene with min (up) log fold change
SRP016059_results_final_sorted[which.min(SRP016059_results_final_sorted$log2FoldChange),]
plotCounts(SRP016059_dds, gene=which.min(SRP016059_results_final_sorted$log2FoldChange), intgroup="experiment_title")

##################################
# annotations
##################################
# check & split up the rownames of the results object, which contain ENSEMBL gene ids, separated by the plus sign, +.
grep("\\+",rownames(SRP016059_results_final_sorted))
library("AnnotationDbi")
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")
# This is the organism annotation package (“org”) for Homo sapiens (“Hs”), organized as an AnnotationDbi database package (“db”), using Entrez Gene IDs (“eg”) as primary key. 
# get a list of all available key types
columns(org.Hs.eg.db)
# use the mapIds function to add individual columns to our results table. We provide the row names of our results table as a key, and specify that keytype=ENSEMBL. The column argument tells the mapIds function which information we want, and the multiVals argument tells the function what to do if there are multiple possible values for a single input value.
SRP016059_results_final_sorted$symbol <- mapIds(org.Hs.eg.db,
                                                keys=row.names(SRP016059_results_final_sorted),
                                                column="SYMBOL",
                                                keytype="ENSEMBL",
                                                multiVals="first")
SRP016059_results_final_sorted$entrez <- mapIds(org.Hs.eg.db,
                                                keys=row.names(SRP016059_results_final_sorted),
                                                column="ENTREZID",
                                                keytype="ENSEMBL",
                                                multiVals="first")
SRP016059_results_final_sorted$entrez <- mapIds(org.Hs.eg.db,
                                                keys=row.names(SRP016059_results_final_sorted),
                                                column="GENENAME",
                                                keytype="ENSEMBL",
                                                multiVals="first")
SRP016059_results_final_sorted$entrez <- mapIds(org.Hs.eg.db,
                                                keys=row.names(SRP016059_results_final_sorted),
                                                column="ENTREZID",
                                                keytype="ENSEMBL",
                                                multiVals="first")

SRP016059_results_final_sorted

# Exporting results
SRP016059_results_final_sorted_DataFrame <- as.data.frame(SRP016059_results_final_sorted)
save(SRP016059_results_final_sorted_DataFrame,file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/SRP016059_results_final_sorted_DataFrame.rda")
write.csv(SRP016059_results_final_sorted_DataFrame, file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/SRP016059_results_final_sorted.csv")
