# Data: SRP016059 / GSE41476 data with 3 gastric tumors and 2 gastric normal tissues. 
# Goal: to identify genes that are differentially expressed in gastric tumor tissue as compared to normal tissue
# DifferentialExpressionAnalysis: using DESeq2 that fits a negative binomial glm (generalized linear model) and estimates for dispersion

##################################
## load packages
##################################
# if a new package from bioconductor has to be installed. write following two commands
# source("https://bioconductor.org/biocLite.R")
# biocLite("package name")
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)

##################################
## Differential Gene Expression Analysis
##################################
load("/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/GSE41476_SummarizedOverlap.rda")
#rowRanges(GSE41476_SummarizedOverlap)
#colData(GSE41476_SummarizedOverlap)
dim(assay(GSE41476_SummarizedOverlap))
head(assay(GSE41476_SummarizedOverlap))

# The last variable in the design is used by default for building results tables. make sure the "control" or "untreated" level is the first level, such that log fold changes will be treated over control, and not control over treated.
GSE41476_SummarizedOverlap$experiment_title <- factor(GSE41476_SummarizedOverlap$experiment_title,levels=c("normal","tumor"))
#GSE41476_SummarizedOverlap$experiment_title <- relevel(GSE41476_SummarizedOverlap$experiment_title, ref="normal")
levels(GSE41476_SummarizedOverlap$experiment_title)

# Creating a DESeqDataSet object.
# The *DESeqDataSet* object is just an extension of the *SummarizedExperiment* object, with a few changes. The matrix in `assay` is now accessed with `counts` and the elements of this matrix are required to be non-negative integers (0,1,2,...).
GSE41476_dds <- DESeqDataSet(GSE41476_SummarizedOverlap, design= ~ experiment_title)
class(GSE41476_dds)
#run the *DESeq2* model. 
GSE41476_dds <- DESeq(GSE41476_dds)
save(GSE41476_dds,file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/GSE41476_dds.rda")

### CAN SKIP ################################################
## use threshold of 0.05 for adjusted Pvalue
GSE41476_results <- results(GSE41476_dds,alpha=0.05)
table(GSE41476_results$padj < 0.05)
summary(GSE41476_results,alpha=0.05)
GSE41476_results
dim(GSE41476_results)
dim(subset(GSE41476_results,padj<0.05,))
#GSE41476_results_final <- subset(GSE41476_results,padj<0.05,)
#GSE41476_results_final_sorted <- GSE41476_results_final[order(-abs(GSE41476_results_final$log2FoldChange)),]
################################################ ################################################

# use threshold of 0.05 for adjusted Pvalue AND Foldchange greater than 1.5
GSE41476_results_1.5fold <- results(GSE41476_dds,alpha=0.05,lfcThreshold=0.585)
table(GSE41476_results_1.5fold$padj < 0.05)
summary(GSE41476_results_1.5fold,alpha=0.05)
GSE41476_results_1.5fold <- subset(GSE41476_results_1.5fold,padj<0.05 & abs(log2FoldChange)>0.585)
GSE41476_results_1.5fold_sorted <- GSE41476_results_1.5fold[order(-abs(GSE41476_results_1.5fold$log2FoldChange)),]
save(GSE41476_results_1.5fold_sorted,file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/1.5fold/GSE41476_results_1.5fold_sorted.rda")
dim(GSE41476_results_1.5fold)

##################################
### Visualizing results
##################################
# MA plot
summary(GSE41476_results$log2FoldChange)
plotMA(GSE41476_results,alpha=0.05,main=c("alpha=0.05","lfcThreshold=0"),ylim=c(-21,15))

# plot dispersion estimates
par(mfrow=c(1,1))
plotDispEsts(GSE41476_dds)
#A p-value histogram:
hist(GSE41476_results$pvalue[GSE41476_results$baseMean > 1], col="grey", border="white", xlab="", ylab="", main="")
# Examine the counts for the top gene, sorting by p-value:
plotCounts(GSE41476_dds, gene=which.min(GSE41476_results$padj), intgroup="experiment_title")
GSE41476_results[which.min(GSE41476_results$padj),]
# gene with max padj value 
plotCounts(GSE41476_dds, gene=which.max(GSE41476_results$padj), intgroup="experiment_title")
GSE41476_results[which.max(GSE41476_results$padj),]

# gene with max (up) log fold change
gene_max_up <- GSE41476_results_1.5fold_sorted[which.max(GSE41476_results_1.5fold_sorted$log2FoldChange),]
plotCounts(GSE41476_dds, gene=rownames(gene_max_up), intgroup="experiment_title")
# gene with max (down) log fold change
gene_max_down <- GSE41476_results_1.5fold_sorted[which.max(abs(GSE41476_results_1.5fold_sorted$log2FoldChange)),]
plotCounts(GSE41476_dds, gene=rownames(gene_max_down), intgroup="experiment_title")


##################################
# annotations
##################################
# check & split up the rownames of the results object, which contain ENSEMBL gene ids, separated by the plus sign, +.
grep("\\+",rownames(GSE41476_results_1.5fold_sorted))
# This is the organism annotation package (“org”) for Homo sapiens (“Hs”), organized as an AnnotationDbi database package (“db”), using Entrez Gene IDs (“eg”) as primary key. 
# get a list of all available key types
columns(org.Hs.eg.db)
# use the mapIds function to add individual columns to our results table. We provide the row names of our results table as a key, and specify that keytype=ENSEMBL. The column argument tells the mapIds function which information we want, and the multiVals argument tells the function what to do if there are multiple possible values for a single input value.
GSE41476_results_1.5fold_sorted$symbol <- mapIds(org.Hs.eg.db,
                                               keys=row.names(GSE41476_results_1.5fold_sorted),
                                               column="SYMBOL",
                                               keytype="ENSEMBL",
                                               multiVals="first")
GSE41476_results_1.5fold_sorted$entrez <- mapIds(org.Hs.eg.db,
                                               keys=row.names(GSE41476_results_1.5fold_sorted),
                                               column="ENTREZID",
                                               keytype="ENSEMBL",
                                               multiVals="first")
GSE41476_results_1.5fold_sorted$genename <- mapIds(org.Hs.eg.db,
                                               keys=row.names(GSE41476_results_1.5fold_sorted),
                                               column="GENENAME",
                                               keytype="ENSEMBL",
                                               multiVals="first")


GSE41476_results_1.5fold_sorted

# Exporting results
GSE41476_results_final_sorted_DataFrame <- as.data.frame(GSE41476_results_1.5fold_sorted)
save(GSE41476_results_final_sorted_DataFrame,file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/GSE41476_results_1.5fold_sorted_DataFrame.rda")
write.csv(GSE41476_results_final_sorted_DataFrame, file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/GSE41476_results_1.5fold_sorted.csv")

