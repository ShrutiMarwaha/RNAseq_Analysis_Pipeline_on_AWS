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

# threshold of 0.05 for adjusted Pvalue
GSE41476_results <- results(GSE41476_dds,alpha=0.05)
table(GSE41476_results$padj < 0.05)
summary(GSE41476_results,alpha=0.05)
dim(GSE41476_results)
dim(subset(GSE41476_results,padj<0.05,))

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
GSE41476_results_DataFrame <- as.data.frame(GSE41476_results_1.5fold_sorted)

# TO DO: find why so many genes do not have a corresponding entrez id
dim(GSE41476_results_DataFrame)
missing <- is.na(GSE41476_results_DataFrame$entrez)
length(GSE41476_results_DataFrame$entrez[missing])

#function to remove rows which are not associated with any gene 
RemoveRows <- function(MatrixName,ColumnName){ #ColumnName is name of column containing gene id
  # remove rows that contain "NA" in any column
  complete <- complete.cases(MatrixName)
  Matrix_complete <- MatrixName[complete,]
  # remove rows that contain "" in column for gene id
  nonempty_gene_ids <- Matrix_complete[,ColumnName]!=""
  Matrix2 <- Matrix_complete[nonempty_gene_ids,]
  
  # if GeneID contains more than 1 id, remove additional ones
  if (length(grep("/",Matrix2[,ColumnName])>0))
  {
    GeneIds <- unlist(lapply(Matrix2[,ColumnName],function(x) {
      sub("[ ]/.+","",x)
    }))
    Matrix2[,ColumnName] <- GeneIds
  }   
  
  return(Matrix2)
}

DEG_gse41476 <- RemoveRows(GSE41476_results_DataFrame,"entrez")
gse41476_up <- unique(subset(DEG_gse41476,log2FoldChange>0,)$entrez)
gse41476_down <- unique(subset(DEG_gse41476,log2FoldChange<0,)$entrez)
length(gse41476_up)
length(gse41476_down)
save(DEG_gse41476,file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/DEG_gse41476.rda")

####   Gene Signature Based Drug Repurposing: Connectivity Map /LINCS analysis #########
# use biomaRt package to convert entrez ids to hgu1331a afffy probeset ids.}
library("biomaRt")
listMarts(host="www.ensembl.org")
database <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
listDatasets(database)[grep("sapiens",listDatasets(database)$description,),]

## Filters (one or more) that should be used in the query. 
filters <- listFilters(database)
filters[grep("entrez",filters$description,ignore.case=T),]

## attribites are values that you are interested in to retrieve
attributes <- listAttributes(database)
grep("hg.?u.?133.?a",attributes$description,ignore.case=T)
grep("entrez",attributes$description,ignore.case=T)
attributes[c(58,59,105,106),]


# function to convert entrez ids to hgu1331a afffy probeset ids
GenerateCmapInputList <- function(MatrixName,ColumnName){ #ColumnName is name of column containing gene id
  # selecting columns
  CmapInput <- MatrixName[,c(ColumnName,"log2FoldChange","padj")]
  CmapInput <- CmapInput[order(-abs(CmapInput$log2FoldChange)),]
  
  DEG_EntrezGeneID <- as.character(unique(CmapInput[,ColumnName]))
  DEG_probesetIDs <- getBM(attributes=c('entrezgene','affy_hg_u133a'), filters = 'entrezgene', values = DEG_EntrezGeneID, mart = database, uniqueRows=T)
  DEG_probesetIDs <- subset(DEG_probesetIDs,affy_hg_u133a!="")
  DEG_probesetIDs_FC <- merge(x=CmapInput[,c(ColumnName,"log2FoldChange")],y=DEG_probesetIDs,by.x=ColumnName,by.y="entrezgene",sort=F,all.y=F)
  DEG_probesetIDs_FC <- DEG_probesetIDs_FC[!duplicated(DEG_probesetIDs_FC$affy_hg_u133a),]
  DEG_probesetIDs_FC <- DEG_probesetIDs_FC[order(-abs(DEG_probesetIDs_FC$log2FoldChange)),]
  
  return(DEG_probesetIDs_FC)
}

# Generate Cmap Input List 
GSE41476_CmapInput <- GenerateCmapInputList(DEG_gse41476,"entrez")
GSE41476_CmapInput_Up <- subset(GSE41476_CmapInput,log2FoldChange>0,)$affy_hg_u133a
GSE41476_CmapInput_Down <- subset(GSE41476_CmapInput,log2FoldChange<0,)$affy_hg_u133a
## count no.of up and Down DEGs to get idea about iterations you will like to run 
length(GSE41476_CmapInput_Up)
length(GSE41476_CmapInput_Down)
# select top 500 up and down genes for standard cmap input
write.table(GSE41476_CmapInput_Up[1:500],file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/Cmap/GSE41476_CmapInput_Up.txt.grp",sep="\t",quote=F,col.names=F,row.names=F)
write.table(GSE41476_CmapInput_Down[1:500],file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/Cmap/GSE41476_CmapInput_Down.txt.grp",sep="\t",quote=F,col.names=F,row.names=F)

setwd("/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/Cmap/Iteration")
## save probe ids for Up DEGs
for(i in seq(from=1,by=100,to=500))
{
  cMapInputUp <- GSE41476_CmapInput_Up[1:(i+99)]
  write.table(cMapInputUp,file=sprintf("GSE41476CmapInputUpTop%i.txt.grp",(i+99)),sep="\t",quote=F,col.names=F,row.names=F)
}

## save probe ids for Down DEGs
for(i in seq(from=1,by=100,to=500))
{
  cMapInputDown <- GSE41476_CmapInput_Down[1:(i+99)]
  write.table(cMapInputDown,file=sprintf("GSE41476CmapInputDownTop%i.txt.grp",(i+99)),sep="\t",quote=F,col.names=F,row.names=F)
}
## check the last files where there can be NAs if no.of DEgs is not exactly equal to iteration interval
##########################################################

## run iterations in cmap: 100up, 100 dn; 200up, 200 dn; 300up, 300 dn; 400up, 400 dn; 500up, 500 dn;
# provide the iteration interval on which cmap results have been named
Iteration <- c(200,400,600,800,1000)
## following list will store compound names from each iteration
CompoundList <- list()

library("xlsx")
for(i in seq_along(Iteration))  
{
  ## each cmap result file has been named as "gse41476CmapResult" followed by the number of input genes
  CmapResults <- read.xlsx2(file=sprintf("gse41476CmapResultTop%i.xls",Iteration[i]),sheetIndex=1,colClasses=c(c("numeric","character"),rep("numeric",4)))
  ## select compounds which have p.value < 0.05 and negative enrichment score
  CompoundList[[i]] <- subset(CmapResults,p<0.05 & enrichment<0,cmap.name,drop=T)
}
## find common compounds between hits from each iteration
Reduce(intersect,list(CompoundList[[1]],CompoundList[[2]],CompoundList[[3]],CompoundList[[4]],CompoundList[[5]]))
CommonCompoundsGse41476 <- Reduce(intersect,list(CompoundList[[2]],CompoundList[[3]],CompoundList[[4]],CompoundList[[5]]))
save(CommonCompoundsGse41476,file="./CommonCompoundsGse41476.rda")
