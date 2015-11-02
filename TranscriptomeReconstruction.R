# Transcriptome Reconstruction: i.e., maps the aligned reads to genes. it uses .bam files (output of read alignment)

# if a new package from bioconductor has to be installed. write following two commands
# source("https://bioconductor.org/biocLite.R")
# biocLite("package name")

##################################
## load files
##################################
# create a sample-table which contains SRR, GSM ids and whether it is control or treatment sample
SRP016059_sample_table <- read.table("/Users/shruti/Dropbox/SHRUTIM/Rscripts/SRP016059/SRP016059.txt",sep="\t",header=T,stringsAsFactors=T)
SRP016059_sample_table <- SRP016059_sample_table[1:5,c("run_accession","experiment_title")]
pattern_matched <- regexpr("normal|tumor",SRP016059_sample_table$experiment_title)
condition <- regmatches(SRP016059_sample_table$experiment_title,pattern_matched)
SRP016059_sample_table[,2] <- condition
# load .bam files
SRP016059_bam_files <- file.path("/Users/shruti/Dropbox/SHRUTIM/Rscripts/SRP016059/ReadAlignments/",paste0(SRP016059_sample_table$run_accession, "Aligned.sortedByCoord.out.bam"))
# load gtf file
gtf_file <- "/Users/shruti/Dropbox/SHRUTIM/Rscripts/RNAseq/HumanGenome/Release82/gtf/Homo_sapiens.GRCh38.82.gtf"

##################################
## Transcriptome Reconstruction
##################################
# create an *Rsamtools* variable which wraps our BAM files, and create a transcript database from the GTF file. ignore the warning about `matchCircularity`. 
library(Rsamtools)
# BamFileList() provides a convenient way of managing a list of BamFile instances.
bam.list <- BamFileList(SRP016059_bam_files) # use yield size argument if it is too slow

# makeTxDbFromGFF makes a TxDb object from transcript annotations available as a GFF3 or GTF file.
library(GenomicFeatures)
# Homo_sapiens.GRCh38.82.gtf_txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
# saveDb(Homo_sapiens.GRCh38.82.gtf_txdb,"/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/Homo_sapiens.GRCh38.82.gtf_txdb.sqlite")
Homo_sapiens.GRCh38.82.gtf_txdb <- loadDb("/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/Homo_sapiens.GRCh38.82.gtf_txdb.sqlite",packageName="GenomicFeatures")
# Finally, we make a *GRangesList* which contains the exons for each gene.
exons.by.gene <- exonsBy(Homo_sapiens.GRCh38.82.gtf_txdb, by="gene")
length(exons.by.gene)
exons.by.gene[[1]]
summary(elementLengths(exons.by.gene))

# The following code chunk creates a *SummarizedExperiment* containing the counts for the reads in each BAM file (columns) for each gene in `exons.by.gene` (the rows). 
library(GenomicAlignments)
SRP016059_SummarizedOverlap <- summarizeOverlaps(exons.by.gene, bam.list,
                                                 mode="Union",
                                                 singleEnd=FALSE,
                                                 ignore.strand=TRUE,
                                                 fragments=TRUE)
# We add the sample_table as column data. we know the order is correct, because the `bam.list` was constructed from a column of `sample.table`.
colData(SRP016059_SummarizedOverlap) <- DataFrame(SRP016059_sample_table) # this step is IMPORTANT
#colnames(SRP016059_SummarizedOverlap) <- SRP016059_sample_table$run_accession
colnames(SRP016059_SummarizedOverlap) <- (colData(SRP016059_SummarizedOverlap))$run_accession
save(SRP016059_SummarizedOverlap,file="/Users/shruti/Dropbox/SHRUTIM/Rscripts/RnaSeq/GSE41476/SRP016059_SummarizedOverlap.rda")
