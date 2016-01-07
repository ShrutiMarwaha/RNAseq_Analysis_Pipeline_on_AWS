# Commands to perform Read Alignment Using STAR (RNA-seq aligner) on AWS (Amazon Web Services)
 
#############################
# Connect to AWS
#############################
# cd to the directory where AWS pen file is stored
cd ~/Documents/AWS
# eb2 machine configuration - r3.2xlarge: 8 cores, 61 GB ram, 160 GB storage
# launch EC2 instance through gui and connect by terminal. click "Connect" and copy ssh command
ssh -i "shruti.pem" ec2-user@54.183.33.24
aws configure # provde the AWS Access Key and Secret key


#############################
# Install software - STAR 
#############################
wget https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz
tar -xzf STAR_2.4.2a.tar.gz
cd STAR-STAR_2.4.2a
ls /home/ec2-user/STAR-STAR_2.4.2a/bin/Linux_x86_64_static
export PATH=/home/ec2-user/STAR-STAR_2.4.2a/bin/Linux_x86_64_static:$PATH
echo $PATH
cd ..


#############################
# Download human reference-genome & gtf files
#############################
# create a folder "data" in ec2 where you want to store files 
mkdir -p data/release82/reference-genome
cd data/release82/reference-genome 
# DOWNLOAD human reference-genome from ensemble - http://uswest.ensembl.org/info/data/ftp/index.html
curl -O ftp://ftp.ensembl.org/pub/release-82/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ## reference genome
# DOWNLOAD human gtf file that contains annotation
cd data
mkdir gtf
cd gtf
curl -O ftp://ftp.ensembl.org/pub/release-82/gtf/homo_sapiens/Homo_sapiens.GRCh38.82.gtf.gz ## contains annotation

# create a bucket in S3 to store reference-genome and gtf files. You can create sub-folders within the bucket using gui or use boto
aws s3 mb s3://human-genome
aws s3 cp /home/ec2-user/data/release82 s3://human-genome --recursive


#############################
# Generate Genome Index
#############################
# if you have already saved genome index, skip the above and following step and load it from S3 using following 2 commands. Jump to read Mapping step.
# aws s3 ls s3://human-genome/release82/genome-index/
# aws s3 cp s3://human-genome/release82/genome-index/ /home/ec2-user/data/genome-index/ --recursive

# copy gtf and reference-genome files from S3 if you have already stored them.
# aws s3 cp s3://human-genome/release82 /home/ec2-user/data/release82 --recursive
cd data/release82
gunzip ./reference-genome/*.fa.gz
gunzip ./gtf/*.gtf.gz
# create output directory where genome indices will be stored
mkdir genome-index
--runMode genomeGenerate \ # directs STAR to run genome indices generation job
--runThreadN \ # how many threads to be used. depends on number of cores available.
--genomeDir \ # path to the direcotry where genome indices are stored
--genomeFastaFiles \ # path of directory containing one or more FASTA files with the genome reference sequences.
--sjdbGTFfile \ # path to the file with annotated transcripts in the standard GTF format
time STAR --runMode genomeGenerate --runThreadN 8 --genomeDir genome-index --genomeFastaFiles reference-genome/*.fa  --sjdbGTFfile gtf/Homo_sapiens.GRCh38.82.gtf --sjdbOverhang 100
# copy files from ec2 to s3, cd to the parent directory of the folder (genome-index) which you want to copy
aws s3 sync genome-index s3://human-genome/release82/genome-index/


#############################
# Read Mapping
#############################
mkdir -p data/ReadMappingInputFiles
cd ReadMappingInputFiles/
# copy the processed fastq files (that have passed QC) from s3 to ec2. Files which did not need processing can also be copied directly. I have not copied them again together in a new folder to reduce S3 storage cost.
aws s3 ls s3://gse41476/ProcessedFastqFiles/
aws s3 ls s3://gse41476/RawFastqFiles/
aws s3 cp s3://gse41476/ProcessedFastqFiles/ /home/ec2-user/data/ReadMappingInputFiles/ --recursive --exclude "*" --include "*.fastq.gz"
aws s3 cp s3://gse41476/RawFastqFiles/ /home/ec2-user/data/ReadMappingInputFiles/ --recursive --exclude "SRR585574_*.fastq.gz"
#mkdir fastq
#aws s3 cp s3://gse41476/ProcessedFastqFiles/ ./fastq --recursive --exclude "*" --include "*.fastq.gz"
#aws s3 cp s3://gse41476/ ./fastq_files --recursive --exclude "*" --include "*.fastq.gz"
gunzip *.fastq.gz


# create a folder for storing read alignment for each sample. If there are too many samples, run a loop. You can run all files togther too. Separate paired end read by space and reads from two different samples by a comma. 
mkdir -p data/AlignedReads/
cd AlignedReads/
mkdir SRR585570 SRR585571 SRR585572 SRR585573 SRR585574

cd SRR585570
--readFilesIn \ # path of directory containing sequences to be mapped (RNA-seq Fastq files). 
--outSAMtype BAM SortedByCoordinate \ output should be sorted .bam file
screen
#time STAR --runThreadN 8 --genomeDir ../genome-index/ --readFilesIn ../fastq/SRR585570_1.fastq ../fastq/SRR585570_2.fastq --outSAMtype BAM SortedByCoordinate 
time STAR --runThreadN 8 --genomeDir /home/ec2-user/data/genome-index/ --readFilesIn /home/ec2-user/data/ReadMappingInputFiles/SRR585570_1.fastq /home/ec2-user/data/ReadMappingInputFiles/SRR585570_2.fastq --outSAMtype BAM SortedByCoordinate 

# Repeat this for each sample.
cd ../SRR585571
time STAR --runThreadN 8 --genomeDir /home/ec2-user/data/genome-index/ --readFilesIn /home/ec2-user/data/ReadMappingInputFiles/SRR585571_1.fastq /home/ec2-user/data/ReadMappingInputFiles/SRR585571_2.fastq --outSAMtype BAM SortedByCoordinate 
cd ../SRR585572
time STAR --runThreadN 8 --genomeDir /home/ec2-user/data/genome-index/ --readFilesIn /home/ec2-user/data/ReadMappingInputFiles/SRR585572_1.fastq /home/ec2-user/data/ReadMappingInputFiles/SRR585572_2.fastq --outSAMtype BAM SortedByCoordinate 
cd ../SRR585573
time STAR --runThreadN 8 --genomeDir /home/ec2-user/data/genome-index/ --readFilesIn /home/ec2-user/data/ReadMappingInputFiles/SRR585573_1.fastq /home/ec2-user/data/ReadMappingInputFiles/SRR585573_2.fastq --outSAMtype BAM SortedByCoordinate 
cd ../SRR585574
time STAR --runThreadN 8 --genomeDir /home/ec2-user/data/genome-index/ --readFilesIn /home/ec2-user/data/ReadMappingInputFiles/SRR585574_1_trimmed.fastq /home/ec2-user/data/ReadMappingInputFiles/SRR585574_2_trimmed.fastq --outSAMtype BAM SortedByCoordinate 

# now copy read alignment output for each sample to S3. repeat this for each sample.
aws s3 cp SRR585570 s3://gse41476/AlignedReads/SRR585570

aws s3 cp AlignedReads/ s3://gse41476/AlignedReads/ --recursive
