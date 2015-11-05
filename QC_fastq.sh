# 1st Step in RNA-Seq Analysis Pipleline: Quality control for fastq files

#############################
# Connect to AWS
#############################
# cd to the directory where AWS PEM file is stored
cd ~/Documents/AWS

# configure aws credentials
aws configure # provde the AWS Access Key, Secret key & region.

# ec2 instance configuration - r3.2xlarge: 8 cores, 61 GB ram, 160 GB storage
# launch EC2 instance through GUI/command line. click "Connect" and copy ssh command
aws ec2 run-instances --image-id ami-d5ea86b5 --count 1 --instance-type r3.xlarge --key-name shruti

# ssh into the ec2 instance.
ssh -i "shruti.pem" ec2-user@<IP_ADDRESS_OF_EC2_HOST>



#############################
# Install software 
#############################
# FastQC: make sure you change the version number to the latest one.
curl -O http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/fastqc_v0.11.3.zip
unzip fastqc_v0.11.3.zip

# FASTX-Toolkit
curl -O http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
tar -xvjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2


#############################
# Download fastq files from ENA http://www.ebi.ac.uk/ena
#############################
# create a folder in ec2 where you want to store files 
mkdir fastq_files
cd fastq_files
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR58557[0-4]/SRR58557[0-4]_[1-2].fastq.gz
# cross check size of files and ensure you did not get error "Recv failure: Operation timed out"
ls -lgh
# copy the fastq files from s3 to ec2
# aws s3 ls s3://gse41476/
# mkdir fastq_files
# aws s3 cp s3://gse41476/ ./fastq_files --recursive --exclude "*" --include "*.fastq.gz"


###########################
# Run FastQC
###########################
# create output directory to store FastQC results
mkdir /home/ec2-user/fastq_files/FastQC_output
# ./../FastQC/fastqc # to open interactive mode
/home/ec2-user/FastQC/fastqc -help #./../FastQC/fastqc -help # command line mode
/home/ec2-user/FastQC/fastqc --threads 4 --extract -o FastQC_output *.fastq.gz
# copy the data from ec2 to s3 bucket
aws s3 cp /home/ec2-user/fastq_files/FastQC_output/ s3://gse41476/FastQC_output/ --recursive

# analyze the FastQC results using the python script - fastqc_validator.py. This will tell which tests failed/warning. You need to manually inspect them and use biological knowledge and decide further action
aws s3 cp s3://rna-seq-qc/FastQC_result_analyzer.py /home/ec2-user
python /home/ec2-user/FastQC_result_analyzer.py -d /home/ec2-user/fastq_files/FastQC_output -s FAIL


###########################
# Run FastX-tool-kit to remove bases with low sequence quality. 
###########################
mkdir TrimmedFastq
# unzip files that you need to run
gunzip SRR585574_*.fastq.gz

# since it is Illumina 1.5, Qulaity score Encoding is Phred+64. https://en.wikipedia.org/wiki/FASTQ_format#Encoding. [-t N] = Trim N nucleotides from the end of the read.
fastx_trimmer -Q64 -t 3 -i SRR585574_1.fastq -o ./TrimmedFastq/SRR585574_1_trimmed.fastq 
fastx_trimmer -Q64 -t 3 -i SRR585574_2.fastq -o ./TrimmedFastq/SRR585574_2_trimmed.fastq 
# run fastqc again with the trimmed files and see if all tests are pass and ok
gzip ./TrimmedFastq/*.trimmed.fastq 

# add all processed fastq files in one folder and as compressed files
mkdir ProcessedFastqFiles
mv ./TrimmedFastq/*_trimmed.fastq.gz ./ProcessedFastqFiles/
mv *.fastq.gz ./ProcessedFastqFiles/

# save files to S3
aws s3 cp /home/ec2-user/fastq_files/ProcessedFastqFiles s3://gse41476/ProcessedFastqFiles/ --recursive

# Next (2nd) Step in RNA-Seq Analysis Pipleline: Read alignment




