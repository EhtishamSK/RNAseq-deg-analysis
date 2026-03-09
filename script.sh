# RNA-seq Differential Expression Analysis Pipeline
# Author: Ehtisham Khokhar
# This project explores gene expression variation between two chile pepper cultivars
# (resistant and susceptible) using RNA-seq data collected from root samples 24 hours after inoculation.
# The analysis was performed between March and April 2024 using Linux and R.

# STEP 1: CREATE PROJECT DIRECTORIES

# Create working directories to organize the analysis
mkdir final_project
mkdir pepper_genome
mkdir resistant_72
mkdir susceptible_72
# Directory descriptions
# pepper_genome   -> reference genome files
# resistant_72    -> RNA-seq reads for resistant cultivar
# susceptible_72  -> RNA-seq reads for susceptible cultivar

# STEP 2: DOWNLOAD PEPPER REFERENCE GENOME
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/capsicum_annuum/dna/Capsicum_annuum.ASM51225v2.dna.toplevel.fa.gz

# STEP 3: DOWNLOAD GENOME ANNOTATION
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/plants/gtf/capsicum_annuum/Capsicum_annuum.ASM51225v2.57.chr.gtf.gz

# STEP 4: OBTAIN RNA-seq DATA
# RNA-seq data were provided through a company download portal.
# Because wget was not supported, the files were downloaded manually
# on Windows and transferred to the Linux server using FileZilla.
# Replicate information:
# Black = Replication 1
# Green = Replication 2
# Red = Replication 3

# STEP 5: UNZIP FASTQ FILES
# Decompress all sequencing files
gzip -d *.gz

# STEP 6: VERIFY FORWARD AND REVERSE READ FILES
# Checking line counts ensures forward and reverse files match

wc -l 24-root-resistent-black_S23_L001_R1_001.fastq
wc -l 24-root-resistent-black_S23_L001_R2_001.fastq
wc -l 24-root-resistent-green_S22_L001_R1_001.fastq
wc -l 24-root-resistent-green_S22_L001_R2_001.fastq
wc -l 24-root-resistent-red_S21_L001_R1_001.fastq
wc -l 24-root-resistent-red_S21_L001_R2_001.fastq

wc -l 24-root-susceptible-black_S8_L001_R1_001.fastq
wc -l 24-root-susceptible-black_S8_L001_R2_001.fastq
wc -l 24-root-susceptible-green_S7_L001_R1_001.fastq
wc -l 24-root-susceptible-green_S7_L001_R2_001.fastq
wc -l 24-root-susceptible-red_S6_L001_R1_001.fastq
wc -l 24-root-susceptible-red_S6_L001_R2_001.fastq

# Example outputs (line counts should match between pairs)
# 28156508 24-root-resistent-black_S23_L001_R1_001.fastq
# 28156508 24-root-resistent-black_S23_L001_R2_001.fastq
# 34092008 24-root-resistent-green_S22_L001_R1_001.fastq
# 34092008 24-root-resistent-green_S22_L001_R2_001.fastq

# STEP 7: QUALITY CONTROL WITH FASTQC
mkdir FastQC_output
fastqc 24-root-resistent-black_S23_L001_R1_001.fastq -o FastQC_output/
fastqc 24-root-resistent-black_S23_L001_R2_001.fastq -o FastQC_output/
# Repeat for all replicates (black, green, red) and susceptible samples

# STEP 8: TRIM RNA-seq READS USING TRIMMOMATIC
# trimming command for resistant replicate 1 (black)

trimmomatic PE -phred33 \
24-root-resistent-black_S23_L001_R1_001.fastq \
24-root-resistent-black_S23_L001_R2_001.fastq \
24-root-resistent-black_S23_L001_R1_001_paired.fq \
24-root-resistent-black_S23_L001_R1_001_unpaired.fq \
24-root-resistent-black_S23_L001_R2_001_paired.fq \
24-root-resistent-black_S23_L001_R2_001_unpaired.fq \
ILLUMINACLIP:/home/cdb3ny/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:80
# Repeat for all replicates (black, green, red) and susceptible samples

# Example output
# Input Read Pairs: 7039127
# Both Surviving: 6541429 (92.93%)
# Forward Only Surviving: 275086 (3.91%)
# Reverse Only Surviving: 93730 (1.33%)
# Dropped: 128882 (1.83%)

# STEP 9: GENERATE STAR GENOME INDEX
# Pepper genome is large and requires additional RAM
nohup STAR \
--runMode genomeGenerate \
--runThreadN 1 \
--genomeDir . \
--genomeFastaFiles Capsicum_annuum.ASM51225v2.dna.toplevel.fa \
--limitGenomeGenerateRAM 32958142048

# Example STAR output
# Nov 11 08:27:50 started STAR run
# Nov 11 11:22:23 finished successfully

## STEP 10: ALIGN RNA-seq READS USING STAR
# Attempted running multiple replicates together but encountered segmentation fault
# Likely due to memory limitations
# Final strategy: align each replicate individually
STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir pepper_genome \
--readFilesIn resistant24/24-root-resistent-black_S23_L001_R1_001_paired.fq \
resistant24/24-root-resistent-black_S23_L001_R2_001_paired.fq \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix resistant_black24

# Example STAR alignment output
# started mapping
# finished mapping
# started sorting BAM
# finished successfully
# Repeat alignment for green and red replicates
# Repeat entire process for susceptible cultivar samples

# STEP 11: VERIFY BAM FILE FORMAT
# Inspect BAM file
samtools view output_black24.bam | less
# Confirm BAM header
samtools view -H resistant_black24Aligned.sortedByCoord.out.bam

## STEP 12: COUNT READS USING FEATURECOUNTS
featureCounts \
-a pepper_genome/Capsicum_annuum.ASM51225v2.57.chr.gtf \
-o black24_resistant_readcountoutput.txt \
-p -B -C -Q 10 -T 8 \
resistant_black24Aligned.sortedByCoord.out.bam

# Example output summary
# Total alignments : 6540466
# Successfully assigned alignments : 4180520 (63.9%)

featureCounts \
-a pepper_genome/Capsicum_annuum.ASM51225v2.57.chr.gtf \
-o black24_susceptible_readcountoutput.txt \
-p -B -C -Q 10 -T 8 \
susceptible_black24Aligned.sortedByCoord.out.bam

# Example output summary
# Total alignments : 7128574
# Successfully assigned alignments : 4668524 (65.5%)


## STEP 13: DIFFERENTIAL EXPRESSION ANALYSIS WITH DESeq2 (R)
# The following steps were performed in R

Install required packages
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")

Load libraries
library(DESeq2)
library(apeglm)
library(dplyr)
library(ggplot2)

# Load count matrix
count_data <- read.csv("all_counts.csv", header=TRUE, row.names=1)
col_data <- read.csv("info.csv", header=TRUE)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=count_data,
                               colData=col_data,
                               design=~condition)

# Remove low-expression genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run differential expression
dds <- DESeq(dds)
res <- results(dds)

# Export results
write.csv(res, "DEG_results.csv")

# Extract significantly upregulated genes
upregulated <- subset(res, padj < 0.01 & log2FoldChange > 0)
write.csv(upregulated, "upregulated_genes.csv")

# MA plot
plotMA(res)

# Volcano-style visualization
ggplot(as.data.frame(res), aes(x=log10(baseMean), y=log2FoldChange)) +
geom_point()
