# RNA-seq Differential Gene Expression Analysis 

## Overview

This repository contains a bioinformatics pipeline for analyzing **RNA-seq data** to identify **differentially expressed genes (DEGs)** between two chile pepper (*Capsicum annuum*) cultivars:
- Resistant cultivar
- Susceptible cultivar

RNA-seq reads were obtained from **root samples collected 24 hours after pathogen inoculation**.  
The analysis pipeline performs genome alignment, read counting, and differential gene expression analysis.

This project was conducted between **March and April 2024** using **Linux and R** for graduate coursework.

## Pipeline Summary

The analysis workflow includes the following steps:

1. Project directory setup  
2. Reference genome download  
3. Genome annotation download  
4. RNA-seq data acquisition  
5. Decompression of FASTQ files  
6. Verification of sequencing reads  
7. Quality control with FastQC  
8. Read trimming using Trimmomatic  
9. Genome indexing with STAR  
10. RNA-seq read alignment  
11. BAM file verification  
12. Gene expression quantification using featureCounts  
13. Differential expression analysis using DESeq2  

## Software Requirements

The following tools are required to run the pipeline:

| Software | Purpose |
|--------|--------|
| FastQC | Quality control of sequencing reads |
| Trimmomatic | Adapter removal and read trimming |
| STAR | RNA-seq read alignment |
| SAMtools | BAM file processing |
| featureCounts | Gene-level read counting |
| R | Statistical analysis |
| DESeq2 | Differential gene expression analysis |


## Running the Pipeline
Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/RNAseq-deg-analysis.git
cd RNAseq-deg-analysis
```

Run the pipeline script:

```bash
bash script.sh
```

The script performs preprocessing, alignment, and read counting for RNA-seq data.

## RNA-seq Data

The dataset consists of **paired-end RNA-seq reads** representing **three biological replicates** for each cultivar.

Replicate labeling used during sequencing:

| Color | Replicate |
|-----|-----|
| Black | Replicate 1 |
| Green | Replicate 2 |
| Red | Replicate 3 |

## Quality Control

Raw sequencing reads are evaluated using **FastQC** to assess:

- Per-base sequence quality  
- Adapter contamination  
- GC content distribution  
- Sequence duplication levels  

## Differential Expression Analysis

Gene expression quantification is performed using **featureCounts**, followed by **differential expression analysis in R using DESeq2**.

The DESeq2 workflow includes:

- Filtering low-expression genes  
- Normalization of read counts  
- Statistical testing for differential expression  
- Identification of significantly upregulated genes  

Significant genes are identified using:

```
adjusted p-value < 0.01
```

## Output Files

Typical outputs generated during the analysis include:

| File | Description |
|-----|-----|
| `Aligned.sortedByCoord.out.bam` | Sorted RNA-seq alignment file |
| `*_readcountoutput.txt` | Gene-level read counts |
| `DEG_results.csv` | Differential gene expression results |
| `upregulated_genes.csv` | Significantly upregulated genes |


## Author
**Ehtisham S. Khokhar**

Graduate Bioinformatics Coursework  
RNA-seq Differential Gene Expression Analysis
