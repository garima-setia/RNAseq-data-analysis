#!/bin/bash

#download reference genome sequences (FASTA files) and annotations (GTF file)

wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz

#install star
conda install -c bioconda star

#unzip reference sequence files
gunzip *.gz
mkdir ref

#generating genome indexes files
STAR --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.107.gtf --runThreadN 16

#import fastq sequences into fastq directory
mkdir mapped
cd fastq

#aligning reads: STAR maps the reads to the genome, and writes output files
for file in *.fastq; do STAR --runMode alignReads --genomeDir ../ref/ --outSAMtype BAM SortedByCoordinate --readFilesIn ${file} --runThreadN 12 --outFileNamePrefix ../mapped/${file}; done
mkdir bams
mv mapped/*.bam bams/

conda install -c bioconda subread

#generate count table
featureCounts -a Homo_sapiens.GRCh38.107.gtf -o count.out -T 8 bams/*.bam