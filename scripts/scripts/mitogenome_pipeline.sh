#!/bin/bash

# Quality filtering
fastp -i sample_R1.fastq -I sample_R2.fastq \
-o clean_R1.fastq -O clean_R2.fastq

# Alignment
bwa index reference.fasta
bwa mem reference.fasta clean_R1.fastq clean_R2.fastq > aligned.sam

# SAM to BAM
samtools view -Sb aligned.sam > aligned.bam
samtools sort aligned.bam -o aligned_sorted.bam
samtools index aligned_sorted.bam
