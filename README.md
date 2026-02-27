# MitoGenome_Ngsdata_analysis
Comprehensive mitochondrial genome assembly, annotation, and phylogenetic analysis of a deep-sea Dumbo Octopus (Grimpoteuthis sp.) collected from the Andaman Seamounts during  internship. The study integrates reference-based assembly using NGS data to investigate genetic variation, evolutionary relationships, and species identification.


# Project Overview
- Assemble complete mitochondrial genome
- Perform gene annotation
- Validate species identity
- Investigate evolutionary relationships

# Study Organism

*Dumbo Octopus (Grimpoteuthis sp.)  
Deep-sea cephalopod species collected from Andaman Seamounts.

# Analysis steps

Raw Illumina FASTQ  
Quality Filtering (fastp)  
Quality Assessment (FastQC)  
Reference-based Alignment (BWA)  
SAM to BAM Processing (SAMtools)  
Consensus Mitogenome Generation  
Annotation (MITOS2)  
BLAST Validation  
Phylogenetic Analysis  

# Step-by-Step Procedure

# Quality Filtering using fastp

```bash
fastp -i sample_R1.fastq -I sample_R2.fastq \
-o clean_R1.fastq -O clean_R2.fastq \
-h fastp_report.html -j fastp_report.json
```

#  Quality Check using FastQC

```bash
fastqc clean_R1.fastq clean_R2.fastq
```

# Reference-based Alignment using BWA

Index reference:

```bash
bwa index reference.fasta
```

Align reads:

```bash
bwa mem reference.fasta clean_R1.fastq clean_R2.fastq > aligned.sam
```

# SAM/BAM Processing using SAMtools

Convert SAM to BAM:

```bash
samtools view -Sb aligned.sam > aligned.bam
```

Sort BAM:

```bash
samtools sort aligned.bam -o aligned_sorted.bam
```

Index BAM:

```bash
samtools index aligned_sorted.bam
```

# Generate Consensus Mitogenome

```bash
samtools mpileup -uf reference.fasta aligned_sorted.bam | \
bcftools call -c | \
vcfutils.pl vcf2fq > consensus.fastq
```

Convert to FASTA:

```bash
seqtk seq -a consensus.fastq > mitogenome_assembly.fasta
```

# Annotation using MITOS2

- Uploaded `mitogenome_assembly.fasta` to MITOS2 web server  
- Selected invertebrate mitochondrial genetic code  
- Downloaded annotated genome and GFF files  
