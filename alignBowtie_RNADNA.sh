#!/bin/bash
#
#SBATCH --job-name=alignBowtie_RNADNA
#SBATCH --output=job_errors/alignBowtie_RNADNA.%j.out
#SBATCH --error=job_errors/alignBowtie_RNADNA.%j.err
#SBATCH --time=3:59:00
#SBATCH -p dpetrov,hns,owners,normal
#SBATCH --mem=10000
#SBATCH -c 16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cweiss19@stanford.edu

threads=16

#cat Oligos_libraryA.fasta Oligos_libraryB.fasta > Oligos_library_joint.fasta
#cat David-assoc1_S1_R1_001.fastq.gz David-assoc2_S2_R1_001.fastq.gz David-assoc3_S3_R1_001.fastq.gz > David-assoc_R1_001.fastq.gz
#cat David-assoc1_S1_R2_001.fastq.gz David-assoc2_S2_R2_001.fastq.gz David-assoc3_S3_R2_001.fastq.gz > David-assoc_R2_001.fastq.gz
#cat David-assoc1_S1_R3_001.fastq.gz David-assoc2_S2_R3_001.fastq.gz David-assoc3_S3_R3_001.fastq.gz > David-assoc_R3_001.fastq.gz

#fastafile = 'oligo2bc.fasta'
#fastafile = 'Oligos_library_joint_noDups.fasta'
#fastq_R1 = 'David-assoc2_S2_R1_001.fastq.gz'
#fastq_R3 = 'David-assoc2_S2_R3_001.fastq.gz'

cell='ES'
rep='3'


bowtie2-build -f /scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/alignments/oligos_to_barcodes_comb_30_3.fasta Aligned_RNA_${cell}_rep${rep}_comb
time bowtie2 \
    -x Aligned_RNA_${cell}_rep${rep}_comb \
    -1 /scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/comb_fastq/Neand_${cell}_RNA_rep${rep}_R1_comb.fastq.gz \
    -2 /scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/comb_fastq/Neand_${cell}_RNA_rep${rep}_R3_comb.fastq.gz \
    --very-sensitive \
    --threads ${threads} \
    | samtools view \
        -@ ${threads} \
        -b \
        - \
    > Aligned_RNA_${cell}_rep${rep}_comb.bam


bowtie2-build -f /scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/alignments/oligos_to_barcodes_comb_30_3.fasta Aligned_DNA_${cell}_rep${rep}_comb
time bowtie2 \    -x Aligned_DNA_${cell}_rep${rep}_comb \
    -1 /scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/comb_fastq/Neand_${cell}_gDNA_rep${rep}_R1_comb.fastq.gz \
    -2 /scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/comb_fastq/Neand_${cell}_gDNA_rep${rep}_R3_comb.fastq.gz \
    --very-sensitive \
    --threads ${threads} \
    | samtools view \
        -@ ${threads} \
        -b \
        - \
    > Aligned_DNA_${cell}_rep${rep}_comb.bam
