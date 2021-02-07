#!/bin/bash
#
#SBATCH --job-name=check.sh
#SBATCH --output=job_errors/check.%j.out
#SBATCH --error=job_errors/check.%j.err
#SBATCH --time=1:59:00
#SBATCH -p owners,hns,dpetrov,normal
#SBATCH --mem=30000
#SBATCH -c 16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cweiss19@stanford.edu

threads=16

cell='NPC'
RNA_or_DNA='RNA'
rep='3'
S_num='18'

path1='/scratch/users/cweiss19/MPRA/DATA/initial_RNA_DNA_sequencing/fastq'
path2='/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/fastq1'
path3='/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/fastq2'
path4='/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/fastq3'
path5='/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/fastq4'
path6='/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/fastq5'
path7='/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/comb_fastq'

echo ${cell} ${RNA_or_DNA} ${rep}

cat ${path1}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz ${path2}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz ${path3}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz ${path4}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz ${path5}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz ${path6}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz > ${path7}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_R1_comb.fastq.gz 
cat ${path1}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz ${path2}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz ${path3}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz ${path4}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz ${path5}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz ${path6}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz > ${path7}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_R2_comb.fastq.gz
cat ${path1}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz ${path2}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz ${path3}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz ${path4}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz ${path5}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz ${path6}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz > ${path7}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_R3_comb.fastq.gz

echo fastq_initial
zcat ${path1}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz | wc -l
zcat ${path1}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz | wc -l
zcat ${path1}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz | wc -l

echo fastq1
zcat ${path2}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz | wc -l
zcat ${path2}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz | wc -l
zcat ${path2}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz | wc -l

echo fastq2
zcat ${path3}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz | wc -l
zcat ${path3}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz | wc -l
zcat ${path3}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz | wc -l

echo fastq3
zcat ${path4}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz | wc -l
zcat ${path4}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz | wc -l
zcat ${path4}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz | wc -l

echo fastq4
zcat ${path5}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz | wc -l
zcat ${path5}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz | wc -l
zcat ${path5}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz | wc -l

echo fastq5
zcat ${path6}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R1_001.fastq.gz | wc -l
zcat ${path6}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R2_001.fastq.gz | wc -l
zcat ${path6}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_S${S_num}_R3_001.fastq.gz | wc -l

echo combined fastq
zcat ${path7}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_R1_comb.fastq.gz | wc -l
zcat ${path7}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_R2_comb.fastq.gz | wc -l
zcat ${path7}/Neand_${cell}_${RNA_or_DNA}_rep${rep}_R3_comb.fastq.gz | wc -l






