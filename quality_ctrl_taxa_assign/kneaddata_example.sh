#!/bin/bash
#SBATCH --job-name kneaddata_S70 
#SBATCH --partition general
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=04:00:00 
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kmvanden@iu.edu

cd /N/scratch/kmvanden

gzip -d kneaddata_input/KMV_460_S70_R2_001.fastq.gz

sed 's/ 2/ 1/g' kneaddata_input/KMV_460_S70_R2_001.fastq > kneaddata_input/KMV_460_m_S70_R2_001.fastq

gzip -d kneaddata_input/KMV_460_S70_R1_001.fastq.gz

source activate kneaddata

kneaddata --input1 kneaddata_input/KMV_460_S70_R1_001.fastq --input2 kneaddata_input/KMV_460_m_S70_R2_001.fastq --reference-db human_genome_db --output kneaddata_460_output --trimmomatic /N/u/kmvanden/Carbonate/miniconda3/bin/Trimmomatic-0.33 --run-trim-repetitive --fastqc /N/scratch/kmvanden/fastqc/FastQC/fastqc --remove-intermediate-output -t 20

gzip kneaddata_460_output/KMV_460_S70_R1_001_kneaddata_paired_1.fastq

gzip kneaddata_460_output/KMV_460_S70_R1_001_kneaddata_paired_2.fastq


