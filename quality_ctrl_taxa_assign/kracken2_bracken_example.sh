#!/bin/bash
#SBATCH --job-name k2b_KMV_460_S70.sh 
#SBATCH --partition general
#SBATCH -A general
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=00:20:00 
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kmvanden@iu.edu

cd /N/scratch/kmvanden

source activate kraken2

gzip -d kneaddata_output/KMV_460_S70_R1_001_kneaddata_paired_1.fastq.gz

gzip -d kneaddata_output/KMV_460_S70_R1_001_kneaddata_paired_2.fastq.gz

kraken2 --db stdk2_db --threads 20 --report-zero-counts --report-minimizer-data --minimum-hit-groups 4 --confidence 0.10 --report k2_report/KMV_460_S70.k2report --paired kneaddata_output/KMV_460_S70_R1_001_kneaddata_paired_1.fastq kneaddata_output/KMV_460_S70_R1_001_kneaddata_paired_2.fastq > k2_output/KMV_460_S70.kraken2

bracken -d stdk2_db -i k2_report/KMV_460_S70.k2report -r 100 -t 10 -l S -o bracken_output/KMV_460_S70.bracken -w bracken_report/KMV_460_S70.breport

kreport2mpa.py --report bracken_report/KMV_460_S70.breport --display-header --no-intermediate-ranks --read_count --output bracken_mpa/KMV_460_S70_bracken_mpa.txt

