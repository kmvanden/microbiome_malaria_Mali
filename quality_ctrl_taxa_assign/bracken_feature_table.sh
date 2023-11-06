#!/bin/bash

cd /N/scratch/kmvanden/bracken_mpa

conda activate kraken2

combine_mpa.py --input KMV_*_S*_bracken_mpa.txt --output merged_table.txt

grep -E "s__|Classification" merged_table.txt > feature_table.txt
