#!/bin/bash

cd /Users/kmvanden/Desktop/bracken_bray_curtis

conda activate qiime2-2022.11

biom convert -i deseq2_norm_med_table.txt -o deseq2_norm_med_table.biom --to-hdf5

qiime tools import --input-path deseq2_norm_med_table.biom --type 'FeatureTable[Frequency]' --output-path deseq2_norm_med_table.qza

qiime diversity-lib bray-curtis --i-table deseq2_norm_med_table.qza --o-distance-matrix bray_curtis_deseq2_norm_med_table.qza

qiime diversity pcoa --i-distance-matrix bray_curtis_deseq2_norm_med_table.qza --o-pcoa bray_curtis_deseq2_norm_med_table.pcoa.qza

qiime emperor plot --i-pcoa bray_curtis_deseq2_norm_med_table.pcoa.qza --m-metadata-file metadata.txt --o-visualization bray_curtis_deseq2_norm_med_table.qzv

qiime diversity beta-group-significance --i-distance-matrix bray_curtis_deseq2_norm_med_table.qza --m-metadata-file metadata.txt --m-metadata-column outcome --p-method permanova --p-pairwise --p-permutations 999 --o-visualization bray_curtis_deseq2_norm_med_table_permanova.qzv

