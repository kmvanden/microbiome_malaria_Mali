#!/bin/bash

cd /Users/kmvanden/Desktop/human_16S

source activate qiime2-2022.11

biom convert -i otu_table_may_ten_css.txt -o otu_table_may_ten_css.biom --to-hdf5

qiime tools import --input-path otu_table_may_ten_css.biom --type 'FeatureTable[Frequency]' --output-path otu_table_may_ten_css.qza

qiime diversity-lib bray-curtis --i-table otu_table_may_ten_css.qza --o-distance-matrix bray_curtis_otu_table_may_ten_css.qza

qiime diversity pcoa --i-distance-matrix bray_curtis_otu_table_may_ten_css.qza --o-pcoa bray_curtis_otu_table_may_ten_css.pcoa.qza

qiime emperor plot --i-pcoa bray_curtis_otu_table_may_ten_css.pcoa.qza --m-metadata-file metadata_may_ten.txt --o-visualization bray_curtis_otu_table_may_ten_css.qzv

qiime diversity beta-group-significance --i-distance-matrix bray_curtis_otu_table_may_ten_css.qza --m-metadata-file metadata_may_ten.txt --m-metadata-column outcome1 --p-method permanova --p-pairwise --p-permutations 999 --o-visualization bray_curtis_otu_table_may_ten_css_outcome1_permanova.qzv

qiime diversity beta-group-significance --i-distance-matrix bray_curtis_otu_table_may_ten_css.qza --m-metadata-file metadata_may_ten.txt --m-metadata-column outcome2 --p-method permanova --p-pairwise --p-permutations 999 --o-visualization bray_curtis_otu_table_may_ten_css_outcome2_permanova.qzv
