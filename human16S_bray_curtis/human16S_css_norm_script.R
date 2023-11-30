# CSS normalization of 16S human data
# metagenomeSeq


# set working directory
setwd("/N/scratch/kmvanden/human16S")

# load libraries
library("metagenomeSeq") 


# read in feature table and metadata
# check that sample order is the same
otu <- read.table("otu_table_may_ten.txt")
metadata <- read.table("metadata_may_ten.txt")
sum(rownames(metadata)==colnames(otu))

# change metadata into an annotated data frame
metadata = AnnotatedDataFrame(metadata)

# create an MRexperiment object
metaSeq = newMRexperiment(otu, phenoData = metadata)

# CSS normalization of the data
p = cumNormStatFast(metaSeq) 
metaSeq = cumNorm(metaSeq, p = p)

# export normalized count matrices to a csv file (normalized counts)
otu_table_may_ten_css = MRcounts(metaSeq, norm = TRUE, log = TRUE)
otu_table_may_ten_css = as.data.frame(otu_table_may_ten_css)
write.csv(otu_table_may_ten_css, row.names= TRUE, file = "otu_table_may_ten_css.csv")


### Session Info ###
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux
# 
# Matrix products: default
# BLAS/LAPACK: /geode2/soft/hps/rhel7/intel/19.5/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin/libmkl_rt.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] metagenomeSeq_1.36.0 RColorBrewer_1.1-2   glmnet_4.1-2         Matrix_1.5-4         limma_3.50.3        
# [6] Biobase_2.58.0       BiocGenerics_0.38.0 
# 
# loaded via a namespace (and not attached):
#   [1] pillar_1.9.0       compiler_4.1.1     bitops_1.0-7       iterators_1.0.13   tools_4.1.1       
# [6] lifecycle_1.0.3    tibble_3.2.1       gtable_0.3.0       lattice_0.20-44    pkgconfig_2.0.3   
# [11] rlang_1.1.1        foreach_1.5.1      cli_3.6.1          rstudioapi_0.15.0  dplyr_1.1.2       
# [16] caTools_1.18.2     gtools_3.9.2       generics_0.1.0     vctrs_0.6.3        locfit_1.5-9.6    
# [21] grid_4.1.1         tidyselect_1.2.0   glue_1.6.2         R6_2.5.1           Wrench_1.12.0     
# [26] fansi_0.5.0        survival_3.2-13    ggplot2_3.4.2      magrittr_2.0.3     gplots_3.1.1      
# [31] scales_1.2.1       codetools_0.2-18   matrixStats_0.60.1 splines_4.1.1      shape_1.4.6       
# [36] colorspace_2.0-2   KernSmooth_2.23-20 utf8_1.2.2         munsell_0.5.0     

