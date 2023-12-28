# preprocessing of metabolomics data
# metaboanalystR


# set working directory
setwd("/N/scratch/kmvanden/metabo_sPLSDA")

# load libraries
library(MetaboAnalystR)


# construct data object and read in data
# one feature with a single value across samples was removed from analysis
mSet1 = InitDataObjects(data.type = "pktable", anal.type = "stat", paired = FALSE)
mSet1 = Read.TextData(mSetObj = mSet1, filePath = "CHEMID_peak_int.txt", 
                      format = "colu", lbl.type = "disc")
mSet1 = SanityCheckData(mSet1)


# missing/zero values replace by 1/5 the smallest positive value for each variable
mSet1 = ReplaceMin(mSet1)
mSet1 = SanityCheckData(mSet1)

# 40% of features filtered by IQR
mSet1 = FilterVariable(mSetObj = mSet1, 
                       filter = "iqr", filter.cutoff = 40, 
                       qcFilter = "F", privileged = F)

# normalization by median and log10 transformation
mSet1 = PreparePrenormData(mSet1)
mSet1 = Normalization(mSetObj = mSet1, 
                      rowNorm = "MedianNorm", transNorm = "LogNorm",
                      scaleNorm = "NULL", ratio = FALSE, ratioNum = 20)

### exporting data
mSet_filt_norm = mSet1$dataSet$norm
mSet_filt_norm = t(mSet_filt_norm)
write.csv(as.data.frame(mSet_filt_norm_trans), file = "metabo_filt_norm.csv")


### SESSION INFO ###
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux
# 
# Matrix products: default
# BLAS:   /geode2/soft/hps/rhel7/r/4.2.1/lib64/R/lib/libRblas.so
# LAPACK: /geode2/soft/hps/rhel7/r/4.2.1/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] httr_1.4.6           ggplot2_3.4.3        memoise_2.0.1        Rserve_1.8-11        pheatmap_1.0.12      MetaboAnalystR_4.0.0
# 
# loaded via a namespace (and not attached):
#   [1] fgsea_1.22.0         colorspace_2.1-0     class_7.3-22         siggenes_1.70.0      rstudioapi_0.15.0    listenv_0.9.0       
# [7] farver_2.1.1         ggrepel_0.9.3        bit64_4.0.5          prodlim_2023.03.31   fansi_1.0.4          lubridate_1.9.2     
# [13] codetools_0.2-19     splines_4.2.1        cachem_1.0.8         impute_1.70.0        scrime_1.3.5         glasso_1.11         
# [19] jsonlite_1.8.5       pROC_1.18.2          Cairo_1.6-0          caret_6.0-94         compiler_4.2.1       Matrix_1.5-4.1      
# [25] fastmap_1.1.1        lazyeval_0.2.2       limma_3.52.4         cli_3.6.1            crmn_0.0.21          htmltools_0.5.5     
# [31] tools_4.2.1          igraph_1.4.3         gtable_0.3.3         glue_1.6.2           reshape2_1.4.4       dplyr_1.1.2         
# [37] fastmatch_1.1-3      Rcpp_1.0.10          Biobase_2.56.0       vctrs_0.6.3          multtest_2.52.0      nlme_3.1-162        
# [43] iterators_1.0.14     timeDate_4022.108    gower_1.0.1          stringr_1.5.0        globals_0.16.2       timechange_0.2.0    
# [49] lifecycle_1.0.3      gtools_3.9.4         future_1.32.0        edgeR_3.38.4         MASS_7.3-60          scales_1.2.1        
# [55] ipred_0.9-14         pcaMethods_1.88.0    parallel_4.2.1       RColorBrewer_1.1-3   qs_0.25.5            curl_5.0.1          
# [61] gridExtra_2.3        rpart_4.1.19         stringi_1.7.12       RSQLite_2.3.1        foreach_1.5.2        caTools_1.18.2      
# [67] BiocGenerics_0.42.0  hardhat_1.3.0        BiocParallel_1.30.4  lava_1.7.2.1         rlang_1.1.1          pkgconfig_2.0.3     
# [73] bitops_1.0-7         lattice_0.21-8       purrr_1.0.2          recipes_1.0.8        htmlwidgets_1.6.2    labeling_0.4.2      
# [79] bit_4.0.5            tidyselect_1.2.0     parallelly_1.36.0    plyr_1.8.8           magrittr_2.0.3       R6_2.5.1            
# [85] gplots_3.1.3         generics_0.1.3       DBI_1.1.3            pillar_1.9.0         withr_2.5.0          survival_3.5-5      
# [91] nnet_7.3-19          tibble_3.2.1         future.apply_1.11.0  KernSmooth_2.23-21   utf8_1.2.3           ellipse_0.4.5       
# [97] RApiSerialize_0.1.2  plotly_4.10.2        locfit_1.5-9.7       grid_4.2.1           data.table_1.14.8    blob_1.2.4          
# [103] ModelMetrics_1.2.2.2 digest_0.6.31        tidyr_1.3.0          RcppParallel_5.1.7   stats4_4.2.1         munsell_0.5.0       
# [109] stringfish_0.15.8    viridisLite_0.4.2 




