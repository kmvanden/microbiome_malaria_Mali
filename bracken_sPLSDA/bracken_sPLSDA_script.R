# sparse partial least squares discriminant analysis
# mixOmics
# Kraken2/Bracken feature table reduced to the top 350 taxa based on median counts
# feature table normalized according to DESeq2 defaults - relative log expression


# set working directory
setwd("/N/scratch/kmvanden/metagen_sPLSDA")

# load libraries
library(mixOmics)


# read in feature table and metadata
# check that sample order is the same
X_short = read.table("deseq2_norm_med_table.txt", header = TRUE)
X_short = as.matrix(X_short) # convert to matrix
X_short = t(X_short) # transpose matrix
class(X_short)
dim(X_short)

Y = read.table("metadata.txt", header = TRUE)
sum(row.names(X_short)==row.names(Y))
Y = Y$outcome
Y = as.factor(Y)
class(Y)
summary(Y)


### preliminary PCA analysis - variables zero centered
pca.mali = pca(X_short, ncomp = 5, scale = TRUE)
plotIndiv(pca.mali, group = Y, ind.names = FALSE, comp = c(1,2),
          legend = TRUE, title = "PCA of metagenomics data")


### evaluate classification performance of PLS-DA - using all variables
# run perf() with 12-fold cross-validation repeated 50 times for 5 components
plsda.mali = plsda(X_short, Y, ncomp = 5)
perf.plsda.mali = perf(plsda.mali, validation = "Mfold", folds = 12, nrepeat = 50)
plot(perf.plsda.mali)

perf.plsda.mali$choice.ncomp # optimal number of components by prediction distance and error rate 
# error rate is at a minimum with one component

perf.plsda.mali$error.rate # error rate for each component


### variable selection for sparse method
# estimation of classification error rate with respect to the number of variables
# run tune.splsda() with 12-fold cross-validation repeated 50 times for 5 components 
list.keepX = c(1:10, seq(15, 350, 5))
tune.splsda.mali = tune.splsda(X_short, Y, ncomp = 5, validation = "Mfold",
                               folds = 12, nrepeat = 50, dist = "centroids.dist",
                               test.keepX = list.keepX)

### plot mean classification error rate on each component across all tested keepX values
pdf("metagen_sPLSDA_tune.pdf", width = 10, height = 10)
plot(tune.splsda.mali, sd = TRUE)
dev.off()

tune.splsda.mali$choice.ncom$ncomp # optimal number of components
# [1] 1
tune.splsda.mali$choice.keepX # optimal keepX parameters according to minimal error rate
# comp1 comp2 comp3 comp4 comp5 
# 160    85     8     6     1 


### final model and performance
# optimal number of components
ncomp = 2 # ncomp increased to 2 to allow visualization of model using mixOmics functions

# optimal number of variables 
select.keepX = tune.splsda.mali$choice.keepX[1:ncomp]

# final model
splsda.mali = splsda(X_short, Y, ncomp = ncomp, keepX = select.keepX)

# model performance
# performance of model assessed with perf() using 12-fold cross-validation repeated 50 times
perf.splsda.mali = perf(splsda.mali, folds = 12, validation = "Mfold",
                        dist = "centroids.dist", nrepeat = 50)

perf.splsda.mali$error.rate # error rate for each component
# overall error rate for component 1: 0.2968571
perf.splsda.mali$error.rate.class # error rate per class
# error rate for resistant + component 1: 0.4775000
# error rate for susceptible + component 1: 0.2026087


### variable selection and stability
# frequency at which the same taxa are selected across the folds
pdf("feature_stability_comp1.pdf", width = 5, height = 5)
stable.comp1 = perf.splsda.mali$features$stable$comp1
barplot(stable.comp1, xlab = "variables selected across CV folds",
        ylab = "stability frequency", main = "feature stability for comp 1")
dev.off()


### selected variables for each component with their frequency and loading values
# output the taxa selected for a given component and their loading values 
# concatenate those results with the feature stability 
select.name = selectVar(splsda.mali, comp = 1)$name # extract name of selected variable
stability = perf.splsda.mali$features$stable$comp1[select.name] # extract the stability values from perf()
freq_loadings = cbind(selectVar(splsda.mali, comp = 1)$value, stability) # frequency and loading values
write.csv(as.data.frame(freq_loadings), file = "sPLSDA_freq_loadings.csv")


# sPLS-DA sample plot
pdf("plot_indiv_comp1_2.pdf", width = 5, height = 5)
plotIndiv(splsda.mali, comp = c(1,2),
          ind.names = TRUE, legend = TRUE, 
          star = FALSE, ellipse = TRUE, ellipse.level = 0.95, # 95% confidence ellipses
          col = c("dodgerblue3", "firebrick3"), # color of the groups
          cex = 3, pch = c(1,1), point.lwd = 0.5, # size, shape and line width of points
          title = "sPLS-DA, comp 1 + 2")
dev.off()


# correlation circle plots
pdf("plot_corr_circle_noname.pdf", width = 10, height = 10)
plotVar(splsda.mali, comp = c(1,2), var.names = FALSE,
        style = "graphics", legend = TRUE,
        pch = 1, cex = 0.5, col = "deepskyblue3", # shape, size and color of the points
        cutoff = 0.35 # correlation coefficient cut-off
)
dev.off()

pdf("plot_corr_circle_large.pdf", width = 30, height = 30) # to see names of taxa clearly
plotVar(splsda.mali, comp = c(1,2), var.names = TRUE,
        style = "graphics", legend = TRUE,
        pch = 1, cex = 0.5, col = "deepskyblue3", # shape, size and color of the points
        cutoff = 0.35 # correlation coefficient cut-off
)
dev.off()


# plot of loading weights for each variable
pdf("plot_loadings_comp1.pdf", width = 10, height = 7)
plotLoadings(splsda.mali, 
             comp = 1, # component of interest
             contrib = "max", # color of bar corresponds to group with the maximal expression levels
             method = "median", # criterion to assess the contribution
             ndisplay = 30, # how many of the most important variables are plotted
             legend.color = c("dodgerblue3", "firebrick3"), # color of the outcomes
             legend = FALSE
)
dev.off()


# AUROC plot
# AUC calculated from training cross-validation sets and averaged
pdf("auroc_plot_comp1.pdf", width = 5, height = 5)
auc.mali = auroc(splsda.mali, roc.comp = 1)
dev.off()
# AUC for comp 1 is 0.8243
# p-value: 9.44e-06


##### session info #####
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux
# 
# Matrix products: default
# BLAS:   /geode2/soft/hps/rhel7/r/4.2.1/lib64/R/lib/libRblas.so
# LAPACK: /geode2/soft/hps/rhel7/r/4.2.1/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] mixOmics_6.20.0     ggplot2_3.4.3       lattice_0.21-8      MASS_7.3-60         BiocManager_1.30.20 phyloseq_1.40.0    
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-162                bitops_1.0-7                matrixStats_1.0.0           bit64_4.0.5                 RColorBrewer_1.1-3         
# [6] httr_1.4.6                  GenomeInfoDb_1.32.4         tools_4.2.1                 utf8_1.2.3                  R6_2.5.1                   
# [11] vegan_2.6-4                 DBI_1.1.3                   BiocGenerics_0.42.0         mgcv_1.8-42                 colorspace_2.1-0           
# [16] permute_0.9-7               rhdf5filters_1.8.0          ade4_1.7-22                 withr_2.5.0                 tidyselect_1.2.0           
# [21] gridExtra_2.3               DESeq2_1.36.0               bit_4.0.5                   compiler_4.2.1              cli_3.6.1                  
# [26] Biobase_2.56.0              DelayedArray_0.22.0         labeling_0.4.2              scales_1.2.1                genefilter_1.78.0          
# [31] stringr_1.5.0               digest_0.6.31               XVector_0.36.0              pkgconfig_2.0.3             MatrixGenerics_1.8.1       
# [36] fastmap_1.1.1               rlang_1.1.1                 rstudioapi_0.15.0           RSQLite_2.3.1               farver_2.1.1               
# [41] generics_0.1.3              jsonlite_1.8.5              BiocParallel_1.30.4         dplyr_1.1.2                 RCurl_1.98-1.12            
# [46] magrittr_2.0.3              GenomeInfoDbData_1.2.8      biomformat_1.24.0           Matrix_1.5-4.1              Rcpp_1.0.10                
# [51] munsell_0.5.0               S4Vectors_0.34.0            Rhdf5lib_1.18.2             fansi_1.0.4                 ape_5.7-1                  
# [56] lifecycle_1.0.3             stringi_1.7.12              SummarizedExperiment_1.26.1 zlibbioc_1.42.0             rhdf5_2.40.0               
# [61] plyr_1.8.8                  grid_4.2.1                  blob_1.2.4                  parallel_4.2.1              ggrepel_0.9.3              
# [66] crayon_1.5.2                Biostrings_2.64.1           splines_4.2.1               multtest_2.52.0             annotate_1.74.0            
# [71] KEGGREST_1.36.3             locfit_1.5-9.7              pillar_1.9.0                igraph_1.4.3                GenomicRanges_1.48.0       
# [76] corpcor_1.6.10              geneplotter_1.74.0          reshape2_1.4.4              codetools_0.2-19            stats4_4.2.1               
# [81] XML_3.99-0.14               glue_1.6.2                  data.table_1.14.8           png_0.1-8                   vctrs_0.6.3                
# [86] foreach_1.5.2               purrr_1.0.2                 tidyr_1.3.0                 gtable_0.3.3                cachem_1.0.8               
# [91] xtable_1.8-4                RSpectra_0.16-1             survival_3.5-5              rARPACK_0.11-0              tibble_3.2.1               
# [96] iterators_1.0.14            AnnotationDbi_1.58.0        memoise_2.0.1               IRanges_2.30.1              ellipse_0.4.5              
# [101] cluster_2.1.4  

