# multi-block sparse partial least squares discriminant analysis
# data integration analysis for biomarker discovery using latent components 
# mixOmics - DIABLO

# metagenomics data
# Kraken2/Bracken feature table reduced to the top 350 taxa based on median counts 
# feature table normalized according to DESeq2 defaults - relative log expression

# metabolomics data
# missing/zero values replace by 1/5 of the minimum positive value for each variable
# 40% filtered out by interquartile range
# normalization by median plus log10 transformation

# set working directory
setwd("/N/scratch/kmvanden/diablo")

#load libraries
library(mixOmics)


# read in feature tables and metadata
#check that sample order is the same
metabo = read.table("metabo_filt_norm.txt", header = TRUE)
metabo = as.matrix(metabo) # convert to matrix
metabo = t(metabo) # transpose matrix
class(metabo)
dim(metabo)

metagen = read.table("deseq2_norm_med_red_table.txt", header = TRUE)
metagen = as.matrix(metagen) # convert to matrix
metagen = t(metagen) # transpose the matrix
class(metagen)
dim(metagen)

Y = read.table("metadata_red.txt", header = TRUE)
sum(row.names(metabo)==row.names(Y))
sum(row.names(metagen)==row.names(Y))
Y = Y$outcome
Y = as.factor(Y)
class(Y)
summary(Y)

# store metabolomic and metagenomic data as a list
X = list(metabo = metabo, metagen = metagen)


### creation of design matrix (design)
# perform PLS with one component 
# calculate cross-correlations between components associated with each data set
res.pls.mali = pls(X$metagen, X$metabo, ncomp = 1)
cor(res.pls.mali$variates$X, res.pls.mali$variates$Y)
#           comp1
# comp1 0.8789482

design = matrix(0.88, ncol = length(X), nrow = length(X), dimnames = list(names(X), names(X)))
diag(design) = 0


### choose number of components to use
# assess global performance of model without variable selection
# run perf() with 8-fold cross-validation repeated 100 times for 5 components
diablo.mali = block.plsda(X, Y, ncomp = 5, design = design)
set.seed(20230711)
perf.diablo.mali = perf(diablo.mali, validation = "Mfold", folds = 8, nrepeat = 100)

pdf("numb_comp.pdf", width = 5, height = 5)
plot(perf.diablo.mali)
dev.off()

# set ncomp value - optimal number of components by prediction distance and error rate
ncomp = perf.diablo.mali$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# error rate is at a minimum with two components


### variable selection for sparse method
# estimation of classification error rate with respect to the number of variables
# run tune.block.splsda() with 8-fold cross-validation repeated 100 times for 1 component
test.keepX = list(metagen = seq(20,340,20), metabo = seq(20,780,20))
set.seed(20230711)
tune.diablo.mali = tune.block.splsda(X, Y, ncomp = ncomp,
                                     test.keepX = test.keepX, design = design,
                                     validation = "Mfold", folds = 8, nrepeat = 100,
                                     BPPARAM = BiocParallel::MulticoreParam(workers = 5),
                                     dist = "centroids.dist")

list.keepX = tune.diablo.mali$choice.keepX
# $metabo
# [1] 60 20
# 
# $metagen
# [1] 20 20


### final model and performance
# optimal number of components 
ncomp = 2

# optimal number of variables
list.keepX = list(metabo = c(60,20), metagen = c(20,20))

# final model
diablo.mali = block.splsda(X, Y, ncomp = ncomp,
                           keepX = list.keepX, design = design)

# model performance
# performance of model assessed with perf() using 8-fold cross-validation repeated 100 times
perf.diablo.mali = perf(diablo.mali, folds = 8, validation = "Mfold",
                        dist = "centroids.dist", nrepeat = 100)

perf.diablo.mali$error.rate # error rate for each component
# $metabo
# centroids.dist
# comp1      0.4079167
# comp2      0.4095833
# 
# $metagen
# centroids.dist
# comp1      0.2950000
# comp2      0.3422917

perf.diablo.mali$error.rate.per.class # error rate per class
# $metabo
# $metabo$centroids.dist
# comp1     comp2
# resistant   0.4937500 0.5037500
# susceptible 0.3220833 0.3154167
# 
# 
# $metagen
# $metagen$centroids.dist
# comp1     comp2
# resistant   0.3966667 0.4258333
# susceptible 0.1933333 0.2587500

# variables and their loading weights
selectVar(diablo.mali, block = "metagen", comp = 1)
selectVar(diablo.mali, block = "metagen", comp = 2)
selectVar(diablo.mali, block = "metabo", comp = 1)
selectVar(diablo.mali, block = "metabo", comp = 2)


# check correlation between metagenomics and metabolomics data
plotDiablo(diablo.mali, ncomp = 1) 
plotDiablo(diablo.mali, ncomp = 2)


# sample plot for each component separately
pdf("sample_plot_comp1_2.pdf", width = 10, height = 5)
plotIndiv(diablo.mali, ind.names = FALSE, legend = TRUE,
          X.label = "component 1", Y.label = "component 2",
          ellipse = TRUE, ellipse.level = 0.95, # 95% confidence ellipses
          col = c("dodgerblue3", "firebrick3"), # color of the groups
          cex = 3, pch = c(1,1), point.lwd = 0.5 # size, shape and line width of points
          )
dev.off()


# correlation circle plots
pdf("plot_corr_circle_noname.pdf", width = 5, height = 5)
plotVar(diablo.mali, var.names = FALSE, 
        style = "graphics", legend = TRUE,
        pch = c(1,1), cex = c(0.5,0.5), # shape and size of the points
        col = c("deepskyblue3", "mediumorchid4"), # color of the points
        cutoff = 0.30 # correlation coefficient cut-off
)
dev.off()

pdf("plot_corr_circle_large.pdf", width = 30, height = 30)
plotVar(diablo.mali, var.names = TRUE, # to see names of metabolites and taxa clearly
        style = "graphics", legend = TRUE,
        pch = c(1,1), cex = c(0.5,0.5), # shape and size of the points
        col = c("deepskyblue3", "mediumorchid4"), # color of the points
        cutoff = 0.30 # correlation coefficient cut-off
        )
dev.off()


# plot loading weights
pdf("plot_loadings_metagen_comp1.pdf", width = 10, height = 10)
plotLoadings(diablo.mali, 
             comp = 1, # component of interest
             contrib = "max", # color of bar corresponds to group with the maximal expression levels
             method = "median", # criterion to asses the contribution
             ndisplay = 20, # how many of the most important variables are plotted
             block = "metagen", # which block to show
             legend.color = c("dodgerblue3", "firebrick3") # color of the outcomes
             )
dev.off()

pdf("plot_loadings_metagen_comp2.pdf", width = 10, height = 10)
plotLoadings(diablo.mali, 
             comp = 2, # component of interest
             contrib = "max", # color of bar corresponds to group with the maximal expression levels
             method = "median", # criterion to asses the contribution
             ndisplay = 20, # how many of the most important variables are plotted
             block = "metagen", # which block to show 
             legend.color = c("dodgerblue3", "firebrick3") # color of the outcomes
)
dev.off()


pdf("plot_loadings_metabo_comp1.pdf", width = 10, height = 10)
plotLoadings(diablo.mali, 
             comp = 1, # component of interest
             contrib = "max", # color of bar corresponds to group with the maximal expression levels
             method = "median", # criterion to asses the contribution
             ndisplay = 40, # how many of the most important variables are plotted
             block = "metabo", # which block to show
             legend.color = c("dodgerblue3", "firebrick3") # color of the outcomes
)
dev.off()

pdf("plot_loadings_metabo_comp2.pdf", width = 10, height = 10)
plotLoadings(diablo.mali, 
             comp = 2, # component of interest
             contrib = "max", # color of bar corresponds to group with the maximal expression levels
             method = "mean", # criterion to asses the contribution
             ndisplay = 20, # how many of the most important variables are plotted
             block = "metabo", # which block to show 
             legend.color = c("dodgerblue3", "firebrick3") # color of the outcomes
)
dev.off()


### plot relevance network
# taxa names shortened to 10 characters or fewer for better visualization in the relevance network
# metabolites given 6 digit code for better visualization in the relevance network

# read in feature tables with shortened names
# check that sample order is the same
metabo_short = read.table("metabo_filt_norm_short.txt", header = TRUE)
metabo_short = as.matrix(metabo_short) # convert to matrix
metabo_short = t(metabo_short) # transpose matrix
class(metabo_short)
dim(metabo_short)

metagen_short = read.table("deseq2_norm_med_red_short.txt", header = TRUE)
metagen_short = as.matrix(metagen_short) # convert to matrix
metagen_short = t(metagen_short) # transpose the matrix
class(metagen_short)
dim(metagen_short)

sum(row.names(metabo_short)==row.names(metagen_short))

# update X variable
X_short = list(metabo_short = metabo_short, metagen_short = metagen_short)

# update matrix
design_short = matrix(0.88, ncol = length(X_short), nrow = length(X_short), dimnames = list(names(X_short), names(X_short)))
diag(design_short) = 0

# update components 
ncomp_short = 2

# update variable selection
list.keepX_short = list(metabo_short = c(60,20), metagen_short = c(20,20))

# final model
diablo.mali_short = block.splsda(X_short, Y, ncomp = ncomp_short,
                                 keepX = list.keepX_short, design = design_short)

# plot relevance network
network(diablo.mali_short, blocks = c(1,2),
        # save = "pdf", name.save = "plot_network_short_corr035_comp1_2",
        color.node = c("lightskyblue", "plum3"), # color of the metabo_short and metagen_short nodes
        cex.node.name = 0.7, # font size for node labels
        block.var.names = TRUE, # FALSE = remove node names
        # interactive = TRUE, # interactive correlation coefficient scroll bar
        cutoff = 0.35, # correlation coefficient cutoff
        color.edge = c("red4","red3","red2","green2","green3","green4"), # color of edges (represents positive or negative correlations)
        lwd.edge = 1.5, # width of the edges
        comp = list(metagen_short=c(1,2),metabo_short=c(1,2)) # which components to include
)


# AUROC plots
auc.diablo.mali.metabo = auroc(diablo.mali, roc.block = 1, roc.comp = 2, print = TRUE)
# AUC = 0.724
# p-value = 0.007815
auc.diablo.mali.metagen = auroc(diablo.mali, roc.block = 2, roc.comp = 2, print = TRUE)
# AUC = 0.7778
# p-value = 0.0009698


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
#   [1] mixOmics_6.20.0 ggplot2_3.4.3   lattice_0.21-8  MASS_7.3-60    
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10         RSpectra_0.16-1     plyr_1.8.8          pillar_1.9.0        compiler_4.2.1      RColorBrewer_1.1-3  tools_4.2.1        
# [8] lifecycle_1.0.3     tibble_3.2.1        gtable_0.3.3        pkgconfig_2.0.3     rlang_1.1.1         Matrix_1.5-4.1      igraph_1.4.3       
# [15] cli_3.6.1           rstudioapi_0.15.0   ggrepel_0.9.3       parallel_4.2.1      gridExtra_2.3       stringr_1.5.0       withr_2.5.0        
# [22] dplyr_1.1.2         generics_0.1.3      vctrs_0.6.3         grid_4.2.1          tidyselect_1.2.0    glue_1.6.2          ellipse_0.4.5      
# [29] R6_2.5.1            fansi_1.0.4         rARPACK_0.11-0      BiocParallel_1.30.4 farver_2.1.1        purrr_1.0.2         tidyr_1.3.0        
# [36] reshape2_1.4.4      corpcor_1.6.10      magrittr_2.0.3      scales_1.2.1        codetools_0.2-19    matrixStats_1.0.0   colorspace_2.1-0   
# [43] labeling_0.4.2      utf8_1.2.3          stringi_1.7.12      munsell_0.5.0      

