# Differential Abundance Analysis 
# DESeq2 | Corncob | MaAsLin2 | ALDEx2
# feature table reduced to taxa with a median value greater than zero in either resistant or susceptible group


# set working directory
setwd("/N/scratch/kmvanden/mouse16S_DAA")

### DESEQ2 ###

# load libraries
library("DESeq2")
library("apeglm")
library("ashr")

# read in feature table and metadata
# check that sample order is the same
otu = read.table("otu_table_full_med.txt") 
dim(otu)

metadata = read.table("metadata.txt")
sum(colnames(otu)==rownames(metadata))


### DESEQ2 ANALYSIS

# create a DESeq2DataSet object
dds = DESeqDataSetFromMatrix(countData = otu,
                             colData = metadata,
                             design = ~ outcome)

dds # summary of DESeq2DataSet object
 
# differential expression analysis
# dispersion trend was not well captured by fitType "parametric", and was changed to local regression fit
dds = DESeq(dds, fitType = "local")

# extract results from DESeq analysis
res = results(dds)
# order results by padj value
resOrdered = res[order(res$padj),]

sum((resOrdered$padj < 0.05)[!is.na((resOrdered$padj < 0.05))]) # 30 differentially abundant taxa (padj < 0.05)
sum((resOrdered$padj < 0.10)[!is.na((resOrdered$padj < 0.10))]) # 34 differentially abundant taxa (padj < 0.10)


### MA PLOTS
# scatter plot of log2 fold changes versus the average expression signal
plotMA(res, ylim=c(-20,20))

# original DESeq2 shrinkage estimator (adaptive Normal distribution as prior)
resLFC_normal = lfcShrink(dds, coef=2, type="normal")
plotMA(resLFC_normal, main="normal", ylim = c(-10,10))

# adaptive t prior shrinkage estimator (from apeglm package)
resLFC_apeglm = lfcShrink(dds, coef=2, type="apeglm")
plotMA(resLFC_apeglm, main="apeglm", ylim=c(-10,10))

# adaptive shrinkage estimator (from ashr package)
resLFC_ashr = lfcShrink(dds, coef=2, type="ashr")
plotMA(resLFC_ashr, main="ashr", ylim=c(-20,20))


### INDIVIDUAL PLOTS
# normalizes counts by estimated size factors + adds pseudocount of 1/2 for log scale plotting
plotCounts(dds, gene = "Bacteroides_intestinalis", intgroup = "outcome")
plotCounts(dds, gene = "Bacteroides_cellulosilyticus", intgroup = "outcome")
plotCounts(dds, gene = "Bacteroides_ovatus", intgroup = "outcome")
plotCounts(dds, gene = "Blautia_faecis", intgroup = "outcome")
plotCounts(dds, gene = "Turicibacter_sp._LA61", intgroup = "outcome")
plotCounts(dds, gene = "Gemmiger_formicilis", intgroup = "outcome")


### EXPORTING DATA | fold change and padj values
# exporting log2 fold change and padj to csv files 
write.csv(as.data.frame(resOrdered), file = "mouse16S_deseq2_fold_change_padj.csv")


### EXPORTING DATA | normalized counts
dds = estimateSizeFactors(dds)
ddsCounts = counts(dds, normalized = TRUE)
write.csv(as.data.frame(ddsCounts), file = "mouse16S_deseq2_counts.csv")


### CORNCOB ###

# load libraries 
library("phyloseq")
library("corncob")
library("janitor")
library("tidyverse") 

# read in feature table and metadata
# check that sample order is the same
otu = read.table("otu_table_full_med.txt") 
dim(otu)

metadata = read.table("metadata.txt")
sum(colnames(otu)==rownames(metadata))

# build otu_table and sample_data for use in phyloseq
otu = otu_table(otu, taxa_are_rows = TRUE)
metadata = sample_data(metadata)


### CORNCOB ANALYSIS

# create phyloseq object
phylo = phyloseq(otu, metadata)
phylo # summary of phyloseq object


### TESTING AND PLOTTING A SINGLE TAXA
# fit our model to a single taxa with the covariate "outcome"
corncob = bbdml(formula = Bacteroides_intestinalis ~ outcome,
                phi.formula = ~ outcome,
                data = phylo)

# plot the data a relative abundance scale 
plot(corncob, total = TRUE, B = 1000, color = "outcome")

summary(corncob) # summary of model output


### TESTING AND PLOTTING ALL TAXA
# test and plot all taxa at a padj of 0.05
corncob_all005 = differentialTest(formula = ~ outcome,
                                  phi.formula = ~ outcome,
                                  formula_null = ~ 1,
                                  phi.formula_null = ~ outcome,
                                  test = "Wald", boot = FALSE,
                                  data = phylo,
                                  fdr_cutoff = 0.05)

# list of all significant taxa (3 at padj < 0.05)
corncob_all005$significant_taxa

# lists the model and associated coefficients for significant taxa
corncob_all005$significant_models[[1]]

# plot the model coefficients of our results for significant taxa
plot(corncob_all005)

#test and plot all taxa at a padj of < 0.10
corncob_all010 = differentialTest(formula = ~ outcome,
                                  phi.formula = ~ outcome,
                                  formula_null = ~ 1,
                                  phi.formula_null = ~ outcome,
                                  test = "Wald", boot = FALSE,
                                  data = phylo,
                                  fdr_cutoff = 0.10)

# list of all significant taxa (7 at padj of 0.10)
corncob_all010$significant_taxa

# lists the model and associated coefficients for significant taxa
corncob_all010$significant_models[[1]]

# plot the model coefficients of our results for significant taxa
plot(corncob_all010)


### EXPORTING DATA | padj values
# https://taylorreiter.github.io/2022-08-29-From-raw-metagenome-reads-to-taxonomic-differential-abundance-with-sourmash-and-corncob/

corncob_df = function(corncob_results){
  df = data.frame()
  for(i in 1:length(corncob_results$significant_taxa)){
    taxa_name = corncob_results$significant_taxa[i]
    taxa_df = corncob_results$significant_models[[i]]$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("covariate") %>%
      mutate(taxa_name = taxa_name)
    df = bind_rows(df, taxa_df)
  }
  df = clean_names(df) %>%
    select(covariate, estimate, std_error, t_value, fdr = pr_t, taxa_name)
  return(df)
}

# use the function to filter taxa with padj < 0.05
corncob_padj005 = corncob_df(corncob_all005) %>%
  filter(fdr < 0.05) %>%
  filter(!grepl("Intercept", covariate))
write.csv(as.data.frame(corncob_padj005), file = "mouse16S_corncob_padj005.csv")

# use the function to filter taxa with padj < 0.10
corncob_padj010 = corncob_df(corncob_all010) %>%
  filter(fdr < 0.10) %>%
  filter(!grepl("Intercept", covariate))
write.csv(as.data.frame(corncob_padj010), file = "mouse16S_corncob_padj010.csv")


### MAASLIN2 ###

# load libraries
library("Maaslin2")


# read in feature table and metadata
# check that sample order is the same
otu = read.table("otu_table_full_med.txt") 
dim(otu)

metadata = read.table("metadata.txt")
sum(colnames(otu)==rownames(metadata))


### CREATE MAASLIN2 OBJECT
# significant features with padj < 0.05
fit_data005 = Maaslin2(input_data = otu, input_metadata = metadata, output = "maaslin2_output005",
                       normalization = "TSS", transform = "LOG", standardize = FALSE, max_significance = 0.05, 
                       fixed_effects = "outcome")

# significant features with padj < 0.10
fit_data010 = Maaslin2(input_data = otu, input_metadata = metadata, output = "maaslin2_output010",
                       normalization = "TSS", transform = "LOG", standardize = FALSE, max_significance = 0.10, 
                       fixed_effects = "outcome")


### EXPORTING DATA | normalized and transformed counts
# code for TSS normalization and LOG transformation

# load required packages
for (lib in c('vegan', 'chemometrics', 'car', 'metagenomeSeq', 'edgeR')) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# normalization
normalizeFeatures = function(features, normalization) {
  if (normalization == 'TSS')
  {
    features = TSSnorm(features)
  }
  
  if (normalization == 'CLR')
  {
    features = CLRnorm(features)
  }
  
  if (normalization == 'CSS')
  {
    features = CSSnorm(features)
  }
  
  if (normalization == 'TMM')
  {
    features = TMMnorm(features)
  }
  
  if (normalization == 'NONE')
  {
    features = features
  }
  
  return(features)
}

# apply TSS normalization to the dataset
TSSnorm = function(features) {
  # convert to matrix from data frame
  features_norm = as.matrix(features)
  dd = colnames(features_norm)
  
  # TSS normalization
  features_TSS =
    vegan::decostand(
      features_norm,
      method = "total",
      MARGIN = 1,
      na.rm = TRUE)
  
  # convert back to data frame
  features_TSS = as.data.frame(features_TSS)
  
  # rename the true positive features
  colnames(features_TSS) <- dd
  
  
  # return
  return(features_TSS)
}

# log transformation
LOG = function(x) {
  y = replace(x, x == 0, min(x[x>0]) / 2)
  return(log10(y))
}

# transformation
transformFeatures = function(features, transformation) {
  if (transformation == 'LOG')     {
    features = apply(features, 2, LOG)
  }
  
  if (transformation == 'LOGIT')     {
    features = apply(features, 2, LOGIT)
  }
  
  if (transformation == 'AST')     {
    features = apply(features, 2, AST)
  }
  
  return(features)
}

otu_tss = normalizeFeatures(features = otu, normalization = 'TSS') # TSS normalization
otu_tss_log = transformFeatures(features = otu_tss, transformation = 'LOG') # LOG transformation
write.csv(as.data.frame(otu_tss_log), file = "mouse16S_maaslin2_tss_log_counts.csv")


### ALDEX2 ###

# load libraries
library("ALDEx2")

# read in feature table and metadata
# check that sample order is the same
otu = read.table("otu_table_full_med.txt") 
dim(otu)

metadata = read.table("metadata.txt")
sum(colnames(otu)==rownames(metadata))

# convert metadata into a vector of the outcome groups
metadata = metadata$outcome

### ALDEX2 ANALYSIS

### create ALDEx2 object and perform differential abundance analysis
ald = aldex(otu, metadata)

# ALDEx2 object ordered by Wilcox corrected p-value
aldOrdered = ald[order(ald$wi.eBH),]


### MA AND MW PLOTS
# Bland-Altman (MA) plot and effect (MW) plot
par(mfrow = c(1,2))
aldex.plot(ald, type="MA", test="wilcox", xlab="log-ratio abundance", ylab="difference")
aldex.plot(ald, type="MW", test="wilcox", xlab="dispersion", ylab="difference")


### EFFECT SIZE AND VOLCANO PLOTS
# effect size plot
plot(ald$effect, ald$we.ep, log="y", cex=0.7, col=rgb(0,0,1,0.2),
     pch=19, xlab="Effect size", ylab="P value", main="Effect size plot")
points(ald$effect, ald$we.eBH, cex=0.7, col=rgb(1,0,0,0.2),
       pch=19)
abline(h=0.05, lty=2, col="grey")
legend(15,1, legend=c("P value", "BH-adjusted"), pch=19, col=c("blue", "red"))

# volcano plot
plot(ald$diff.btw, ald$we.ep, log="y", cex=0.7, col=rgb(0,0,1,0.2),
     pch=19, xlab="Difference", ylab="P value", main="Volcano plot")
points(ald$diff.btw, ald$we.eBH, cex=0.7, col=rgb(1,0,0,0.2),
       pch=19)
abline(h=0.05, lty=2, col="grey")


### EXPORTING DATA | effect size and padj values
write.csv(as.data.frame(ald), file = "mouse16S_aldex2_effect_padj.csv")
ald_clr = aldex.clr(otu, metadata, mc.samples = 128)
ald_clr = getMonteCarloInstances(ald_clr)

write.csv(as.data.frame(ald_clr), file = "mouse16S_aldex2_clr_counts.csv")


### session info ###
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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ALDEx2_1.26.0               zCompositions_1.4.0-1       truncnorm_1.0-8             NADA_1.6-1.1               
# [5] survival_3.2-13             MASS_7.3-54                 edgeR_3.36.0                metagenomeSeq_1.36.0       
# [9] RColorBrewer_1.1-2          glmnet_4.1-2                Matrix_1.5-4                limma_3.50.3               
# [13] car_3.0-12                  carData_3.0-4               chemometrics_1.4.4          rpart_4.1-15               
# [17] vegan_2.6-4                 lattice_0.20-44             permute_0.9-7               Maaslin2_1.8.0             
# [21] lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0               dplyr_1.1.2                
# [25] purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1               
# [29] ggplot2_3.4.2               tidyverse_2.0.0             janitor_2.2.0               corncob_0.3.1              
# [33] phyloseq_1.42.0             ashr_2.2-63                 apeglm_1.16.0               DESeq2_1.34.0              
# [37] SummarizedExperiment_1.24.0 Biobase_2.58.0              MatrixGenerics_1.6.0        matrixStats_0.60.1         
# [41] GenomicRanges_1.44.0        GenomeInfoDb_1.30.1         IRanges_2.26.0              S4Vectors_0.30.2           
# [45] BiocGenerics_0.38.0        
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.2.1          VGAM_1.1-7               plyr_1.8.8               igraph_1.4.0            
# [5] splines_4.1.1            BiocParallel_1.26.2      lpsymphony_1.22.0        invgamma_1.1            
# [9] foreach_1.5.1            lars_1.2                 SQUAREM_2021.1           fansi_0.5.0             
# [13] magrittr_2.0.3           checkmate_2.1.0          memoise_2.0.1            cluster_2.1.2           
# [17] tzdb_0.1.2               Biostrings_2.60.2        annotate_1.72.0          bdsmatrix_1.3-4         
# [21] timechange_0.2.0         colorspace_2.0-2         blob_1.2.2               crayon_1.4.1            
# [25] RCurl_1.98-1.4           jsonlite_1.8.7           genefilter_1.76.0        biglm_0.9-2.1           
# [29] iterators_1.0.13         ape_5.5                  glue_1.6.2               hash_2.2.6.3            
# [33] registry_0.5-1           gtable_0.3.0             zlibbioc_1.38.0          XVector_0.32.0          
# [37] DelayedArray_0.20.0      RcppZiggurat_0.1.6       Rhdf5lib_1.20.0          shape_1.4.6             
# [41] DEoptimR_1.0-9           abind_1.4-5              scales_1.2.1             mvtnorm_1.1-2           
# [45] som_0.3-5.1              DBI_1.1.1                Rcpp_1.0.11              xtable_1.8-4            
# [49] emdbook_1.3.13           mclust_5.4.7             proxy_0.4-26             bit_4.0.4               
# [53] httr_1.4.7               getopt_1.20.3            gplots_3.1.1             pkgconfig_2.0.3         
# [57] XML_3.99-0.7             farver_2.1.0             nnet_7.3-16              locfit_1.5-9.6          
# [61] utf8_1.2.2               tidyselect_1.2.0         labeling_0.4.2           rlang_1.1.1             
# [65] reshape2_1.4.4           AnnotationDbi_1.56.2     munsell_0.5.0            tools_4.1.1             
# [69] cachem_1.0.6             cli_3.6.1                generics_0.1.0           RSQLite_2.2.15          
# [73] pls_2.8-2                ade4_1.7-17              biomformat_1.22.0        fastmap_1.1.0           
# [77] bit64_4.0.5              robustbase_0.93-8        caTools_1.18.2           KEGGREST_1.34.0         
# [81] pbapply_1.4-3            nlme_3.1-153             slam_0.1-48              ROI_1.0-1               
# [85] compiler_4.1.1           rstudioapi_0.15.0        png_0.1-7                e1071_1.7-8             
# [89] lpSolveAPI_5.5.2.0-17.7  geneplotter_1.72.0       pcaPP_1.9-74             stringi_1.7.8           
# [93] multtest_2.50.0          vctrs_0.6.3              ROI.plugin.lpsolve_1.0-2 pillar_1.9.0            
# [97] lifecycle_1.0.3          trust_0.1-8              rhdf5filters_1.6.0       optparse_1.7.3          
# [101] data.table_1.14.8        bitops_1.0-7             irlba_2.3.3              R6_2.5.1                
# [105] detectseparation_0.3     KernSmooth_2.23-20       codetools_0.2-18         gtools_3.9.2            
# [109] Wrench_1.12.0            rhdf5_2.38.1             withr_2.5.0              GenomeInfoDbData_1.2.7  
# [113] mgcv_1.8-36              hms_1.1.3                grid_4.1.1               coda_0.19-4             
# [117] class_7.3-19             Rfast_2.0.7              snakecase_0.11.0         logging_0.10-108        
# [121] mixsqp_0.3-48            bbmle_1.0.25             numDeriv_2016.8-1.1

