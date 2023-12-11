# Network Analysis
# NetCoMi
# Kraken2/Bracken feature table reduced to the top 350 taxa based on median counts


# set working directory
setwd("/N/scratch/kmvanden/bracken_netcomi")


# load libraries
library("NetCoMi")
library("phyloseq")


# read in feature table and metadata
# check that sample order is the same
otu_med = read.table("otu_table_med.txt")
otu_med = as.matrix(otu_med) # convert to matrix
otu_med = t(otu_med) # transpose matrix
dim(otu_med)

metadata = read.table("metadata.txt")
dim(metadata)
sum(row.names(otu_med)==row.names(metadata))

group = metadata$outcome
group = as.factor(group)


### SPARCC ###

### network construction

sparcc_network_med = netConstruct(otu_med, group = group,
                                  measure = "sparcc", 
                                  zeroMethod = "none", # zero method is included in the measure for SparCC
                                  normMethod = "none", # no normalization method is used for SparCC
                                  sparsMethod = "t-test", 
                                  seed = 20230711)

# get edge list
write.csv(sparcc_network_med$edgelist1, "sparcc_network_med_edges_res.csv", row.names = TRUE)
write.csv(sparcc_network_med$edgelist2, "sparcc_network_med_edges_sus.csv", row.names = TRUE)


### network analysis - degree centrality

sparcc_analysis_med_deg = netAnalyze(sparcc_network_med, 
                                     clustMethod = "cluster_fast_greedy", 
                                     hubPar = "degree", 
                                     hubQuant = 0.95) 

sparcc_analysis_summary_med_deg = summary(sparcc_analysis_med_deg, groupNames = c("resistant", "susceptible"))
sparcc_analysis_summary_med_deg_hubs = sparcc_analysis_summary_med_deg$hubs
write.csv(sparcc_analysis_summary_med_deg_hubs, "sparcc_med_deg_hubs.csv", row.names = TRUE)


### network plots - degree centrality

pdf("sparcc_plot_med_deg_total.pdf", width = 60, height = 30)

plot(sparcc_analysis_med_deg,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     # nodeFilter = "highestDegree", nodeFilterPar = 100,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

pdf("sparcc_plot_med_deg_node100.pdf", width = 60, height = 30)

plot(sparcc_analysis_med_deg,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestDegree", nodeFilterPar = 100,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

pdf("sparcc_plot_med_deg_node50.pdf", width = 60, height = 30)

plot(sparcc_analysis_med_deg,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestDegree", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()


### network comparison with permutation - degree centrality

sparcc_compare_med_deg = netCompare(sparcc_analysis_med_deg,
                                    permTest = TRUE, nPerm = 1000,
                                    adjust = "adaptBH", trueNullMethod = "convest",
                                    storeAssoPerm = TRUE, fileStoreAssoPerm = "sparcc_assoperm_med_deg",
                                    verbose = TRUE, seed = 20230711)

### summary of network comparison

sparcc_compare_med_deg_summary = summary(sparcc_compare_med_deg,
                                         groupNames = c("resistant", "susceptible"))

sparcc_taxa_deg_res = sparcc_compare_med_deg_summary$properties$deg1 # degree centrality measures for the resistant network
write.csv(sparcc_taxa_deg_res, "sparcc_taxa_deg_res.csv", row.names = TRUE)

sparcc_taxa_deg_sus = sparcc_compare_med_deg_summary$properties$deg2 # degree centrality measures for the susceptible network
write.csv(sparcc_taxa_deg_sus, "sparcc_taxa_deg_sus.csv", row.names = TRUE)

sparcc_taxa_clust_res = sparcc_compare_med_deg_summary$properties$clust1 # cluster ids for the resistant network 
write.csv(sparcc_taxa_clust_res, "sparcc_taxa_clust_res.csv", row.names = TRUE)

sparcc_taxa_clust_sus = sparcc_compare_med_deg_summary$properties$clust2 # cluster ids for the susceptible network 
write.csv(sparcc_taxa_clust_sus, "sparcc_taxa_clust_sus.csv", row.names = TRUE)


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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] phyloseq_1.42.0 NetCoMi_1.1.0   SpiecEasi_1.1.2
# 
# loaded via a namespace (and not attached):
#   [1] filematrix_1.3         VGAM_1.1-7             colorspace_2.0-2       dynamicTreeCut_1.63-1  htmlTable_2.2.1       
# [6] corpcor_1.6.9          XVector_0.32.0         base64enc_0.1-3        rstudioapi_0.15.0      lavaan_0.6-9          
# [11] bit64_4.0.5            AnnotationDbi_1.56.2   fansi_0.5.0            mvtnorm_1.1-2          codetools_0.2-18      
# [16] splines_4.1.1          mnormt_2.0.2           doParallel_1.0.16      impute_1.68.0          cachem_1.0.6          
# [21] rootSolve_1.8.2.2      knitr_1.34             glasso_1.11            ade4_1.7-17            Formula_1.2-4         
# [26] jsonlite_1.8.7         WGCNA_1.72-1           cluster_2.1.2          GO.db_3.14.0           png_0.1-7             
# [31] compiler_4.1.1         httr_1.4.7             backports_1.2.1        Matrix_1.5-4           fastmap_1.1.0         
# [36] cli_3.6.1              htmltools_0.5.5        tools_4.1.1            igraph_1.4.0           gtable_0.3.0          
# [41] glue_1.6.2             GenomeInfoDbData_1.2.7 reshape2_1.4.4         dplyr_1.1.2            Rcpp_1.0.11           
# [46] Biobase_2.58.0         vctrs_0.6.3            Biostrings_2.60.2      SPRING_1.0.4           rhdf5filters_1.6.0    
# [51] multtest_2.50.0        preprocessCore_1.56.0  ape_5.5                nlme_3.1-153           iterators_1.0.13      
# [56] psych_2.1.6            xfun_0.26              fastcluster_1.2.3      stringr_1.5.0          rbibutils_2.2.3       
# [61] lifecycle_1.0.3        irlba_2.3.3            gtools_3.9.2           zlibbioc_1.38.0        MASS_7.3-54           
# [66] scales_1.2.1           orca_1.1-1             doSNOW_1.0.19          parallel_4.1.1         biomformat_1.22.0     
# [71] rhdf5_2.38.1           huge_1.3.5             RColorBrewer_1.1-2     pbapply_1.4-3          memoise_2.0.1         
# [76] gridExtra_2.3          ggplot2_3.4.2          rpart_4.1-15           latticeExtra_0.6-29    stringi_1.7.8         
# [81] RSQLite_2.2.15         S4Vectors_0.30.2       corrplot_0.90          pcaPP_1.9-74           foreach_1.5.1         
# [86] checkmate_2.1.0        permute_0.9-7          BiocGenerics_0.38.0    BiocParallel_1.26.2    shape_1.4.6           
# [91] GenomeInfoDb_1.30.1    Rdpack_2.1.2           rlang_1.1.1            pkgconfig_2.0.3        bitops_1.0-7          
# [96] matrixStats_0.60.1     lattice_0.20-44        Rhdf5lib_1.20.0        htmlwidgets_1.5.4      bit_4.0.4             
# [101] tidyselect_1.2.0       plyr_1.8.8             magrittr_2.0.3         R6_2.5.1               IRanges_2.26.0        
# [106] snow_0.4-4             generics_0.1.0         Hmisc_4.5-0            mixedCCA_1.6.2         DBI_1.1.1             
# [111] pillar_1.9.0           foreign_0.8-81         mgcv_1.8-36            abind_1.4-5            survival_3.2-13       
# [116] KEGGREST_1.34.0        RCurl_1.98-1.4         nnet_7.3-16            tibble_3.2.1           pulsar_0.3.10         
# [121] crayon_1.4.1           fdrtool_1.2.16         utf8_1.2.2             tmvnsim_1.0-2          jpeg_0.1-9            
# [126] qgraph_1.9.3           grid_4.1.1             pbivnorm_0.6.0         data.table_1.14.8      blob_1.2.2            
# [131] vegan_2.6-4            digest_0.6.27          stats4_4.1.1           munsell_0.5.0          glmnet_4.1-2  

