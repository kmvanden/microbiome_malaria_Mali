# preprocessing of shotgun metagenomics feature table

# set working directory
setwd("/N/scratch/kmvanden/bracken")

# read in feature table 
otu = read.table("feature_table.txt", header = TRUE, row.names = NULL)

# remove bracken report suffix
colnames(otu) = gsub(".breport", "", colnames(otu))

# shorten taxa names to genus_species
otu[,1] = gsub(".*s__", "",otu[,1])

# set taxa names as row names
otu_table = otu[,-1]
rownames(otu_table) = otu[,1]

# calculate median counts of rows and sort 
otu_table$median = apply(otu_table, 1, median)
otu_table = otu_table[order(otu_table$median, decreasing = TRUE),]

# retain only the top 350 taxa by median counts
otu_table = otu_table[1:350,]

# remove the median column
otu_table = otu_table[,-71]

# export table
write.csv(as.data.frame(otu_table), file = "otu_table_med.csv")


### SESSION INFO ###
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux 8.8 (Ootpa)
# 
# Matrix products: default
# BLAS:   /geode2/soft/hps/rhel8/r/gnu/4.2.1_X11/lib64/R/lib/libRblas.so
# LAPACK: /geode2/soft/hps/rhel8/r/gnu/4.2.1_X11/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.2.1  tools_4.2.1     rstudioapi_0.14

