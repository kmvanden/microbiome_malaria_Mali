# preprocessing of mouse 16S rRNA sequencing feature table

# set working directory
setwd("/N/scratch/kmvanden/16S_mouse")

# read in libraries
library(dplyr)

# read in feature table 
otu = read.table("otu_table.txt", header = TRUE, row.names = NULL)

# read in taxonomy file
taxonomy = read.table("taxonomy.txt", header = TRUE, sep = "\t")

# read in the metadata
metadata = read.table("metadata.txt", row.names = NULL)
metadata = metadata[order(metadata$outcome),]
resistant = metadata$row.names[1:20]
susceptible = metadata$row.names[21:32]

# merge tables by otu_id and replace otu_id by taxa name
otu_table = merge(otu, taxonomy, by = "otu_id", all.x = TRUE)
otu_table[,1] = otu_table[,34]
otu_table = otu_table[,-34]

# shorten and format taxa names to genus_species
otu_table[,1] = gsub(".*s__", "", otu_table[,1])
otu_table[,1] = gsub(";+", "", otu_table[,1])
otu_table[,1] = gsub("t__", "", otu_table[,1])
otu_table[,1] = gsub(" ", "_", otu_table[,1])
otu_table[,1] = gsub("_sp_", "_sp._", otu_table[,1])

# set taxa names as row names
otu_table_row = otu_table[,-1]
rownames(otu_table_row) = otu_table[,1]
otu_table = otu_table_row

# calculate median counts of rows for resistant and susceptible mice and sort
otu_table$med_res = apply(otu_table[,resistant], 1, median)
otu_table$med_sus = apply(otu_table[,susceptible], 1, median)
otu_table_ordered = otu_table[order(otu_table$med_res, otu_table$med_sus, decreasing = TRUE),]

# retain only taxa that have a median count greater than zero in either the resistant or the susceptible group
otu_table_full_med = otu_table_ordered[otu_table_ordered$med_res > 0 | otu_table_ordered$med_sus > 0,]

dim(otu_table_full_med)

# remove the median column
otu_table_full_med = otu_table_full_med[,-34]
otu_table_full_med = otu_table_full_med[,-33]

# export table
write.csv(as.data.frame(otu_table_full_med), file = "otu_table_full_med.csv")


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
# other attached packages:
#   [1] dplyr_1.1.2
# 
# loaded via a namespace (and not attached):
#   [1] fansi_1.0.3      withr_2.5.0      utf8_1.2.2       R6_2.5.1         lifecycle_1.0.3  magrittr_2.0.3   pillar_1.9.0     rlang_1.1.1     
# [9] cli_3.6.2        rstudioapi_0.14  vctrs_0.6.3      generics_0.1.3   tools_4.2.1      glue_1.6.2       compiler_4.2.1   pkgconfig_2.0.3 
# [17] tidyselect_1.2.0 tibble_3.2.1 

