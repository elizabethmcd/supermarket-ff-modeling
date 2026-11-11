library(tidyverse)
library(DESeq2)

###############################
# Unfermented vs fermented expression profiles
###############################

#################################################
# read in objects and metadata
#################################################

# metadata for unfermented vs fermented pairs
pairs <- read_tsv("metadata/rnaseq/unfermented-vs-fermented-pairs.tsv")

# RDS objects
dds <- readRDS("results/rnaseq/dds_fitted.rds")
gene_info <- readRDS("results/rnaseq/gene_info.rds")
sample_codes <- readRDS("results/rnaseq/sample_codes.rds")




