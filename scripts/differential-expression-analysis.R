library(tidyverse)
library(DESeq2)
library(apeglm)
library(ggrepel)
library(BiocParallel)
library(pheatmap)

###############################
# DESeq differential expression analysis of filtered RNAseq samples
# Filtered to remove singleton samples, and those that have a mapping rate > 97%
# Just under 200 samples remaining after these filters
###############################

# qced sample map 
qced_sample_map <- read_tsv("metadata/rnaseq/2026-10-06-qced-sample-map.tsv") %>% 
  mutate(sample = paste0(sample, "_count")) %>% 
  select(sample, batch, sample_type)

# qced samples
qced_samples <- qced_sample_map %>% 
  pull(sample)

# count table - select for only the filtered, qced samples
count_table <- read_tsv("raw_data/rnaseq/2026-09-25-combined-expression-matrix.tsv") %>% 
  select(gene_id, gene_name, gene_biotype, all_of(qced_samples))

# check that samples in the map and count table match
map_samples <- qced_sample_map$sample
counts_samples <- setdiff(colnames(count_table), c("gene_id", "gene_name", "gene_biotype"))
setdiff(map_samples, counts_samples)
setdiff(counts_samples, map_samples)

# count matrix
count_matrix <- count_table %>% 
  pivot_longer(
    !c(gene_id, gene_name, gene_biotype),
    names_to = "sample",
    values_to = "count"
  ) %>% 
  mutate(count = round(count)) %>% 
  select(gene_id, sample, count) %>% 
  pivot_wider(names_from = sample, values_from = count) %>% 
  column_to_rownames("gene_id") %>% 
  as.matrix

sample_metadata <- qced_sample_map %>% 
  distinct(sample, sample_type, batch) %>% 
  dplyr::slice(match(colnames(count_matrix), sample)) %>% 
  column_to_rownames("sample") %>% 
  mutate(sample_type = factor(make.names(sample_type)),
         batch = factor(batch))

# check metadata and counts matrix match
all(rownames(sample_metadata) == colnames(count_matrix))

# LPS+ positive control treatment as reference for comparisons
# Testing how the fermented food extracts differ from the LPS+ condition
sample_metadata$sample_type <- relevel(sample_metadata$sample_type,
                                       ref="positive.control")

# DESeq object
# Design to test for effects of each sample_type (fermented food), and control effects of each batch
dds <- DESeqDataSetFromMatrix(count_matrix,
                              colData = sample_metadata,
                              design = ~ sample_type + batch)

# remove low counts
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# sample QC checks and compare sample distances
# compare log2(x+1) vs vst transformations of gene counts
vsd <- vst(dds, blind = FALSE)

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)   # vst upward shift for lower count values, on log scale genes with low values are extremely variable, vst accounts for differences between samples more

# plot sample distance matrix
sampleDists <- dist(t(assay(vsd)))

