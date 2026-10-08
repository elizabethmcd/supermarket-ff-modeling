library(tidyverse)
library(readr)
library(readxl)
library(DESeq2)
library(apeglm)
library(ggrepel)
library(BiocParallel)
library(pheatmap)
library(RColorBrewer)
library(glmpca)
library(genefilter)

###############################
# DESeq differential expression analysis of filtered RNAseq samples
# Filtered to remove singleton samples, and those that have a mapping rate > 97%
# Just under 200 samples remaining after these filters
###############################

# qced sample map 
qced_sample_map <- read_tsv("metadata/rnaseq/2026-10-06-qced-sample-map.tsv") %>% 
  mutate(sample = paste0(sample, "_count")) %>% 
  select(sample, batch, sample_code, sample_type)

# qced samples
qced_samples <- qced_sample_map %>% 
  pull(sample)

# count table - select for only the filtered, qced samples
count_table <- read_tsv("raw_data/rnaseq/2026-09-25-combined-expression-matrix.tsv") %>% 
  select(gene_id, gene_name, gene_biotype, all_of(qced_samples))

gene_info <- count_table %>% 
  select(gene_id, gene_name, gene_biotype) %>% 
  dplyr::rename(gene = gene_id)

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
  distinct(sample, sample_type, sample_code, batch) %>% 
  dplyr::slice(match(colnames(count_matrix), sample)) %>% 
  column_to_rownames("sample") %>% 
  mutate(sample_type = factor(make.names(sample_type)),
         batch = factor(batch))

sample_codes <- sample_metadata %>% 
  select(sample_type, sample_code) %>% 
  tibble::remove_rownames() %>% 
  distinct()
  

# check metadata and counts matrix match
all(rownames(sample_metadata) == colnames(count_matrix))

###############################
# Select samples with IFN bioactivity at certain cutoffs
###############################
bioactivity_df <- read_excel("raw_data/bioactivity/Fermented_food_data_mastersheet.xlsx", sheet = "full_data")

# samples that show some anti-inflammatory activity by IFN
bioactivity_ifn_filtered_codes <- bioactivity_df %>% 
  filter(Sample_type == "sample" | Sample_type == "Sample") %>% 
  dplyr::rename(sample_code = Sample_number) %>% 
  filter(normalized_viability > 70) %>% 
  filter(normalized_IFN <= 100) %>% 
  select(sample_code) %>% 
  unique() %>% 
  pull(sample_code)

# median ruxolitinib IFN activity
ruxolitinib_median_IFN <- bioactivity_df %>% 
  filter(Sample_number == "inhibitor2") %>% 
  filter(normalized_viability > 70) %>% 
  summarise(median(normalized_IFN)) %>% 
  pull(`median(normalized_IFN)`)

# samples that show at or below IFN compared to the median ruxolitinib IFN activity
bioactivity_ifn_filtered_ruxo_level_codes <- bioactivity_df %>% 
  filter(Sample_type == "sample") %>% 
  dplyr::rename(sample_code = Sample_number) %>% 
  filter(normalized_viability > 70) %>% 
  group_by(sample_code) %>% 
  summarise(median = median(normalized_IFN)) %>% 
  filter(median <= ruxolitinib_median_IFN) %>% 
  pull(sample_code)

###############################
# DESeq object and basic QC plots and checks
# check for batch effects between the sequencing runs
###############################
# LPS+ positive control treatment as reference for comparisons
# Testing how the fermented food extracts differ from the LPS+ condition
sample_metadata$sample_type <- relevel(sample_metadata$sample_type,
                                       ref="positive.control")
# DESeq object
# Design to test for effects of each sample_type (fermented food), and control effects of each batch
dds <- DESeqDataSetFromMatrix(count_matrix,
                              colData = sample_metadata,
                              design = ~ batch + sample_type) # ordered this way because looks at the last variable in the design

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
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- vsd$batch
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# plot PCA
plotPCA(vsd, intgroup = "batch")

# generalized principal component analysis
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$batch <- dds$batch
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = batch)) +
  geom_point(size =3) + 
  coord_fixed() + 
  ggtitle("glmpca - Generalized PCA")

# MDS plot from the VST data
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = batch)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

# remove batch variation in the PCA plot
mat <- assay(vsd)
mm <- model.matrix(~sample_type, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
pca <- plotPCA(vsd, intgroup = "batch", returnData = TRUE)
plotPCA(vsd, intgroup = "batch")

# pull metadata for certain clusters
pca$far <- pca$PC1 > 30
far_ids <- rownames(pca)[pca$far]
sample_map %>% 
  filter(sample %in% far_ids) # all negative controls cluster

pca$far <- pca$PC2 > 10
far_ids <- rownames(pca)[pca$far]
sample_map %>% 
  filter(sample %in% far_ids)

# do replicates of some samples in the top cluster appear in the main cluster of samples
far_types <- unique(colData(vsd)[far_ids, "sample_type"])
cd <- as.data.frame(colData(vsd))
cd$far <- rownames(cd) %in% far_ids
table(cd$sample_type[cd$sample_type %in% far_types], cd$far[cd$sample_type %in% far_types]) %>% 
  as.data.frame() %>% 
  dplyr::rename(sample_type = Var1, variable = Var2, frequency = Freq) %>% 
  pivot_wider(names_from = variable, values_from = frequency) %>% 
  filter(`TRUE` > 0 | `FALSE` > 0)

# check if replicates agree on PC2
pca <- plotPCA(vsd, intgroup = c("batch", "sample_type"), returnData = TRUE)
main <- pca[pca$PC1 < 30, ]
main$sample_type <- colData(vsd)[rownames(main), "sample_type"]
summary(aov(PC2 ~ sample_type, data = main))  # 93% of variance explained by sample type
summary(aov(PC2 ~ batch, data = main)) # 24% of variance explained by batch -> could be due to what foods ended up in which batch

# median within-sample type spread vs overall spread
within_sd <- tapply(main$PC2, main$sample_type, sd)
median(within_sd, na.rm = TRUE); sd(main$PC2) # replicates of same sample type differ by 1 unit on PC2, whole cluster on PC2 differs by 6

# variation due to sample type vs batch
summary(aov(PC2 ~ sample_type + batch, data = main))

# are replicates in different batches or the same batches
tab <- table(main$sample_type, main$batch)
sum(rowSums(tab > 0) > 1) 

# run the summaries on uncorrected VST data to compare to corrected
vsd_raw <- vst(dds, blind = FALSE)
pca_raw <- plotPCA(vsd_raw, intgroup = "batch", returnData = TRUE)
main_raw <- pca_raw[rownames(main), ]
main_raw$sample_type <- colData(vsd_raw)[rownames(main_raw), "sample_type"]

summary(aov(PC2 ~ sample_type + batch, data = main_raw)) # correction ends up removes a small, but significant batch effect, so keep batch in the ~design
# treat samples from different batches that have small effects carefully

###############################
# Differential expression analysis
###############################

# use BiocParallel for the DESeq analysis on dds object
register(MulticoreParam(workers = parallel::detectCores() - 2))
dds <- DESeq(dds, parallel = TRUE)

# compare all samples to the positive control
levs <- levels(dds$sample_type)
comparisons <- setdiff(levs, c("positive.control", "negative.control"))

# output results df
res_df <- lapply(comparisons, function(num) {
  r <- results(dds, contrast = c("sample_type", num, "positive.control"), alpha = 0.05)
  data.frame(gene           = rownames(r),
             baseMean       = r$baseMean,
             log2FoldChange = r$log2FoldChange,
             lfcSE          = r$lfcSE,
             stat           = r$stat,
             pvalue         = r$pvalue,
             padj           = r$padj,
             sample_type    = num,
             reference      = "positive.control")
})
names(res_df) <- paste0(comparisons, "_vs_positive.control")

all_res <- data.table::rbindlist(res_df) %>% as.data.frame()

write_tsv(all_res, "results/rnaseq/all-ffs-vs-lps-stimulation-alpha05.tsv.gz")

# save the RDS objects
saveRDS(dds, "dds_fitted.rds")
saveRDS(all_res, "all_res_vs_LPS.rds")

# summary on the results
alpha <- 0.05

res_summary <- all_res %>%
  group_by(sample_type, reference) %>%
  summarise(
    n_tested  = sum(baseMean > 0),
    up        = sum(padj < alpha & log2FoldChange > 0, na.rm = TRUE),
    down      = sum(padj < alpha & log2FoldChange < 0, na.rm = TRUE),
    outliers  = sum(baseMean > 0 & is.na(pvalue)),
    low_count = sum(!is.na(pvalue) & is.na(padj)),
    .groups = "drop"
  ) %>%
  mutate(total_DE = up + down) %>%
  arrange(desc(total_DE))

#################################################
# Comparisons of foods vs LPS+ for fold-change
# Filter to foods that show some amount of anti-inflammatory activity from the IFN bioactivity assay
# Specific food volcano plots vs LPS+
#################################################
# join results df with sample codes and gene metadata including gene name and biotype
all_res_df_info <- all_res %>% 
  left_join(gene_info) %>% 
  left_join(sample_codes, by="sample_type") %>% 
  select(gene, gene_name, gene_biotype, sample_type, sample_code, reference, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

# filter for samples with some anti-inflamm activity from IFN bioassay
base_filtered_codes <- c(bioactivity_ifn_filtered_codes, "POS", "PDTC", "ruxolitinib")

base_anti_inflm_foods_df <- all_res_df_info %>% 
  filter(sample_code %in% base_filtered_codes)

# top DE genes for foods vs LPS+ condition for foods that meet the IFN bioactivity cutoff
top_de_genes <- base_anti_inflm_foods_df %>% 
  filter(!is.na(padj), padj < 0.05) %>% 
  dplyr::count(gene, sort = TRUE) %>% 
  slice_head(n=50) %>% 
  pull(gene)

# build genes x comparisons matrix
lfc_mat <- base_anti_inflm_foods_df %>% 
  filter(gene %in% top_de_genes) %>% 
  select(gene, sample_type, log2FoldChange) %>% 
  tidyr::pivot_wider(names_from = sample_type, values_from = log2FoldChange) %>% 
  tibble::column_to_rownames("gene") %>% 
  as.matrix()

# ensemble Gene IDs to gene name
id_to_name <- base_anti_inflm_foods_df %>%
  dplyr::distinct(gene, gene_name) %>%
  { setNames(.$gene_name, .$gene) }

row_labels <- id_to_name[rownames(lfc_mat)]
row_labels[is.na(row_labels) | row_labels == ""] <- rownames(lfc_mat)[is.na(row_labels) | row_labels == ""]

# plot
lim    <- max(abs(lfc_mat), na.rm = TRUE)
lim    <- min(lim, 4)
breaks <- seq(-lim, lim, length.out = 101)
cols   <- colorRampPalette(c("blue", "white", "red"))(100)

pheatmap(pmax(pmin(lfc_mat, lim), -lim),
         color = cols, breaks = breaks,
         labels_row = row_labels,
         fontsize_col = 6,
         fontsize_row = 5)


#################################################
# Function for comparing subset list of foods to drug controls
# Filter to foods that show similar or lower activity IFN compared to drug ruxolitinib
#################################################

compare_to_drugs <- function(dds, foods,
                             drugs = c("Jak.Stat.inhibitor.ruxolitinib",
                                       "NFKB.inhibitor.PDTC")) {
  out <- list()
  for (ref in drugs) {
    dds_ref <- dds
    dds_ref$sample_type <- relevel(dds_ref$sample_type, ref = ref)
    dds_ref <- nbinomWaldTest(dds_ref)        # one-time step per drug
    
    for (f in foods) {
      r <- results(dds_ref, contrast = c("sample_type", f, ref))
      out[[paste0(f, "_vs_", ref)]] <- data.frame(
        gene           = rownames(r),
        baseMean       = r$baseMean,
        log2FoldChange = r$log2FoldChange,
        lfcSE          = r$lfcSE,
        stat           = r$stat,
        pvalue         = r$pvalue,
        padj           = r$padj,
        sample_type    = f,
        reference      = ref)
    }
    rm(dds_ref); gc()
  }
  data.table::rbindlist(out) %>% as.data.frame()
}


