library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(genefilter)

#################################################
# Comparisons of foods vs LPS+ for differentially expressed genes
#################################################

#################################################
# read in objects
#################################################

dds <- readRDS("results/rnaseq/dds_fitted.rds")
gene_info <- readRDS("results/rnaseq/gene_info.rds")
sample_codes <- readRDS("results/rnaseq/sample_codes.rds")

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

# write to RDS object
saveRDS(all_res, "results/rnaseq/all_res_vs_lps.rds")

# filter to padj < 0.05 and write out the DEGs
filtered_degs <- all_res %>% 
  filter(padj < 0.05)

write_tsv(filtered_degs, "results/rnaseq/filtered-05-degs.tsv.gz")

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

write_tsv(res_summary, "results/rnaseq/degs-food-summaries.tsv")

#################################################
# Specific food volcano plots vs LPS+
#################################################
# join results df with sample codes and gene metadata including gene name and biotype
all_res_df_info <- all_res %>% 
  left_join(gene_info) %>% 
  left_join(sample_codes, by="sample_type") %>% 
  select(gene, gene_name, gene_biotype, sample_type, sample_code, reference, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

write_tsv(all_res_df_info, "results/rnaseq/all-ffs-vs-lps-stimulation-alpha05.tsv.gz")

# top DE genes for foods vs LPS+ condition for foods that meet the IFN bioactivity cutoff
top_de_genes <- all_res_df_info %>% 
  filter(!is.na(padj), padj < 0.05) %>% 
  dplyr::count(gene, sort = TRUE) %>% 
  slice_head(n=50) %>% 
  pull(gene)

# build genes x comparisons matrix
lfc_mat <- all_res_df_info %>% 
  filter(gene %in% top_de_genes) %>% 
  select(gene, sample_type, log2FoldChange) %>% 
  tidyr::pivot_wider(names_from = sample_type, values_from = log2FoldChange) %>% 
  tibble::column_to_rownames("gene") %>% 
  as.matrix()

# ensemble Gene IDs to gene name
id_to_name <- all_res_df_info %>%
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
