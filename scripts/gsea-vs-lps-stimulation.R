library(tidyverse)
library(fgsea)
library(msigdbr)
library(pheatmap)

###############################
# GSEA pathway analysis of fermented food extracts vs LPS+ condition
###############################

# read in saved RDS objects
all_res      <- readRDS("all_res_vs_LPS.rds")
gene_info    <- readRDS("gene_info.rds")
sample_codes <- readRDS("sample_codes.rds")

###############################
# Pathway analysis: fgsea on MSigDB Hallmark gene sets
# Ranked by the Wald stat from each fermented food vs the LPS+ comparison
###############################

# msigdbr sets of genes to match the ensemble IDs
hallmark_genes <- msigdbr(species = "Homo sapiens", collection = "H")
pathways_h <- split(hallmark_genes$db_ensembl_gene, hallmark_genes$gs_name)

# make the Wald stat ranks 
# Ranked gene vector for one comparison
make_ranks <- function(df) {
  df <- df[!is.na(df$stat), ]
  sort(setNames(df$stat, df$gene), decreasing = TRUE)
}

# fgsea for every fermented food comparison to LPS+
set.seed(42)
gsea_res <- lapply(split(all_res, all_res$sample_type), function(df) {
  fgsea(pathways = pathways_h,
        stats    = make_ranks(df),
        minSize  = 15,
        maxSize  = 500,
        eps      = 0)
})

gsea_all <- data.table::rbindlist(gsea_res, idcol = "sample_type") %>%
  as.data.frame()

# Leading-edge genes as gene names
id2name <- setNames(gene_info$gene_name, gene_info$gene)
gsea_all$leadingEdge <- sapply(gsea_all$leadingEdge, function(x) {
  nm <- id2name[x]
  nm[is.na(nm) | nm == ""] <- x[is.na(nm) | nm == ""]
  paste(nm, collapse = ";")
})

gsea_all <- gsea_all %>%
  left_join(sample_codes, by = "sample_type") %>%
  dplyr::select(sample_type, sample_code, pathway, NES, pval, padj, size, leadingEdge)

write_tsv(gsea_all, "results/rnaseq/fgsea-hallmark-ffs-vs-lps.tsv")
saveRDS(gsea_all, "gsea_hallmark_vs_LPS.rds")

# make heatmap of NES comparing pathways to fermented foods
# first filter down to the fermented food samples that have some lowered IFN bioactivity that fits the filtered viability > 70%
drugs <- c("Jak.Stat.inhibitor.ruxolitinib", "NFKB.inhibitor.PDTC")

filtered_ifn_samples <- sample_codes %>%
  dplyr::filter(sample_code %in% bioactivity_ifn_filtered_codes) %>%
  dplyr::pull(sample_type) %>%
  unique() %>%
  as.character()

gsea_sub <- gsea_all %>% dplyr::filter(sample_type %in% c(filtered_ifn_samples, drugs))

# matrix of NES, pathway vs sample
nes_mat <- gsea_sub %>%
  dplyr::select(pathway, sample_type, NES) %>%
  tidyr::pivot_wider(names_from = sample_type, values_from = NES) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# put an asterisk where padj < 0.05 for the pathway
sig_mat <- gsea_sub %>%
  dplyr::mutate(star = ifelse(!is.na(padj) & padj < 0.05, "*", "")) %>%
  dplyr::select(pathway, sample_type, star) %>%
  tidyr::pivot_wider(names_from = sample_type, values_from = star) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()
sig_mat <- sig_mat[rownames(nes_mat), colnames(nes_mat)]

# cluster foods together, put the drugs at the end separately
food_cols  <- setdiff(colnames(nes_mat), drugs)
food_order <- food_cols[hclust(dist(t(nes_mat[, food_cols])))$order]
col_order  <- c(food_order, intersect(drugs, colnames(nes_mat)))

nes_lim <- max(abs(nes_mat), na.rm = TRUE)

# heatmap
pheatmap(nes_mat[, col_order],
         color           = colorRampPalette(c("blue", "white", "red"))(100),
         breaks          = seq(-nes_lim, nes_lim, length.out = 101),
         cluster_cols    = FALSE,
         gaps_col        = length(food_order),
         display_numbers = sig_mat[, col_order],
         labels_row      = sub("^HALLMARK_", "", rownames(nes_mat)),
         fontsize_row    = 6,
         fontsize_col    = 6,
         main            = "Hallmark GSEA: NES vs LPS+ (* padj < 0.05)")

