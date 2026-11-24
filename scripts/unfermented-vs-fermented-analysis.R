library(tidyverse)
library(DESeq2)
library(ggrepel)
library(msigdbr)

###############################
# Unfermented vs fermented expression profiles
###############################

#################################################
# read in objects and metadata
#################################################

# metadata for unfermented vs fermented pairs
pairs_metadata <- read_tsv("metadata/rnaseq/2026-11-23-updated-unfermented-vs-fermented-pairs-sheet.tsv", col_names = c("fermented_sample_code", "fermented_sample_name", "unfermented_sample_code", "raw", "pair"), skip = 1)

# RDS objects
dds <- readRDS("results/rnaseq/dds_fitted.rds")
gene_info <- readRDS("results/rnaseq/gene_info.rds") %>% 
  dplyr::rename(gene_id = gene)
sample_codes <- readRDS("results/rnaseq/sample_codes.rds")

# join pair metadata sheet with sample codes to get correct sample_type formatting
sample_code_df <- as.data.frame(sample_codes)

pairs_sample_types_df <- left_join(pairs_metadata, sample_code_df, by=join_by("fermented_sample_code" == "sample_code")) %>% 
  filter(!is.na(sample_type)) %>% 
  dplyr::rename(fermented = sample_type) %>% 
  select(-fermented_sample_code, -fermented_sample_name) %>% 
  left_join(sample_code_df, by=join_by("unfermented_sample_code" == "sample_code")) %>% 
  select(-unfermented_sample_code, -raw) %>% 
  dplyr::rename(raw = sample_type) %>% 
  select(fermented, raw, pair)

# specify unferm vs ferm contrasts from the DESeq object
pairs <- pairs_sample_types_df %>% 
  mutate(across(c(fermented, raw), as.character),
         contrast_id = paste0(pair, "_", fermented, "_vs_", raw))


ferm_res <- pmap(pairs, function(fermented, raw, pair, contrast_id, ...){
  res <- results(dds,
                 contrast = c("sample_type", fermented, raw),
                 alpha=0.05)
  
  res_shr <- lfcShrink(dds,
                       contrast = c("sample_type", fermented, raw),
                       res = res,
                       type = "ashr")
  
  as.data.frame(res_shr) %>% 
    rownames_to_column("gene_id") %>% 
    mutate(pair = pair, fermented = fermented, raw = raw, contrast = contrast_id)
}) %>% set_names(pairs$contrast_id)

# table with gene annotations
ferm_res_df <- bind_rows(ferm_res) %>% 
  left_join(gene_info, by="gene_id") %>% 
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange >  1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "NS"
  ))

# DE genes per pair
de_genes_summary <- ferm_res_df %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>% 
  dplyr::count(contrast, direction = if_else(log2FoldChange > 0, "up_in_fermented", "down_in_fermented"))

# volcano plot and label the top genes
# remove the kefir comparison because doesn't have significant DE genes identified
top_genes <- ferm_res_df %>%
  filter(sig != "NS") %>% 
  filter(contrast != "milk_Straus.organic.kefir_vs_Straus.organic.milk") %>% 
  group_by(contrast) %>%
  slice_min(padj, n = 10, with_ties = FALSE)

ferm_vs_unferm_volcano_plot <- ferm_res_df %>% 
  filter(contrast != "milk_Straus.organic.kefir_vs_Straus.organic.milk") %>% 
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_text_repel(data = top_genes, aes(label = gene_name),
                  size = 3, color = "black", max.overlaps = 20) +
  scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_wrap(~ contrast) +
  labs(x = "log2 Fold Change of Fermented Food vs Raw Counterpart",
       y = "-log10 adjusted p-value",
       color = NULL) +
  theme_bw()

ggsave("figures/fermented-vs-unfermented-volcano-plots.png", ferm_vs_unferm_volcano_plot, units=c("in"), width=10, height=7)

# combine with pathway information
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  distinct()

hallmark_by_gene <- hallmark %>%
  group_by(ensembl_gene) %>%
  summarise(hallmark_pathways = paste(sort(unique(gs_name)), collapse = "; "),
            n_hallmark = n(),
            .groups = "drop")

ferm_res_annot_sig <- ferm_res_df %>%
  left_join(hallmark_by_gene, by = c("gene_id" = "ensembl_gene"))

# plot counts for certain genes
tnf <- gene_info$gene_id[gene_info$gene_id == "ENSG00000232810"]
keep <- unique(c(pairs$fermented, pairs$raw, "positive.control"))

tnf_counts_plot <- plotCounts(dds, gene = tnf, intgroup = "sample_type", returnData = TRUE) %>%
  filter(sample_type %in% keep) %>%
  ggplot(aes(reorder(sample_type, count, median), count)) +
  geom_point() + scale_y_log10() + coord_flip() +
  labs(x = NULL, y = "Normalized TNF counts")

ggsave("figures/tnf-counts-plot-unferm-vs-ferm.png", tnf_counts_plot, width=8, height=10, units=c("in"))
