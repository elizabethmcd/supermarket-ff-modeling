library(tidyverse)
library(ggrepel)

#######################################
# RNA-seq QC checks
#######################################

# sample map and remove singleton samples
sample_map <- read_tsv("metadata/rnaseq/2026-09-25-rnaseq-sample-map.tsv", col_names = TRUE) %>% 
  set_names(c("sample_code", "batch", "sample_id", "lps_addition", "pre_treatment", "sample_type", "main_ingredient", "food_type", "replicate")) %>% 
  mutate(sample = paste0(sample_id, "_count")) %>% 
  select(sample, batch, sample_code, sample_type, replicate)

singleton_samples <- sample_map %>% 
  distinct() %>% 
  group_by(sample_type) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  mutate(sample = gsub("_count", "", sample)) %>% 
  pull(sample)

# drop additional samples - samples not from this experiment, the chemical control comparisons
singleton_samples <- c(singleton_samples, "44R9PF_17", "44R9PF_18", "44R9PF_19", "44R9PF_20", "44R9PF_21")

filtered_sample_map <- sample_map %>% 
  mutate(sample = gsub("_count", "", sample)) %>% 
  filter(!sample %in% singleton_samples)

# read in mapping stats and gene biotype counts CSVs
mapping_stats_dir <- "raw_data/rnaseq/mapping_stats_qc"

mapping_files <- list.files(mapping_stats_dir, pattern = "-mapping-stats-reads.csv", full.names=TRUE)

# combine all mapping stats CSVs together, filter out the singleton samples
combined_mapping_stats_df <- mapping_files %>% 
  set_names(~ str_remove(basename(.x), "-mapping-stats-reads.csv")) %>% 
  map(\(f) read_csv(f, show_col_types = FALSE) %>% dplyr::rename(sample = 1)) %>%
  list_rbind(names_to = "batch") %>% 
  mutate(
    sample = case_when(
      str_starts(sample, "sample") ~ str_replace(sample, "^sample", batch),
      str_starts(sample, "FF") ~ paste0(batch, "_", row_number()),
      .default = paste0(batch, "_", row_number())
    ),
    .by = batch
  ) %>% 
  relocate(batch, sample) %>% 
  filter(!sample %in% singleton_samples)

setdiff(combined_mapping_stats_df$sample, filtered_sample_map$sample)
setdiff(filtered_sample_map$sample, combined_mapping_stats_df$sample)

# combine all gene biotype count CSVs together, filter out singleton samples
gene_counts_dir <- "raw_data/rnaseq/gene_biotype_counts"

gene_counts_files <- list.files(gene_counts_dir, pattern = "-summary.csv", full.names=TRUE)

combined_gene_counts_df <- gene_counts_files %>% 
  set_names(~ str_remove(basename(.x), "-gene-biotype-5plus_reads-summary.csv")) %>% 
  map(\(f) read_csv(f, show_col_types = FALSE) %>% dplyr::rename(sample = 1)) %>%
  list_rbind(names_to = "batch") %>% 
  mutate(
    sample = case_when(
      str_starts(sample, "sample") ~ str_replace(sample, "^sample", batch),
      str_starts(sample, "FF") ~ paste0(batch, "_", row_number()),
      .default = paste0(batch, "_", row_number())
    ),
    .by = batch
  ) %>% 
  relocate(batch, sample) %>% 
  filter(!sample %in% singleton_samples)

# join sample map info, mapping stats, and gene biotype counts for the filtered, non-singleton samples
combined_qc_stats_df <- left_join(filtered_sample_map, combined_mapping_stats_df) %>% 
  left_join(combined_gene_counts_df) %>% 
  dplyr::rename(uniquely_mapped = `Uniquely Mapped`, multi_mapped = `Multi-mapped`, unmapped = `Unmapped`) %>% 
  mutate(
    total_reads = uniquely_mapped + multi_mapped + unmapped, 
    percent_mapped = (uniquely_mapped + multi_mapped) / total_reads * 100
    )

# plot percent mapped, ordered by sample and colored by batch
plot_df <- combined_qc_stats_df %>% 
  arrange(sample_type, sample, batch) %>% 
  mutate(sample = fct_inorder(sample))

mapping_plot <- ggplot(plot_df, aes(x=sample, y=percent_mapped, fill=batch)) +
  geom_col() +
  labs(x="Sample",
       y="Percent Reads Mapped (%)",
       fill="Batch",
       title = "Percent Reads Mapped to Reference Genome per Sample") +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = c(0,0))

# plot count of protein coding genes
gene_counts_plot <- ggplot(plot_df, aes(x=sample, y=protein_coding, fill=batch)) +
  geom_col() +
  labs(x="Sample",
       y="Count of Protein Coding Genes",
       fill="Batch",
       title = "Counts of Protein Coding Genes per Sample") +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = c(0,0))

# plot mapping rate vs count of protein coding genes
pc_mapping_comps_plot <- ggplot(plot_df, aes(x=protein_coding, y=percent_mapped, color=batch)) +
  geom_point(size=3) +
  labs(
    x="Count of Protein Coding Genes in a Sample",
    y="Percentage of Reads Mapped to Reference Genomes in a Sample",
    color="Batch",
    title="Counts of Protein Coding Genes vs \n Percentage of Reads Mapped to Reference Genomes in a Sample"
  ) +
  scale_color_brewer(palette="Paired") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand=c(0,0))

# save plots
ggsave("figures/mapping-stats-plot.png", mapping_plot, width=18, height=8, units=c("in"))
ggsave("figures/protein-coding-gene-counts-plot.png", gene_counts_plot, width=18, height=8, units=c("in"))
ggsave("figures/pc-mapping-comps-plot.png", pc_mapping_comps_plot, width=18, height=8, units=c("in"))

#######################################
# RNA-seq QC checks
#######################################
# check the M6PXYP run which had lower mapping % for a lot of samples, but the coverage across protein coding genes might look ok just some magnitude lower than the other samples
# compare to a run where all samples look fine - 4M38F8
counts <- read_tsv("raw_data/rnaseq/2026-09-25-combined-expression-matrix.tsv") %>% 
  select(gene_id, gene_name, gene_biotype, ends_with("_cpm")) %>% 
  pivot_longer(
    cols=
      !c(gene_id, gene_name, gene_biotype),
    names_to = "sample",
    values_to = "cpm"
  ) %>% 
  mutate(sample = gsub("_cpm", "", sample),
         batch = str_sub(sample, 1, 6))

# share of counts in each gene_biotype
biotype_df <- counts %>% 
  group_by(sample, batch, gene_biotype) %>% 
  dplyr::summarise(frac = sum(cpm, na.rm = TRUE)/ 1e6, .groups = "drop")

biotype_df %>% 
  group_by(sample) %>% 
  dplyr::summarise(total = sum(frac)) %>% 
  print(n=57)

# keep only biotypes above 0.5% mapping rate
keep_types <- biotype_df %>% 
  group_by(gene_biotype) %>% 
  filter(max(frac) >= 0.005) %>% 
  pull(gene_biotype) %>% 
  unique()

biotype_plot_df <- biotype_df %>% 
  mutate(gene_biotype = if_else(gene_biotype %in% keep_types, gene_biotype, "other")) %>% 
  group_by(sample, batch, gene_biotype) %>% 
  dplyr::summarise(frac = sum(frac), .groups = "drop")

# stacked bar plot for gene biotypes per sample
biotype_bar_plot <- ggplot(biotype_plot_df, aes(x=sample, y=frac, fill=gene_biotype)) +
  geom_col() +
  facet_grid(~ batch, scales = "free_x", space = "free_x") + 
  labs(x = NULL, y= "Fraction of CPM", fill = "Biotype") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + 
  theme(axis.text.x = element_text(
    angle = 90, vjust = 0.5, hjust=1
  )) +
  theme(legend.position = "bottom")

# for each biotype, make a boxlot of fraction of counts going to that biotype - facet by biotype and colored by batch
biotype_box_plot <- ggplot(biotype_plot_df, aes(x=batch, y=frac, color = batch)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.5, size = 2) +
  facet_wrap(~ gene_biotype, scales = "free_y") +
  labs(x=NULL, y="Fraction of CPM") +
  theme_bw() +
  theme(legend.position = "none")

# save figures
ggsave("figures/gene-biotype-mapping-counts-frac-bar-plot.png", biotype_bar_plot, width=16, height=8, units=c("in"))
ggsave("figures/gene-biotype-mapping-counts-frac-box-plot.png", biotype_box_plot, width=14, height=8, units=c("in"))

#######################################
# Filter out samples with low protein-coding counts, below 10,000, remove the PDTC samples, and remove the positive control acids because we don't have bioactivity data for them
#######################################

# remove samples with less than 10,000 protein coding gene count, remove PDTC samples, remove acids, and remove singletons if remaining
qced_sample_map <- combined_qc_stats_df %>% 
  filter(protein_coding > 10000) %>% 
  filter(!sample_type %in% c("NFKB inhibitor PDTC", "0.0125% acetic acid", "0.025% lactic acid", "4-Hydroxyphenyllactic acid", "Indole-3-lactic acid", "Phenyllactic acid", "amazake")) %>% 
  distinct() %>% 
  group_by(sample_type) %>% 
  filter(n() > 1) %>% 
  ungroup()

write_tsv(qced_sample_map, "metadata/rnaseq/2026-10-26-qced-sample-map.tsv")  
