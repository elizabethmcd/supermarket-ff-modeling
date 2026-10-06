library(tidyverse)

#######################################
# RNA-seq QC checks
#######################################

# sample map and remove singleton samples
sample_map <- read_tsv("metadata/rnaseq/2026-09-25-rnaseq-sample-map.tsv", col_names = TRUE) %>% 
  set_names(c("sample_code", "batch", "sample_id", "lps_addition", "pre_treatment", "sample_type", "main_ingredient", "food_type", "replicate")) %>% 
  mutate(sample = paste0(sample_id, "_count")) %>% 
  select(sample, batch, sample_type, replicate)

singleton_samples <- sample_map %>% 
  distinct() %>% 
  group_by(sample_type) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  mutate(sample = gsub("_count", "", sample)) %>% 
  pull(sample)

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

# calculate median mapping rate and protein coding gene count
median_mapping_rate <- combined_qc_stats_df %>% 
  summarise(median(percent_mapped))

median_protein_coding_count <- combined_qc_stats_df %>% 
  summarise(median(protein_coding))

# remove samples with below 97% mapping rate, and then remove samples that are then singletons
qced_sample_map <- combined_qc_stats_df %>% 
  filter(percent_mapped > 97) %>% 
  distinct() %>% 
  group_by(sample_type) %>% 
  filter(n() > 1) %>% 
  ungroup()

write_tsv(qced_sample_map, "metadata/rnaseq/2026-10-06-qced-sample-map.tsv")  
