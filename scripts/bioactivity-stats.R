library(tidyverse)

###############################
# Select samples with IFN bioactivity at certain cutoffs
# Perform statistical tests comparing IFN bioactivity of each fermented food extract to the LPS+ condition
###############################
# bioactivity excel file
bioactivity_df <- read_excel("raw_data/bioactivity/Fermented_food_data_mastersheet.xlsx", sheet = "full_data")

# metadata of paired sampling for bioactivity, metagenomics, metabolomics, RNAseq
metadata_sheet <- read_tsv("metadata/Combined master sheet updated 042126-8ed4f7f4.tsv", col_names = c("sample_code", "sample", "substrate_category", "food_type", "sample_state", "sample_source", "metabolomics_collection", "metagenomics_collection", "bioassay_runs", "rnaseq_runs", "ingredients", "producer_1", "producer_country", "main_ingredient", "dominant_fermentation_chemistry_1", "sample_notes", "all_ingredients", "main_ingredient_organic_status", "producer_2", "dominant_fermentation_chemistry_2", "food_ontology_high", "food_ontology_mid", "food_ontology_specific", "food_ontology_parent_label", "food_ontology_gradparent_label", "food_ontology_lineage", "food_ontology_specific_canonical_id"), skip = 1) %>% 
  select(-producer_2, dominant_fermentation_chemistry_2) %>% 
  dplyr::rename(producer = producer_1, dominant_fermentation_chemistry = dominant_fermentation_chemistry_1)

# filter sample codes to those meeting these conditions:
# metabolomics - yes
# metagenomics - yes
# bioactivity runs - minimum 2
# rnaseq runs -> not for modeling purposes, just extra information/mechanism for specific samples so won't filter by this for now, otherwise would filter by at least 2 runs

filtered_metadata_sheet <- metadata_sheet %>% 
  filter(
    metabolomics_collection == "yes",
    metagenomics_collection == "yes", 
    bioassay_runs >= 2
  ) # 50 samples that have all data types

filtered_metadata_sheet_by_bioactivity_runs <- metadata_sheet %>% 
  filter(bioassay_runs >= 2) # 80 samples that have at least two bioassay runs

filtered_bioactivity_sample_codes <- filtered_metadata_sheet_by_bioactivity_runs %>% 
  pull(sample_code)

# filter bioactivity results to viability > 70% and normalized IFN < 100
ifn_bioactivity_filtered_df <- bioactivity_df %>% 
  dplyr::mutate(Sample_type = tolower(Sample_type)) %>% 
  dplyr::rename(sample_code = Sample_number) %>% 
  filter(normalized_viability > 70)

write_tsv(ifn_bioactivity_filtered_df, "results/bioactivity/filtered-ifn-bioactivity-df.tsv")

# plot bar graph of IFN bioactivity, selecting samples with more than 2 "replicates", color by substrate category
filtered_samples <- ifn_bioactivity_filtered_df %>% 
  filter(sample_code %in% filtered_bioactivity_sample_codes) %>% 
  select(sample_code, normalized_IFN) %>% 
  left_join(filtered_metadata_sheet_by_bioactivity_runs) %>% 
  group_by(sample_code, substrate_category) %>% 
  summarize(median_IFN = signif(median(normalized_IFN), digits=3), n=n()) %>% 
  filter(n > 2) %>% 
  pull(sample_code)

ifn_bar_plot <- ifn_bioactivity_filtered_df %>% 
  filter(sample_code %in% filtered_samples) %>% 
  select(sample_code, normalized_IFN) %>% 
  left_join(filtered_metadata_sheet_by_bioactivity_runs) %>% 
  ggplot(aes(x = reorder(sample_code, -normalized_IFN, FUN = median),
             y = normalized_IFN, fill = substrate_category, color = substrate_category)) +
  stat_summary(fun = median, geom = "col", alpha = 0.5, width = 0.7, linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  labs(x = "sample_code") +
  theme_minimal() +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.5),
    axis.text.x = element_text(angle = 85, hjust = 1))
  

ggsave("figures/ifn-bioactivity-bar-plot-filtered-samples.png", width=13, height=9, units=c("in"))


# samples that show some anti-inflammatory activity by IFN
bioactivity_ifn_filtered_codes <- ifn_bioactivity_filtered_df %>% 
  select(sample_code) %>% 
  unique() %>% 
  pull(sample_code)

# median ruxolitinib IFN activity
ruxolitinib_median_IFN <- ifn_bioactivity_filtered_df %>% 
  filter(sample_code == "inhibitor2") %>% 
  summarise(median(normalized_IFN)) %>% 
  pull(`median(normalized_IFN)`)

# samples that show at or below IFN compared to the median ruxolitinib IFN activity
bioactivity_ifn_filtered_ruxo_level_codes <- ifn_bioactivity_filtered_df %>% 
  group_by(sample_code) %>% 
  summarise(median = median(normalized_IFN)) %>% 
  filter(median <= ruxolitinib_median_IFN) %>% 
  pull(sample_code)

###############################
# Welch T-test for significance to select samples significant different in IFN from LPS+
# BH correction for multiple testing
###############################

bioactivity_stats <- bioactivity_df %>% 
  dplyr::mutate(sample_type = tolower(Sample_type)) %>% 
  dplyr::rename(sample_code = Sample_number) %>% 
  dplyr::group_by(sample_code) %>% 
  summarise(mean_IFN = mean(normalized_IFN), sd_IFN = sd(normalized_IFN), median_IFN = median(normalized_IFN), minimum_IFN = min(normalized_IFN), maximum_IFN = max(normalized_IFN), n = dplyr::n(), .groups = "drop")
