library(tidyverse)
library(Hmisc)

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
  

ggsave("figures/ifn-bioactivity-bar-plot-filtered-samples.png", width=15, height=9, units=c("in"))


###############################
# Welch T-test for significance to select samples significant different in IFN from LPS+
# BH correction for multiple testing
###############################
# basic stats for the filtered samples for quick checks
bioactivity_stats <- ifn_bioactivity_filtered_df %>% 
  filter(sample_code %in% filtered_samples) %>% 
  dplyr::group_by(sample_code) %>% 
  summarise(mean_IFN = mean(normalized_IFN), sd_IFN = sd(normalized_IFN), median_IFN = median(normalized_IFN), minimum_IFN = min(normalized_IFN), maximum_IFN = max(normalized_IFN), n = dplyr::n(), .groups = "drop")

# set of filtered samples and controls
filtered_samples_and_controls <- c(filtered_samples, "positive_control", "inhibitor2")

# filtered df to select set of samples and controls
ifn_normalized_filtered_df <- bioactivity_df %>% 
  dplyr::mutate(sample_type = tolower(Sample_type)) %>% 
  dplyr::rename(sample_code = Sample_number) %>% 
  filter(normalized_viability > 70) %>% 
  select(sample_code, sample_type, normalized_IFN, normalized_viability) %>% 
  filter(sample_code %in% filtered_samples_and_controls)

# LPS+ values
lps_pos_values <- ifn_normalized_filtered_df %>% 
  filter(sample_type == "positive_control") %>% 
  pull(normalized_IFN)

# perform a Welch's 2-sided t-test of each food vs the LPS+ positive control, then perform BH correction
welch_results <- ifn_normalized_filtered_df %>% 
  filter(sample_type != "positive_control") %>% 
  group_by(sample_code, sample_type) %>% 
  summarise(
    n = dplyr::n(),
    mean_IFN = mean(normalized_IFN),
    sd_IFN = sd(normalized_IFN),
    test = list(t.test(normalized_IFN, lps_pos_values)),
    .groups = "drop"
  ) %>% 
  mutate(
    difference_from_LPS = mean_IFN - mean(lps_pos_values),
    CI_low = sapply(test, function(t) t$conf.int[1]),
    CI_high = sapply(test, function(t) t$conf.int[2]),
    t_stat = sapply(test, function(t) unname(t$statistic)),
    df = sapply(test, function(t) unname(t$parameter)),
    p_value = sapply(test, function(t) t$p.value),
    p_adj = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05
  ) %>% 
  select(-test) %>% 
  arrange(p_adj)
 
# plot bar plot of mean normalized IFN with error bars, sorting the bars by highest to lowest mean IFN, and put * on the plot if they are significantly different from the LPS+ condition 
plot_df <- ifn_normalized_filtered_df %>% 
  filter(sample_code %in% filtered_samples) %>% 
  select(sample_code, normalized_IFN) %>% 
  left_join(filtered_metadata_sheet_by_bioactivity_runs, by="sample_code")

sample_order <- plot_df %>% 
  group_by(sample_code) %>% 
  summarise(mean_IFN = mean(normalized_IFN), .groups = "drop") %>% 
  arrange(desc(mean_IFN)) %>% 
  filter(sample_code != "FF0262") %>% 
  pull(sample_code)

plot_df <- plot_df %>% 
  mutate(sample_code = factor(sample_code,
                              levels = sample_order)) %>% 
  filter(sample_code != "FF0262")

signf_df <- welch_results %>% 
  filter(sample_code %in% sample_order) %>% 
  mutate(
    stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    sample_code = factor(sample_code, levels = sample_order)
  ) %>% 
  left_join(
    plot_df %>% group_by(sample_code) %>% 
      summarise(y_pos = max(normalized_IFN) + 6, .groups = "drop"),
    by = "sample_code"
  )

ifn_mean_signf_bar_plot <- ggplot(plot_df, aes(x=sample_code, y=normalized_IFN,
       fill=substrate_category, color = substrate_category)) +
  stat_summary(fun=mean, geom="col", alpha=0.5, width=0.7, linewidth=0.8) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0.25, linewidth=0.6) +
  geom_jitter(width=0.15, size=1.5, alpha=0.8) +
  geom_hline(yintercept=100, linetype="dashed", color="grey40") +
  geom_text(data = signf_df, aes(x=sample_code, y=y_pos, label=stars),
            inherit.aes = FALSE, size=5, vjust=0) +
  labs(x="Sample", y="Normalized IFN (% of LPS+)",
       fill="Substrate", color="Substrate",
       caption = "Mean ± SD. Welch's t-test vs LPS+, BH-adjusted: * q < 0.05, ** q < 0.01, *** q < 0.001") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  guides(
    fill = guide_legend(ncol = 5, title.position = "top"),
    color = guide_legend(ncol = 5, title.position = "top")
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.5),
    axis.text.x = element_text(angle = 85, hjust = 1, size=12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(face="bold", size=14),
    axis.title.y = element_text(face= "bold", size=14),
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.95),
    legend.justification = c(1,1),
    legend.background = element_rect(fill="white", color="grey80"),
    legend.title = element_text(face="bold", size=14),
    legend.text = element_text(size=12)
  )

ggsave("figures/ifn-mean-signf-bar-plot.png", ifn_mean_signf_bar_plot, width=22, height=10, units=c("in"))

# of the filtered samples is IFN bioactivity correlated with viability or not
ifn_viability_plot <- ifn_normalized_filtered_df %>% 
  left_join(filtered_metadata_sheet_by_bioactivity_runs, by="sample_code") %>% 
  ggplot(aes(x=normalized_IFN, y=normalized_viability, color=substrate_category)) +
  geom_point() +
  theme_minimal()

ggsave("figures/ifn-viability-comp-plot.png", ifn_viability_plot, width=10, height=8, units=c("in"))
