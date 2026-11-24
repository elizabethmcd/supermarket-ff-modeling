library(tidyverse)
library(Hmisc)
library(readxl)

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

write_tsv(metadata_sheet, "metadata/cleaned-combined-metadata-sheet.tsv")

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

filtered_samples <- ifn_bioactivity_filtered_df %>% 
  filter(sample_code %in% filtered_bioactivity_sample_codes) %>% 
  select(sample_code, normalized_IFN) %>% 
  left_join(filtered_metadata_sheet_by_bioactivity_runs) %>% 
  group_by(sample_code, substrate_category) %>% 
  summarise(n=n()) %>% 
  filter(n > 2) %>% 
  pull(sample_code)

# filter bioactivity results to viability > 70%, sample codes that have more than 2 replicates, and those with metadata information in the metadata sheet
ifn_bioactivity_filtered_df <- bioactivity_df %>% 
  dplyr::mutate(Sample_type = tolower(Sample_type)) %>% 
  dplyr::rename(sample_code = Sample_number) %>% 
  filter(normalized_viability > 70) %>% 
  filter(sample_code %in% filtered_samples)

write_tsv(ifn_bioactivity_filtered_df, "results/bioactivity/filtered-ifn-bioactivity-df.tsv")

# plot bar graph of IFN bioactivity, selecting samples with more than 2 "replicates", color by substrate category
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
# test for significance to select samples significant different in IFN from LPS+
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

# perform a 2-sided t-test of each food vs the LPS+ positive control, then perform BH correction, 1-sided for less than LPS+
welch_results <- ifn_normalized_filtered_df %>% 
  filter(sample_type != "positive_control") %>% 
  group_by(sample_code, sample_type) %>% 
  summarise(
    n = dplyr::n(),
    mean_IFN = mean(normalized_IFN),
    sd_IFN = sd(normalized_IFN),
    test = list(t.test(normalized_IFN, lps_pos_values, alternative = "less")),
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

write_tsv(welch_results, "results/bioactivity/welch-t-test-signf-ifn.tsv")

# wilcoxon rank sum test with BH correction
# some samples have large outliers, compare Welch t-test to Wilcoxon test for significance from LPS-pos control
# 1-sided for less than LPS+
wilcox_results <- ifn_normalized_filtered_df %>% 
  filter(sample_type != "positive_control") %>% 
  group_by(sample_code, sample_type) %>% 
  summarise(
    n = dplyr::n(),
    median_IFN = median(normalized_IFN),
    IQR_IFN = IQR(normalized_IFN),
    test = list(wilcox.test(normalized_IFN, lps_pos_values, alternative="less")),
    .groups = "drop"
  ) %>% 
  mutate(
    w_stat = sapply(test, function(t) unname(t$statistic)),
    p_value = sapply(test, function(t) unname(t$p.value)),
    p_adj = p.adjust(p_value, method="BH"),
    significant = p_adj < 0.05
    ) %>% 
  select(-test) %>% 
  arrange(p_adj)

write_tsv(wilcox_results, "results/bioactivity/wilcox-signf-test-ifn.tsv")

sum(welch_results$p_adj < 0.05)
sum(wilcox_results$p_adj < 0.05)

# plot bar plot of mean normalized IFN with error bars, sorting the bars by highest to lowest mean IFN, and put * on the plot if they are significantly different from the LPS+ condition - make plots for each of the tests performed above
plot_df <- ifn_normalized_filtered_df %>% 
  filter(sample_code %in% filtered_samples) %>% 
  filter(sample_code != "FF0262") %>%  # drop the amazake sample because we don't have other paired data for it
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
                              levels = sample_order))

welch_signf_df <- welch_results %>% 
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

wilcox_signf_df <- wilcox_results %>% 
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

# order wilcoxon plot bars by median IFN
wilcox_sample_order <- plot_df %>% 
  group_by(sample_code) %>% 
  summarise(median_IFN = median(normalized_IFN), .groups = "drop") %>% 
  arrange(desc(median_IFN)) %>% 
  filter(sample_code != "FF0262") %>% 
  pull(sample_code)

wilcox_plot_df <- plot_df %>% 
  mutate(sample_code = factor(sample_code,
                              levels = wilcox_sample_order))

# mean signf bar plot with * significance from Welch's t-test
welch_ifn_mean_signf_bar_plot <- ggplot(plot_df, aes(x=sample_code, y=normalized_IFN,
       fill=substrate_category, color = substrate_category)) +
  stat_summary(fun=mean, geom="col", alpha=0.5, width=0.7, linewidth=0.8) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0.25, linewidth=0.6) +
  geom_jitter(width=0.15, size=1.5, alpha=0.8) +
  geom_hline(yintercept=100, linetype="dashed", color="grey40") +
  geom_text(data = welch_signf_df, aes(x=sample_code, y=y_pos, label=stars),
            inherit.aes = FALSE, size=5, vjust=0) +
  labs(x="Sample", y="Normalized IFN (% of LPS+)",
       fill="Substrate", color="Substrate",
       caption = "Mean ± SD. One-sided Welch's t-test vs LPS+, BH-adjusted: * q < 0.05, ** q < 0.01, *** q < 0.001") +
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

welch_ifn_mean_signf_bar_plot

ggsave("figures/ifn-mean-signf-bar-plot.png", welch_ifn_mean_signf_bar_plot, width=22, height=10, units=c("in"))

# same bar plot of median normalized IFN, with * from Wilcoxon test
wilcox_ifn_median_signf_bar_plot <- ggplot(wilcox_plot_df, aes(x=sample_code, y=normalized_IFN,
                                                      fill=substrate_category, color = substrate_category)) +
  stat_summary(fun=median, geom="col", alpha=0.5, width=0.7, linewidth=0.8) +
  stat_summary(fun=median,
               fun.min = ~quantile(.x, .25),
               fun.max = ~quantile(.x, .75),
               geom = "errorbar", width=0.25, linewidth=0.6) +
  geom_jitter(width=0.15, size=1.5, alpha=0.8) +
  geom_hline(yintercept=100, linetype="dashed", color="grey40") +
  geom_text(data = wilcox_signf_df, aes(x=sample_code, y=y_pos, label=stars),
            inherit.aes = FALSE, size=5, vjust=0) +
  labs(x="Sample", y="Normalized IFN (% of LPS+)",
       fill="Substrate", color="Substrate",
       caption = "One-sided Wilcoxon rank-sum test vs LPS+, BH-adjusted: * q < 0.05, ** q < 0.01, *** q < 0.001") +
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
    legend.text = element_text(size=12),
    plot.caption = element_text(size=12)
  )

wilcox_ifn_median_signf_bar_plot

ggsave("figures/ifn-wilox-median-signf-bar-plot.png", wilcox_ifn_median_signf_bar_plot, width=22, height=10, units=c("in"))

# of the filtered samples is IFN bioactivity correlated with viability or not
ifn_viability_plot <- ifn_normalized_filtered_df %>% 
  left_join(filtered_metadata_sheet_by_bioactivity_runs, by="sample_code") %>% 
  ggplot(aes(x=normalized_IFN, y=normalized_viability, color=substrate_category)) +
  geom_point() +
  theme_minimal()

ggsave("figures/ifn-viability-comp-plot.png", ifn_viability_plot, width=10, height=8, units=c("in"))

###############################
# Test if certain substrate categories are enriched in samples that have significantly lower IFN than the positive control
###############################

significance_by_substrate_category <- wilcox_results %>% 
  left_join(filtered_metadata_sheet_by_bioactivity_runs) %>% 
  distinct(sample_code, substrate_category, significant) %>% 
  filter(!sample_code %in% c("positive_control", "inhibitor2"))

set.seed(123)
fisher.test(table(significance_by_substrate_category$substrate_category, significance_by_substrate_category$significant), simulate.p.value=TRUE, B=1e5, alternative="less") # p=0.046 on the aggregate table, some evidence that samples with significantly lower IFN from the LPS+ differs across groups, so suppression isn't spread evenly across substrate category groups

enrichment <- purrr::map_dfr(unique(na.omit(significance_by_substrate_category$substrate_category)), function(cat) {
  in_cat <- significance_by_substrate_category$substrate_category == cat
  tab <- table(
    in_cat = factor(in_cat, levels = c(FALSE, TRUE)),
    lower  = factor(significance_by_substrate_category$significant, levels = c(FALSE, TRUE))
  )
  ft <- fisher.test(tab, alternative = "greater")
  tibble(
    category   = cat,
    n_samples  = sum(in_cat, na.rm = TRUE),
    n_lower    = sum(in_cat & significance_by_substrate_category$significant, na.rm = TRUE),
    pct_lower  = round(100 * n_lower / n_samples, 1),
    odds_ratio = unname(ft$estimate),
    p_value    = ft$p.value
  )
}) %>% 
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>% 
  arrange(p_adj) # no category comes out statistically significant
  
