library(tidyverse)

###############################
# Select samples with IFN bioactivity at certain cutoffs
# Perform statistical tests comparing IFN bioactivity of each fermented food extract to the LPS+ condition
###############################
bioactivity_df <- read_excel("raw_data/bioactivity/Fermented_food_data_mastersheet.xlsx", sheet = "full_data")

ifn_bioactivity_filtered_df <- bioactivity_df %>% 
  dplyr::mutate(Sample_type = tolower(Sample_type)) %>% 
  dplyr::rename(sample_code = Sample_number) %>% 
  filter(normalized_viability > 70) %>% 
  filter(normalized_IFN <=100)

write_tsv(ifn_bioactivity_filtered_df, "results/bioactivity/filtered-ifn-bioactivity-df.tsv")

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
