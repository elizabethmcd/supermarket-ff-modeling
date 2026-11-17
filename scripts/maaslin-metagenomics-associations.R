library(maaslin3)
library(tidyverse)

###############################
# Use MaAsLin3 to test associations of microbial species prevalence/abundance from metagenomics data with IFN bioactivity levels
###############################

# prep input data to input_data and input_metadata in format required by Maaslin3

# sample metadata
sample_metadata <- read_tsv("metadata/cleaned-combined-metadata-sheet.tsv")

# bioactivity results
filtered_bioactivity_results <- read_tsv("results/bioactivity/filtered-ifn-bioactivity-df.tsv")

# merge sample metadata with median bioactivity
sample_metadata_bioactivity_table <- left_join(filtered_bioactivity_results, sample_metadata) %>% 
  select(sample_code, sample, substrate_category, normalized_IFN) %>% 
  group_by(sample_code, sample, substrate_category) %>% 
  summarise(median_IFN = median(normalized_IFN)) %>% 
  column_to_rownames("sample_code") %>% 
  mutate(substrate_category = as.factor(substrate_category))

# sylph metagenomic profiling
sylph_results <- read_tsv("results/metagenomics/2026-01-30-seqcoast-supermarket-sweep-profiles.tsv") %>% 
  select(sample_code, species, Sequence_abundance) %>% 
  pivot_wider(
    id_cols = sample_code,
    names_from = species,
    values_from = Sequence_abundance, 
    values_fill = 0, 
    values_fn = sum
  ) %>% 
  column_to_rownames("sample_code")

# maaslin3 model fitting
set.seed(1)
fit_out <- maaslin3(input_data = sylph_results,
                    input_metadata = sample_metadata_bioactivity_table,
                    output = "results/metagenomics",
                    formula = "~ substrate_category",
                    normalization = "TSS",
                    transform = "LOG",
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1
                    )

