library(maaslin3)
library(vegan)
library(tidyverse)

###############################
# Test associations of microbes/metabolites with metadata and IFN bioactivity
###############################

###############################
# PCoA of metagenomics relative abundance colored by substrate category
###############################

# read in sylph profiles and metadata
sylph_profiles <- read_tsv("results/metagenomics/2026-01-30-seqcoast-supermarket-sweep-profiles.tsv") %>% 
  dplyr::rename(relative_abundance = Sequence_abundance) %>% 
  select(sample_code, species, relative_abundance)

sample_metadata <- read_tsv("metadata/cleaned-combined-metadata-sheet.tsv") %>% 
  select(sample_code, substrate_category)

# pivot sylph relative abundance wide to long with species as columns
abund <- sylph_profiles %>% 
  distinct(sample_code, species, .keep_all = TRUE) %>% 
  pivot_wider(names_from = species,
              values_from = relative_abundance,
              values_fill = 0) %>% 
  column_to_rownames("sample_code") %>% 
  as.matrix()

# match metadata in order of the matrix rownames
meta <- sample_metadata %>% 
  filter(sample_code %in% rownames(abund))
abund <- abund[meta$sample_code, , drop=FALSE]

# Bray-Curtis distance for PCoA
bc <- vegdist(abund / 100, method = "bray")
pcoa <- wcmdscale(bc, k=2, eig=TRUE)
eig <- pcoa$eig[pcoa$eig > 0]
pct <- round(100 * eig / sum(eig), 1)

scores_df <- as.data.frame(pcoa$points) %>% 
  setNames(c("PCoA1", "PCoA2")) %>% 
  rownames_to_column("sample_code") %>% 
  left_join(meta, by="sample_code")

# plot with substrate category coloring the points
pal14 <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
           "#42D4F4", "#F032E6", "#BFEF45", "#469990", "#9A6324",
           "#800000", "#808000", "#000075", "#FFD700")

shapes14 <- rep(c(16, 17, 15, 18, 8, 4, 3), length.out = 14)

ggplot(scores_df, aes(PCoA1, PCoA2,
                      colour = substrate_category,
                      shape  = substrate_category)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.4,
               alpha = 0.6, show.legend = FALSE) +
  scale_colour_manual(values = pal14) +
  scale_shape_manual(values = shapes14) +
  labs(x = paste0("PCoA1 (", pct[1], "%)"),
       y = paste0("PCoA2 (", pct[2], "%)"),
       colour = "Substrate category",
       shape  = "Substrate category",   # same title merges the two legends
       title = "Bray–Curtis PCoA of sylph profiles") +
  theme_bw() +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 4)))

adonis2(bc ~ substrate_category, data = meta, permutations = 999)
anova(betadisper(bc, meta$substrate_category))
