library(maaslin3)
library(vegan)
library(tidyverse)
library(patchwork)

###############################
# Test associations of microbes/metabolites with metadata and IFN bioactivity
###############################

# general function for plotting PCoA
plot_pcoa <- function(d, meta, title) {
  p   <- wcmdscale(d, k = 2, eig = TRUE)
  eig <- p$eig[p$eig > 0]
  pct <- round(100 * eig / sum(eig), 1)
  df  <- as.data.frame(p$points) %>%
    setNames(c("PCoA1", "PCoA2")) %>%
    rownames_to_column("sample_code") %>%
    left_join(meta, by = "sample_code")
  plt <- ggplot(df, aes(PCoA1, PCoA2,
                        colour = substrate_category,
                        shape  = substrate_category)) +
    geom_point(size = 3, alpha = 0.9) +
    stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.4,
                 alpha = 0.6, show.legend = FALSE) +
    scale_colour_manual(values = pal_named) +
    scale_shape_manual(values = shapes_named) +
    labs(x = paste0("PCoA1 (", pct[1], "%)"),
         y = paste0("PCoA2 (", pct[2], "%)"),
         colour = "Substrate category", shape = "Substrate category",
         title = title) +
    theme_bw() +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 4)))
  list(plot = plt, pcoa = p)
}

# sample substrate category for uniform palettes 
cats   <- sort(unique(sample_metadata$substrate_category))
pal_named    <- setNames(pal14[seq_along(cats)], cats)
shapes_named <- setNames(shapes14[seq_along(cats)], cats)

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

# binary Jaccard sitance
jac_mg <- vegdist(abund > 0, method = "jaccard")

# plots with different distance metrics on same plot with shared legend
mg_bc  <- plot_pcoa(bc,     meta, "Metagenomics – Bray–Curtis")
mg_jac <- plot_pcoa(jac_mg, meta, "Metagenomics – Jaccard (presence/absence)")

combined_mg_plot <- mg_bc$plot / mg_jac$plot + plot_layout(guides = "collect")

ggsave("figures/metagenomics-pcoa-comp-distances.png", combined_mg_plot, width=10, height=12, units=c("in"))

###############################
# PCoA of metabolomics data using Euclidean distance calculations, filter out rarely detected metabolites
###############################
hilic_tidy <- read_tsv("raw_data/metabolomics/tidy_dfs/hilic_metabolomics_tidy.tsv") %>% 
  filter(!str_detect(sample_id, "^QC")) %>% 
  dplyr::rename(sample_code = sample_id) %>% 
  mutate(feature_id = paste(metabolite_name, Adduct, Polarity, sep="_"))

# pivot long to wide for matrix
hilic_wide <- hilic_tidy %>% 
  select(sample_code, feature_id, metabolite_intensity) %>% 
  pivot_wider(names_from = feature_id,
              values_from = metabolite_intensity) %>% 
  column_to_rownames("sample_code") %>% 
  as.matrix()
hilic_wide[is.na(hilic_wide)] <- 0

metab_met <- sample_metadata %>%
  filter(sample_code %in% rownames(hilic_wide),
         !is.na(substrate_category))
metab_met <- metab_met[match(intersect(rownames(hilic_wide), metab_met$sample_code),
                             metab_met$sample_code), ]
keep <- metab_met$sample_code

# filter to only include feature ID (metabolites) that are present in at least 10% of samples
hilic_wide <- hilic_wide[, colSums(hilic_wide > 0) >=.10 * nrow(hilic_wide)]

# normalize to total signal per sample
normalized <- hilic_wide / rowSums(hilic_wide) * median(rowSums(hilic_wide))

imp <- apply(normalized, 2, function(x) {
  x[x == 0] <- min(x[x > 0]) / 2
  x
})

# log10 transform plus pareto scaling
X <- log10(imp)
X <- sweep(scale(X, center = TRUE, scale = FALSE), 2,
           sqrt(apply(X, 2, sd)), "/")

# euclidean distance and PCoA
d_euc <- dist(X, method = "euclidean")
pcoa  <- wcmdscale(d_euc, k = 2, eig = TRUE)

eig <- pcoa$eig[pcoa$eig > 0]
pct <- round(100 * eig / sum(eig), 1)

scores_df <- as.data.frame(pcoa$points) %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column("sample_code") %>%
  left_join(sample_metadata, by = "sample_code")

d_euc_k <- as.dist(as.matrix(d_euc)[keep, keep])


# binary jaccard on detected features, using the raw matrix
jac_met <- vegdist(hilic_wide[keep, ] > 0, method = "jaccard")

# plot PCoA of HILIC metabolites with same palette as BC PCoA metagenomics plot
pal14 <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
           "#42D4F4", "#F032E6", "#BFEF45", "#469990", "#9A6324",
           "#800000", "#808000", "#000075", "#FFD700")

shapes14 <- rep(c(16, 17, 15, 18, 8, 4, 3), length.out = 14)

ggplot(scores_df, aes(PCoA1, PCoA2, colour = substrate_category, shape = substrate_category)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.4,
               alpha = 0.6, show.legend = FALSE) +
  scale_colour_manual(values = pal14) +
  scale_shape_manual(values = shapes14) +
  labs(x = paste0("PCoA1 (", pct[1], "%)"),
       y = paste0("PCoA2 (", pct[2], "%)"),
       colour = "Substrate category", shape = "Substrate category",
       title = "HILIC metabolome – Euclidean (log10, Pareto-scaled)") +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 4)))

# plot euclidean vs jaccard PCoA
metab_euc <- plot_pcoa(d_euc_k, metab_met, "Metabolomics – Euclidean (log10, Pareto)")
metab_jac <- plot_pcoa(jac_met, metab_met, "Metabolomics – Jaccard (presence/absence)")

combined_metabolomics_pcoa_plots <- metab_euc$plot / metab_jac$plot + plot_layout(guides = "collect")

ggsave("figures/metabolomics-pcoa-comp-distances.png", combined_metabolomics_pcoa_plots, width=10, height=12, units=c("in"))

adonis2(jac_met ~ substrate_category, data = metab_met, permutations = 999)
adonis2(jac_mg ~ substrate_category, data = meta, permutations = 999)
anova(betadisper(jac_met, metab_met$substrate_category))
