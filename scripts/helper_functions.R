library(tidyverse)

# compare each sample to the drugs
compare_to_drugs <- function(dds, foods,
                             drugs = c("Jak.Stat.inhibitor.ruxolitinib",
                                       "NFKB.inhibitor.PDTC"),
                             alpha = 0.05) {
  out <- list()
  for (ref in drugs) {
    for (f in foods) {
      r <- results(dds, contrast = c("sample_type", f, ref), alpha = alpha)
      out[[paste0(f, "_vs_", ref)]] <- data.frame(
        gene = rownames(r), baseMean = r$baseMean,
        log2FoldChange = r$log2FoldChange, lfcSE = r$lfcSE,
        stat = r$stat, pvalue = r$pvalue, padj = r$padj,
        sample_type = f, reference = ref)
    }
  }
  data.table::rbindlist(out) %>% as.data.frame()
}