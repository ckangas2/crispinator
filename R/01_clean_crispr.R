library(tidyverse)

# R/01_clean_crispr.R

# Entry Point 1: Standard Lab MAGeCK
clean_mageck_lab <- function(file_path) {
  read_tsv(file_path, show_col_types = FALSE) %>%
    # REVERT: We explicitly rename to the pipeline standard (Gene, LFC, FDR)
    # We use 'any_of' to be safe if names vary slightly
    dplyr::select(Gene = any_of(c("id", "Gene", "OFFICIAL_SYMBOL")), 
                  LFC = any_of(c("neg|lfc", "LFC", "lfc")), 
                  FDR = any_of(c("neg|fdr", "FDR", "fdr"))) %>%
    mutate(
      Category = case_when(
        FDR < 0.05 & LFC < -1 ~ "Depleted (Essential)",
        FDR < 0.05 & LFC > 1  ~ "Enriched (Resistant)",
        TRUE ~ "Non-Significant"
      ),
      Label = if_else(FDR < 0.01 & abs(LFC) > 1.5, Gene, NA_character_)
    )
}

# Entry Point 2: DepMap Consolidated
clean_depmap_consensus <- function(file_path) {
  read_tsv(file_path, show_col_types = FALSE) %>%
    # REVERT: Explicit rename
    dplyr::select(Gene = any_of(c("id", "Gene", "OFFICIAL_SYMBOL")), 
                  LFC = any_of(c("neg|lfc", "LFC", "lfc")), 
                  FDR = any_of(c("neg|fdr", "FDR", "fdr"))) %>%
    mutate(
      Category = case_when(
        FDR < 0.20 & LFC < -0.5 ~ "Depleted (Essential)",
        FDR < 0.20 & LFC > 0.5  ~ "Enriched (Resistant)",
        TRUE ~ "Non-Significant"
      ),
      Label = if_else(FDR < 0.05 & LFC < -0.8, Gene, NA_character_)
    )
}