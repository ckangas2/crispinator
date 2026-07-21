# R/05_depmap_filter.R
library(tidyverse)

annotate_viral_vs_essential <- function(viral_df, depmap_df) {
  
  # 1. Define 'Essential' in the Host (DepMap Baseline)
  # Chronos < -0.5 is a standard cutoff for "Moderate Essentiality"
  host_essentials <- depmap_df %>%
    filter(`neg|lfc` < -0.5) %>%
    select(Gene = id, DepMap_Score = `neg|lfc`)
  
  # 2. Join with your Viral Data
  annotated_df <- viral_df %>%
    left_join(host_essentials, by = "Gene") %>%
    mutate(
      # Create a new classification column
      Hit_Type = case_when(
        # If it's a hit in your screen AND essential in DepMap
        Category == "Depleted (Essential)" & !is.na(DepMap_Score) ~ "Common Essential (Housekeeping)",
        
        # If it's a hit in your screen but NOT in DepMap
        Category == "Depleted (Essential)" & is.na(DepMap_Score) ~ "Viral-Specific Dependency",
        
        # Keep original categories for everything else
        TRUE ~ Category
      )
    )
  
  return(annotated_df)
}