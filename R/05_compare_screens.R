library(tidyverse)
library(ggrepel)

compare_viral_vs_depmap <- function(viral_data, depmap_path, output_dir) {
  
  # 1. Load DepMap Reference (This IS a raw file, so we select 'id')
  depmap_df <- read_tsv(depmap_path, show_col_types = FALSE) %>%
    dplyr::select(Gene = id, DepMap_LFC = `neg|lfc`)
  
  # 2. Prepare Viral Data (This is ALREADY cleaned, so we select 'Gene')
  # The 'clean_data' target has columns: Gene, LFC, FDR, Category
  viral_hits <- viral_data %>%
    dplyr::select(Gene, Viral_LFC = LFC, Viral_FDR = FDR) # <--- FIXED HERE
  
  # 3. Join
  comparison <- viral_hits %>%
    inner_join(depmap_df, by = "Gene") %>%
    mutate(
      # CLASSIFICATION LOGIC:
      Target_Class = case_when(
        Viral_FDR < 0.05 & Viral_LFC > 0.5 & DepMap_LFC > -0.5 ~ "Ideal Target (Safe)",
        Viral_FDR < 0.05 & Viral_LFC > 0.5 & DepMap_LFC <= -0.5 ~ "Toxic Target (Essential)",
        TRUE ~ "Other"
      ),
      Label = if_else(Target_Class != "Other", Gene, NA_character_)
    )
  
  # 4. Generate The "Quadrant Plot"
  p <- ggplot(comparison, aes(x = DepMap_LFC, y = Viral_LFC, color = Target_Class)) +
    geom_point(alpha = 0.6) +
    
    # Add Quadrant Lines
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20) +
    
    scale_color_manual(values = c("Ideal Target (Safe)" = "#2ca02c", # Green
                                  "Toxic Target (Essential)" = "#d62728", # Red
                                  "Other" = "grey90")) +
    
    labs(
      title = "Target Safety Check: Viral Resistance vs. Host Essentiality",
      x = "Host Fitness (DepMap Score)\n< Left is Toxic | Right is Healthy >",
      y = "Viral Resistance (Log2 Fold Change)\n< Up is Resistant | Down is Sensitized >"
    ) +
    theme_minimal()
  
  # 5. Save
  if(!dir.exists(output_dir)) dir.create(output_dir)
  ggsave(file.path(output_dir, paste0("safety_check_", format(Sys.time(), "%H%M%S"), ".png")), p, width = 8, height = 6)
  
  return(comparison)
}