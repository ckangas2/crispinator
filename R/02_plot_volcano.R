library(tidyverse)
library(ggrepel)

# R/02_plot_volcano.R

plot_volcano <- function(input_df) {
  
  # --- 1. DIRECT COLUMN CHECK (No Renaming) ---
  # We expect: Gene, LFC, FDR (Output from 01_clean_crispr.R)
  
  if (!"Gene" %in% names(input_df)) stop("Plotter Error: Input missing 'Gene' column.")
  if (!"LFC" %in% names(input_df))  stop("Plotter Error: Input missing 'LFC' column.")
  if (!"FDR" %in% names(input_df))  stop("Plotter Error: Input missing 'FDR' column.")
  
  # --- 2. MATH SAFETY (The Infinity Fix) ---
  df <- input_df %>%
    mutate(
      LFC = as.numeric(LFC),
      FDR = as.numeric(FDR),
      # Cap FDR to avoid Infinity (-log10(0) = Inf)
      safe_fdr = pmax(FDR, 1e-300),
      # Cap LFC to avoid log(0) if in Rank Mode
      safe_lfc = pmax(LFC, 1e-20) 
    )
  
  # --- 3. AUTO-DETECT DATA TYPE ---
  # RRA (MAGeCK) scores are 0 to 1. LFC scores are usually -5 to +5.
  vals <- na.omit(df$safe_lfc)
  is_rra <- length(vals) > 0 && min(vals) >= 0 && max(vals) <= 1.05
  
  if (is_rra) {
    message(">>> DETECTED RRA/RANK SCORES. Switching to Rank Mode.")
    
    plot_data <- df %>%
      mutate(
        # RRA Transformation: -log10(Score) puts Best Ranks (0) on the Right
        X_Val = -log10(safe_lfc),
        Y_Val = -log10(safe_fdr),
        
        Category = case_when(
          FDR < 0.05 & X_Val > 3 ~ "Significant Hit", 
          TRUE ~ "Non-Significant"
        ),
        Label = if_else(Category == "Significant Hit", Gene, NA_character_)
      )
    x_label <- "Gene Rank Score (-log10 RRA)"
    
  } else {
    # STANDARD LFC MODE
    plot_data <- df %>%
      mutate(
        X_Val = LFC,
        Y_Val = -log10(safe_fdr),
        
        Category = case_when(
          FDR < 0.05 & abs(X_Val) > 1 ~ "Significant Hit",
          TRUE ~ "Non-Significant"
        ),
        Label = if_else(Category == "Significant Hit", Gene, NA_character_)
      )
    x_label <- "Log2 Fold Change"
  }
  
  # --- 4. PLOTTING ---
  p <- ggplot(plot_data, aes(x = X_Val, y = Y_Val, color = Category)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Non-Significant" = "grey85", "Significant Hit" = "dodgerblue4")) +
    
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20, box.padding = 0.5) +
    
    labs(
      title = if(is_rra) "CRISPR Screen (Rank Mode)" else "CRISPR Screen (LFC Mode)",
      subtitle = "Values capped to avoid Infinity",
      x = x_label,
      y = "-Log10 FDR"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}