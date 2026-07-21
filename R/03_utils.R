# R/03_utils.R

save_plot_timestamp <- function(plot_object, base_prefix, output_dir) {
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
  # Name format: sample_name_type_timestamp.png
  filename <- file.path(output_dir, paste0(base_prefix, "_", timestamp_str, ".png"))
  
  ggsave(filename, plot = plot_object, device = "png", 
         width = 8, height = 6, dpi = 300, bg = "white")
  
  return(filename)
}

save_hits_csv <- function(data, base_prefix, output_dir) {
  hits_df <- data %>%
    dplyr::filter(Category != "Non-Significant") %>%
    dplyr::arrange(FDR) %>%
    dplyr::select(Gene, Category, LFC, FDR, Label)
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- file.path(output_dir, paste0(base_prefix, "_", timestamp_str, ".csv"))
  
  write_csv(hits_df, filename)
  return(filename)
}