# _targets.R
library(targets)
library(tarchetypes)

tar_source() 

# 1. Add 'conflicted' to the packages
tar_option_set(packages = c(
  "tidyverse", "conflicted", "ggrepel", 
  "clusterProfiler", "enrichplot",
  "org.Hs.eg.db", # Human
  "org.Mm.eg.db",  # Mouse
  "httr",      # <--- Explicitly added
  "jsonlite"
))

# 2. RESOLVE CONFLICTS GLOBALLY
# This ensures select() always means dplyr::select() in every script
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

if (file.exists("_targets_config.rds")) {
  config <- readRDS("_targets_config.rds")
} else {
  config <- list(
    input_path = "data/raw/my_first_screen.gene_summary.txt",
    output_dir = "output/figures"
  )
}


list(
  # 1. Define Input
  tar_target(file_input, config$input_path, format = "file"),
  
  # 2. Extract Name (e.g., "biogrid_ORCS_screen_2137")
  tar_target(sample_name, tools::file_path_sans_ext(basename(file_input))),
  
  # 3. Clean Data 
  # (clean_mageck_lab now handles BioGRID names automatically)
  tar_target(
    clean_data,
    if (str_detect(file_input, "depmap")) {
      clean_depmap_consensus(file_input)
    } else {
      clean_mageck_lab(file_input)
    }
  ),
  
  # 4. Volcano Plot (Adaptive RRA detection)
  tar_target(volcano_plot, plot_volcano(clean_data)),
  
  # 5. Save Plot (Respects GUI Output Folder)
  tar_target(
    saved_volcano, 
    save_plot_timestamp(volcano_plot, paste0(sample_name, "_volcano"), config$output_dir), 
    format = "file"
  ),
  
  # 6. Save Hits CSV
  tar_target(
    saved_hits, 
    save_hits_csv(clean_data, paste0(sample_name, "_hits"), config$output_dir), 
    format = "file"
  ),
  
  # 7. Enrichment Analysis (Optional but included)
  tar_target(enrichment_results, perform_enrichment(clean_data)),
  tar_target(pathway_plot, plot_enrichment(enrichment_results)),
  
  tar_target(
    saved_pathway_plot, 
    save_plot_timestamp(pathway_plot, paste0(sample_name, "_enrichment"), config$output_dir), 
    format = "file"
  ),
  
  # 8. Comparison to DepMap (Safety Check)
  tar_target(ref_file, "data/raw/depmap_A549_standardized.txt", format = "file"),
  
  tar_target(
    safety_check,
    compare_viral_vs_depmap(clean_data, ref_file, config$output_dir)
  ),
  
  # 9. Drug Discovery (GraphQL)
  tar_target(drug_candidates, get_drug_targets(clean_data)),
  
  # 10. Save Drug List (Respects GUI Output Folder)
  tar_target(
    saved_drugs,
    command = {
      # Use file.path with config$output_dir to obey the GUI selection
      outfile <- file.path(config$output_dir, paste0(sample_name, "_drugs.csv"))
      
      # Ensure directory exists
      dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
      
      # Write
      write_csv(drug_candidates, outfile)
      
      # Return path for tracking
      outfile 
    },
    format = "file"
  )
)