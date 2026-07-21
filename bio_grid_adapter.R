# R/import_biogrid_recursive.R
library(tidyverse)
library(vroom)

ingest_biogrid_folder <- function(raw_dir = "data/raw") {
  
  message(">>> Scanning for BioGRID subfolders...")
  
  # 1. Find all BioGRID Screen files (Recursively)
  # This looks inside folders like 'BIOGRID-ORCS-7982...'
  screen_files <- list.files(raw_dir, pattern = "screen.tab.txt$", 
                             recursive = TRUE, full.names = TRUE)
  
  if(length(screen_files) == 0) stop("No BioGRID screen files found in subfolders.")
  
  message(">>> Found ", length(screen_files), " potential screen(s). Processing...")
  
  for (screen_path in screen_files) {
    
    # 2. Find the Sibling Index File
    # We assume the index is in the SAME folder as the screen file
    parent_folder <- dirname(screen_path)
    index_path    <- list.files(parent_folder, pattern = "index.tab.txt$", full.names = TRUE)
    
    if (length(index_path) == 0) {
      warning("Skipping ", basename(screen_path), ": No matching index file found.")
      next
    }
    
    # 3. Read Metadata to Identify the Study
    # The index file is tiny, so vroom is instant.
    index_data <- vroom(index_path, show_col_types = FALSE)
    
    # BioGRID Index contains info for the specific screen ID in the filename
    # Extract Screen ID from filename (e.g., "...SCREEN_1431-...")
    screen_id_match <- str_extract(basename(screen_path), "SCREEN_(\\d+)", group = 1)
    
    if (is.na(screen_id_match)) {
      warning("Could not parse Screen ID from filename: ", basename(screen_path))
      next
    }
    
    # Get metadata for THIS specific screen ID
    meta <- index_data %>% filter(SCREEN_ID == screen_id_match)
    
    if (nrow(meta) == 0) {
      warning("Metadata for Screen ", screen_id_match, " not found in index.")
      next
    }
    
    # 4. Generate the Clean "Standard" Name
    # Format: biogrid_Author_CellLine_ID_converted.txt
    clean_author <- str_replace_all(meta$PUBLICATION_FIRST_AUTHOR, "[^[:alnum:]]", "")
    clean_cell   <- str_replace_all(meta$Cell_Line_Name, "[^[:alnum:]]", "")
    clean_virus  <- str_replace_all(meta$Phenotype, "[^[:alnum:]]", "") # e.g. SARSCoV2
    
    new_filename <- paste0("biogrid_", clean_author, "_", clean_cell, "_", clean_virus, "_", screen_id_match, ".txt")
    dest_path    <- file.path(raw_dir, new_filename)
    
    message("\n---------------------------------------------------")
    message("   Study:     ", meta$PUBLICATION_FIRST_AUTHOR, " (", meta$PUBLICATION_YEAR, ")")
    message("   Cell Line: ", meta$Cell_Line_Name)
    message("   Target:    ", meta$Phenotype)
    message("   Score Type:", meta$SCORE_1_TYPE)
    message("   -> Saving to: ", new_filename)
    
    # 5. Load, Clean, and Move
    # This is where we apply the column mapping we discussed earlier
    raw_data <- vroom(screen_path, show_col_types = FALSE)
    
    clean_df <- raw_data %>%
      dplyr::select(
        id = OFFICIAL_SYMBOL,
        # Dynamically map scores. Usually SCORE.1 is effect, SCORE.2 is FDR.
        `neg|lfc` = SCORE.1,
        `neg|fdr` = SCORE.2
      ) %>%
      mutate(
        `neg|lfc` = as.numeric(`neg|lfc`),
        `neg|fdr` = as.numeric(`neg|fdr`),
        `pos|fdr` = 1.0
      ) %>%
      filter(!is.na(`neg|lfc`))
    
    write_tsv(clean_df, dest_path)
  }
  
  message("\n>>> Ingestion Complete. Files are ready for the pipeline.")
}