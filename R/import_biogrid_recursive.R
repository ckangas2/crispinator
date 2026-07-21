# R/import_biogrid_recursive.R
library(tidyverse)
library(vroom)

ingest_biogrid_folder <- function(raw_dir = "data/raw") {
  
  message(">>> Scanning for BioGRID subfolders...")
  screen_files <- list.files(raw_dir, pattern = "screen.tab.txt$", 
                             recursive = TRUE, full.names = TRUE)
  
  if(length(screen_files) == 0) stop("No BioGRID screen files found.")
  
  message(">>> Found ", length(screen_files), " potential screen(s). Processing...")
  
  for (screen_path in screen_files) {
    
    # 1. Load Sibling Index
    parent_folder <- dirname(screen_path)
    index_path    <- list.files(parent_folder, pattern = "index.tab.txt$", full.names = TRUE)
    
    if (length(index_path) == 0) next
    
    # 2. Dynamic Metadata Lookup
    index_data <- vroom(index_path, show_col_types = FALSE)
    cols <- colnames(index_data)
    
    # Find Identifier Columns
    id_col    <- grep("Screen.*ID", cols, ignore.case = TRUE, value = TRUE)[1]
    auth_col  <- grep("Author", cols, ignore.case = TRUE, value = TRUE)[1]
    cell_col  <- grep("Cell.*Line", cols, ignore.case = TRUE, value = TRUE)[1]
    pheno_col <- grep("Phenotype|Condition", cols, ignore.case = TRUE, value = TRUE)[1]
    
    if (is.na(id_col)) { warning("Skipping ", basename(screen_path), ": No ID col."); next }
    
    screen_id_match <- str_extract(basename(screen_path), "SCREEN_(\\d+)", group = 1)
    if (is.na(screen_id_match)) next
    
    # Extract row
    meta <- index_data %>% 
      filter(as.character(.data[[id_col]]) == as.character(screen_id_match))
    
    if (nrow(meta) == 0) { warning("Metadata not found for ", screen_id_match); next }
    
    # --- ROBUST SCORE TYPE DETECTION ---
    # Helper to get type string from index
    get_type_val <- function(n) {
      pat <- paste0("Score.*", n, ".*Type")
      col_name <- grep(pat, cols, ignore.case = TRUE, value = TRUE)[1]
      if(is.na(col_name)) return("UNKNOWN")
      return(as.character(meta[[col_name]]))
    }
    
    s1_type <- get_type_val(1)
    s2_type <- get_type_val(2)
    s3_type <- get_type_val(3)
    s4_type <- get_type_val(4)
    
    message(paste0("    -> Scores Identified: [1:", s1_type, "] [2:", s2_type, "] [3:", s3_type, "] [4:", s4_type, "]"))
    
    # --- THE FIX: ADD RRA AND STARS TO REGEX ---
    # 1. Determine "Effect" Column (X-Axis)
    # Now accepts: Log, Fold, Z-Score, RRA, STARS, Score
    effect_regex <- regex("Log|Fold|Z[- ]?Score|RRA|STARS|Score", ignore_case=T)
    
    lfc_source <- case_when(
      str_detect(s3_type, effect_regex) & !str_detect(s3_type, "FDR|P[- ]?Val") ~ "SCORE.3",
      str_detect(s2_type, effect_regex) & !str_detect(s2_type, "FDR|P[- ]?Val") ~ "SCORE.2",
      str_detect(s1_type, effect_regex) & !str_detect(s1_type, "Rank") ~ "SCORE.1",
      TRUE ~ "UNKNOWN"
    )
    
    # 2. Determine "Significance" Column (Y-Axis)
    sig_regex <- regex("FDR|P[- ]?Value|Significance", ignore_case=T)
    
    fdr_source <- case_when(
      str_detect(s2_type, sig_regex) ~ "SCORE.2",
      str_detect(s3_type, sig_regex) ~ "SCORE.3",
      str_detect(s4_type, sig_regex) ~ "SCORE.4",
      TRUE ~ "UNKNOWN"
    )
    
    message("    -> Mapping Decision: LFC = ", lfc_source, " | FDR = ", fdr_source)
    
    if(lfc_source == "UNKNOWN" || fdr_source == "UNKNOWN") {
      warning("Could not automatically determine score columns for ", screen_id_match, ". Skipping.")
      next
    }
    
    # 3. Clean Filename
    get_val <- function(col, def) if(is.na(col)) def else { v <- meta[[col]]; if(is.null(v) || is.na(v)) def else v }
    
    clean_author <- str_replace_all(get_val(auth_col, "UnkAuth"), "[^[:alnum:]]", "")
    clean_cell   <- str_replace_all(get_val(cell_col, "UnkCell"), "[^[:alnum:]]", "")
    clean_virus  <- str_replace_all(get_val(pheno_col, "UnkVirus"), "[^[:alnum:]]", "")
    
    new_filename <- paste0("biogrid_", clean_author, "_", clean_cell, "_", clean_virus, "_", screen_id_match, ".txt")
    dest_path    <- file.path(raw_dir, new_filename)
    
    # 4. Save
    raw_data <- vroom(screen_path, show_col_types = FALSE)
    
    clean_df <- raw_data %>%
      dplyr::select(
        id = OFFICIAL_SYMBOL,
        `neg|lfc` = all_of(lfc_source),
        `neg|fdr` = all_of(fdr_source)
      ) %>%
      mutate(
        `neg|lfc` = as.numeric(`neg|lfc`),
        `neg|fdr` = as.numeric(`neg|fdr`),
        `pos|fdr` = 1.0
      ) %>%
      filter(!is.na(`neg|lfc`))
    
    write_tsv(clean_df, dest_path)
    message("    -> Saved.\n")
  }
  message(">>> Ingestion Complete.")
}