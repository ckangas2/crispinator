library(tidyverse)
library(vroom)

# --- CONFIGURATION ---
# Change this to any cell line. The script will find the best match.
search_term <- "BLaER1" 
# ---------------------

# Update this path to your local project directory
raw_path <- "C:/Users/Chase Kangas/Documents/Everything_R/CRISPR_SCREEN_PIPELINE/data/raw"

# Identify files (handles versioning like 25Q4 or 26Q1)
effect_file <- list.files(raw_path, pattern = "CRISPRGeneEffect", full.names = TRUE)
dep_file    <- list.files(raw_path, pattern = "CRISPRGeneDependency", full.names = TRUE)
model_file  <- list.files(raw_path, pattern = "Model.csv", full.names = TRUE)

# 1. Hierarchical Search Logic
if(length(model_file) == 0) stop("CRITICAL ERROR: Model.csv missing from data/raw/")

models <- vroom(model_file, show_col_types = FALSE)

# Prepare alphanumeric 'clean' versions for fuzzy matching
search_clean <- str_replace_all(search_term, "[^[:alnum:]]", "")

candidates <- models %>%
  mutate(clean_name = str_replace_all(CellLineName, "[^[:alnum:]]", "")) %>%
  filter(str_detect(clean_name, fixed(search_clean, ignore_case = TRUE)))

# Level 1: Look for an EXACT alphanumeric match (e.g., "A375" finds "A-375")
target_line <- candidates %>% 
  filter(clean_name == search_clean) %>% 
  slice(1)

# Level 2: If no exact match, look for the shortest name containing the term
# (Prevents "A375" from grabbing "A375-SKIN-CJ1" if a shorter match exists)
if (nrow(target_line) == 0 && nrow(candidates) > 0) {
  target_line <- candidates %>%
    mutate(name_length = nchar(CellLineName)) %>%
    arrange(name_length) %>%
    slice(1)
  message(">>> No exact match. Choosing shortest candidate: ", target_line$CellLineName)
}

# ERROR HANDLING: Hard stop if search fails completely
if (nrow(target_line) == 0) {
  stop("\n!!! SEARCH FAILED !!!\nNo matches for '", search_term, "' found in Model.csv.")
}

target_id <- target_line$ModelID
# Create a safe filename (replaces spaces/dashes with underscores)
clean_line_name <- str_replace_all(target_line$CellLineName, "[^[:alnum:]]", "_")

message(">>> Target Found: ", target_line$CellLineName)
message(">>> Model ID:    ", target_id)
message(">>> Tissue:      ", target_line$OncodePrimaryDisease)

# 2. Memory-Efficient Extraction
message("Reading Effect Matrix...")
ge <- vroom(effect_file, show_col_types = FALSE) %>%
  dplyr::filter(...1 == target_id) %>%
  tidyr::pivot_longer(-...1, names_to = "Gene", values_to = "LFC")

message("Reading Dependency Matrix...")
gd <- vroom(dep_file, show_col_types = FALSE) %>%
  dplyr::filter(...1 == target_id) %>%
  tidyr::pivot_longer(-...1, names_to = "Gene", values_to = "Prob")

# 3. Standardize and Clean
message("Finalizing pipeline-ready file...")
depmap_output <- ge %>%
  inner_join(gd, by = c("...1", "Gene")) %>%
  mutate(
    # Clean "AAMP (14)" -> "AAMP" for pathway analysis compatibility
    id = str_remove(Gene, " \\(\\d+\\)$"),
    `neg|lfc` = LFC,
    # Convert Probability (0-1) to Pseudo-FDR (1-Prob)
    `neg|fdr` = 1 - Prob,
    `pos|fdr` = 1.0 
  ) %>%
  dplyr::select(id, `neg|lfc`, `neg|fdr`, `pos|fdr`)

# 4. Save to Project Folder
output_filename <- file.path(raw_path, paste0("depmap_", clean_line_name, "_standardized.txt"))
write_tsv(depmap_output, output_filename)

message("\nSUCCESS: Saved ", nrow(depmap_output), " genes to: ", output_filename)