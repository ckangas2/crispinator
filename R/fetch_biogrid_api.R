library(httr)
library(jsonlite)
library(tidyverse)

BASE_URL <- "https://orcsws.thebiogrid.org"

# 1. SEARCH FUNCTION (Unchanged)
search_by_pubmed <- function(pmid, access_key) {
  message(">>> Querying API for PubMed ID: ", pmid)
  url <- paste0(BASE_URL, "/screens/")
  resp <- GET(url, query = list(accesskey = access_key, format = "json", pubmedID = pmid))
  
  if (status_code(resp) != 200) stop("Search Failed: ", status_code(resp))
  
  raw_text <- content(resp, "text", encoding = "UTF-8")
  if (raw_text == "[]" || raw_text == "") return(NULL)
  
  data <- fromJSON(raw_text) %>% as_tibble()
  cols <- colnames(data)
  
  id_col    <- grep("^id$|screen.*id", cols, ignore.case=T, value=T)[1]
  pheno_col <- grep("phenotype|condition|virus", cols, ignore.case=T, value=T)[1]
  
  return(data %>% select(screen_id = all_of(id_col), virus = any_of(pheno_col)))
}

# 2. FETCH DATA (Fixed Gene Detection)
fetch_screen_data <- function(screen_id, access_key, output_dir = "data/raw") {
  
  message(">>> Downloading Screen ID: ", screen_id)
  url <- paste0(BASE_URL, "/screen/", screen_id)
  
  resp <- GET(url, query = list(accesskey = access_key, format = "tab", header = "yes"))
  if (status_code(resp) != 200) stop("Download Failed")
  
  scores <- read_tsv(content(resp, "text", encoding="UTF-8"), show_col_types = FALSE)
  cols <- colnames(scores)
  
  # --- FIXING THE GENE ID BUG ---
  # 1. Try Regex for Symbol
  gene_col <- grep("OFFICIAL_SYMBOL|GENE_SYMBOL", cols, ignore.case=T, value=T)[1]
  
  # 2. Fallback: If regex failed or picked SCREEN_ID, force Column 3 (Standard BioGRID)
  if (is.na(gene_col) || gene_col == "SCREEN_ID") {
    message("    Regex failed. Defaulting to Column 3 for Gene Symbols.")
    gene_col <- cols[3] 
  }
  
  # LFC / Score Detection
  lfc_col <- grep("SCORE.*1|LFC|Z.*Score|Log", cols, ignore.case=T, value=T)[1]
  fdr_col <- grep("SCORE.*2|FDR|P.*Val", cols, ignore.case=T, value=T)[1]
  
  message("    Mapping: [ID: ", gene_col, "] [Score: ", lfc_col, "]")
  
  clean_df <- scores %>%
    dplyr::select(
      id = all_of(gene_col),
      `neg|lfc` = all_of(lfc_col),
      `neg|fdr` = all_of(fdr_col)
    ) %>%
    mutate(
      `neg|lfc` = as.numeric(`neg|lfc`),
      `neg|fdr` = as.numeric(`neg|fdr`)
    ) %>%
    filter(!is.na(`neg|lfc`))
  
  fname <- paste0("biogrid_ORCS_screen_", screen_id, ".txt")
  out_path <- file.path(output_dir, fname)
  write_tsv(clean_df, out_path)
  return(out_path)
}