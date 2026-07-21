# R/04_enrichment.R

# Load these explicitly inside the function or Ensure they are loaded in _targets.R
# library(clusterProfiler)
# library(org.Mm.eg.db)

perform_enrichment <- function(data) {
  sig_genes <- data %>%
    dplyr::filter(Category == "Depleted (Essential)") %>%
    dplyr::pull(Gene)
  
  if(length(sig_genes) == 0) return(NULL)
  
  # --- AUTODETECT LOGIC ---
  # Check if the genes are mostly uppercase (Human) or Titlecase (Mouse)
  # or check for specific known mouse-only symbols.
  is_human <- any(data$Gene %in% c("TP53", "MYC", "EGFR", "POLR2A"))
  org_db <- if(is_human) "org.Hs.eg.db" else "org.Mm.eg.db"
  
  message(">>> Species Autodetect: Using ", org_db)
  # -------------------------
  
  gene_map <- bitr(sig_genes, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org_db) # Dynamic
  
  go_results <- enrichGO(
    gene = gene_map$ENTREZID,
    OrgDb = org_db,             # Dynamic
    keyType = "ENTREZID",
    ont = "BP", 
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  return(go_results)
}

plot_enrichment <- function(enrich_result) {
  
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    # Return an empty placeholder plot if nothing found
    return(ggplot() + annotate("text", x=1, y=1, label="No Enriched Pathways"))
  }
  
  # Generate Dot Plot
  # Shows top 10 pathways
  p <- dotplot(enrich_result, showCategory = 10) + 
    ggtitle("Enriched Biological Processes (Essential Genes)") +
    theme(axis.text.y = element_text(size = 10)) # Make pathway names readable
  
  return(p)
}