library(tidyverse)
library(httr)
library(jsonlite)

get_drug_targets <- function(clean_df) {
  
  # 1. Filter for Significant Hits
  hits <- clean_df %>%
    filter(FDR < 0.05) %>%
    arrange(FDR) %>%
    head(40) %>% 
    pull(Gene)
  
  if(length(hits) == 0) {
    message(">>> No significant hits found.")
    return(tibble(Gene = character(), Drug = character(), Note = "No hits"))
  }
  
  message(">>> Querying dGIdb (v5 GraphQL) for ", length(hits), " genes...")
  
  # 2. Construct GraphQL Query
  gene_list_json <- toJSON(hits, auto_unbox = FALSE)
  
  query_string <- sprintf('
    {
      genes(names: %s) {
        nodes {
          name
          interactions {
            drug {
              name
              approved
            }
            interactionScore
            interactionTypes {
              type
            }
            sources {
              sourceDbName
            }
          }
        }
      }
    }
  ', gene_list_json)
  
  # 3. Execute POST Request
  url <- "https://dgidb.org/api/graphql"
  
  resp <- POST(
    url,
    body = list(query = query_string),
    encode = "json",
    add_headers(`Content-Type` = "application/json")
  )
  
  if (status_code(resp) != 200) {
    warning("API Request Failed: ", status_code(resp))
    return(NULL)
  }
  
  # 4. Parse Response (THE FIX: simplifyVector = FALSE)
  # This ensures 'nodes' is a List of Objects, not a Data Frame
  parsed <- httr::content(resp, as = "text", encoding = "UTF-8") %>% 
    fromJSON(simplifyVector = FALSE)
  
  nodes <- parsed$data$genes$nodes
  
  if (is.null(nodes) || length(nodes) == 0) {
    message(">>> API returned no data.")
    return(tibble(Gene = hits, Drug = "None"))
  }
  
  # 5. Flatten Nested Lists
  results <- map_dfr(nodes, function(node) {
    
    gene_name <- node$name
    interactions <- node$interactions
    
    if (is.null(interactions) || length(interactions) == 0) return(NULL)
    
    # Iterate through interactions for this gene
    map_dfr(interactions, function(intr) {
      
      # Extract interaction types safely
      types <- if(!is.null(intr$interactionTypes)) {
        map_chr(intr$interactionTypes, ~ .x$type) %>% paste(collapse = ", ")
      } else { NA_character_ }
      
      # Extract sources safely
      srcs <- if(!is.null(intr$sources)) {
        map_chr(intr$sources, ~ .x$sourceDbName) %>% paste(collapse = ", ")
      } else { NA_character_ }
      
      tibble(
        Gene = gene_name,
        Drug = intr$drug$name,
        Approved = if(is.null(intr$drug$approved)) FALSE else intr$drug$approved,
        Interaction = types,
        Score = if(is.null(intr$interactionScore)) 0 else intr$interactionScore,
        Sources = srcs
      )
    })
  })
  
  if (nrow(results) == 0) {
    return(tibble(Gene = hits, Drug = "None"))
  }
  
  # 6. Prioritize Inhibitors
  final_table <- results %>%
    mutate(
      Is_Inhibitor = str_detect(Interaction, regex("inhibitor|antagonist|suppressor", ignore_case = TRUE))
    ) %>%
    arrange(desc(Is_Inhibitor), desc(Score))
  
  message(">>> Found ", nrow(final_table), " interactions.")
  return(final_table)
}