# =============================================================================
# 03_build_networks.R — Build igraph networks and lookup indices
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

#' Preprocess STRING edges: deduplicate, filter by score, handle protein IDs (Bug 2)
preprocess_string_edges <- function(edges, score_cutoff, config) {
  if (is.null(edges)) return(NULL)

  log_message(sprintf("Preprocessing STRING edges: %d raw edges", nrow(edges)))

  # Detect and map protein IDs (Bug 2)
  # If columns look like STRING protein IDs, attempt mapping
  # For now, assume gene symbols are already used (simulation data uses symbols)
  if (is_protein_id(edges$geneA[1]) || is_protein_id(edges$geneB[1])) {
    log_message("WARNING: STRING protein IDs detected. Ensure gene symbol mapping is applied.", "WARN")
  }

  # Standardize symbols
  edges$geneA <- standardize_gene_symbols(edges$geneA)
  edges$geneB <- standardize_gene_symbols(edges$geneB)
  edges <- edges[!is.na(edges$geneA) & !is.na(edges$geneB), ]

  # Remove self-loops
  if (config$networks$string$remove_self_loops) {
    edges <- edges[edges$geneA != edges$geneB, ]
  }

  # Collapse duplicate gene pairs (keep max score) - Bug 16 fix
  if (config$networks$string$collapse_duplicates_by == "max_score") {
    edges <- edges %>%
      dplyr::mutate(
        g1 = pmin(geneA, geneB),
        g2 = pmax(geneA, geneB)
      ) %>%
      dplyr::group_by(g1, g2) %>%
      dplyr::summarise(combined_score = max(combined_score, na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::rename(geneA = g1, geneB = g2) %>%
      as.data.frame()
  }

  # Filter by score
  edges <- edges[edges$combined_score >= score_cutoff, ]

  log_message(sprintf("STRING edges after preprocessing: %d", nrow(edges)))
  edges
}

#' Build STRING functional association graph
build_string_graph <- function(string_edges, config) {
  if (is.null(string_edges)) {
    log_message("No STRING edges provided, skipping STRING graph construction", "WARN")
    return(NULL)
  }

  edges_clean <- preprocess_string_edges(string_edges,
                                         config$networks$string$score_cutoff,
                                         config)

  if (nrow(edges_clean) == 0) {
    log_message("No STRING edges passed filter, returning NULL graph", "WARN")
    return(NULL)
  }

  edges_use <- edges_clean %>%
    dplyr::mutate(
      edge_confidence  = combined_score / 1000,
      distance_weight  = 1 / edge_confidence
    )

  # Remove duplicates
  edges_use <- edges_use %>%
    dplyr::distinct(geneA, geneB, .keep_all = TRUE)

  graph <- igraph::graph_from_data_frame(
    edges_use[, c("geneA", "geneB", "edge_confidence", "distance_weight")],
    directed = FALSE
  )

  log_message(sprintf("STRING graph: %d vertices, %d edges",
                      igraph::vcount(graph), igraph::ecount(graph)))
  graph
}

#' Build directional signaling graph from OmniPath/SIGNOR edges
build_signaling_graph <- function(signaling_edges, config) {
  if (is.null(signaling_edges)) {
    log_message("No signaling edges provided, skipping signaling graph construction", "WARN")
    return(NULL)
  }

  log_message(sprintf("Building signaling graph from %d edges", nrow(signaling_edges)))

  edges_use <- signaling_edges %>%
    dplyr::mutate(
      edge_confidence = ifelse(is.na(confidence), 1, as.numeric(confidence)),
      distance_weight = 1 / edge_confidence
    )

  # Ensure column names
  if (!"source_gene" %in% names(edges_use) || !"target_gene" %in% names(edges_use)) {
    log_message("Signaling edges missing source_gene/target_gene columns, skipping", "ERROR")
    return(NULL)
  }

  # Standardize
  edges_use$source_gene <- standardize_gene_symbols(edges_use$source_gene)
  edges_use$target_gene <- standardize_gene_symbols(edges_use$target_gene)

  # Remove self-loops
  if (config$networks$signaling$remove_self_loops) {
    edges_use <- edges_use[edges_use$source_gene != edges_use$target_gene, ]
  }

  if (nrow(edges_use) == 0) return(NULL)

  graph <- igraph::graph_from_data_frame(
    edges_use[, c("source_gene", "target_gene", "edge_confidence", "distance_weight")],
    directed = TRUE
  )

  log_message(sprintf("Signaling graph: %d vertices, %d edges",
                      igraph::vcount(graph), igraph::ecount(graph)))
  graph
}

#' Build Reactome pathway index: gene -> pathway lookup (and reverse)
build_reactome_index <- function(reactome_df) {
  if (is.null(reactome_df)) return(NULL)

  if (!"gene" %in% names(reactome_df) || !"pathway" %in% names(reactome_df)) {
    log_message("Reactome data missing 'gene' or 'pathway' columns, skipping", "WARN")
    return(NULL)
  }

  reactome_df$gene <- standardize_gene_symbols(reactome_df$gene)
  reactome_df <- reactome_df[!is.na(reactome_df$gene), ]

  gene_to_pathways <- split(reactome_df$pathway, reactome_df$gene)
  pathway_to_genes <- split(reactome_df$gene, reactome_df$pathway)

  log_message(sprintf("Reactome index: %d genes, %d pathways",
                      length(gene_to_pathways), length(pathway_to_genes)))
  list(gene_to_pathways = gene_to_pathways, pathway_to_genes = pathway_to_genes)
}

#' Build CORUM complex index
build_complex_index <- function(corum_df) {
  if (is.null(corum_df)) return(NULL)

  if (!"gene" %in% names(corum_df) || !"complex_id" %in% names(corum_df)) {
    log_message("CORUM data missing required columns, skipping", "WARN")
    return(NULL)
  }

  corum_df$gene <- standardize_gene_symbols(corum_df$gene)
  corum_df <- corum_df[!is.na(corum_df$gene), ]

  gene_to_complexes <- split(corum_df$complex_id, corum_df$gene)
  complex_to_genes  <- split(corum_df$gene, corum_df$complex_id)

  log_message(sprintf("CORUM index: %d genes, %d complexes",
                      length(gene_to_complexes), length(complex_to_genes)))
  list(gene_to_complexes = gene_to_complexes, complex_to_genes = complex_to_genes)
}

#' Build SL pair set for O(1) lookup
build_sl_pair_set <- function(sl_db) {
  if (is.null(sl_db)) return(NULL)

  if (!"geneA" %in% names(sl_db) || !"geneB" %in% names(sl_db)) {
    log_message("SL DB missing geneA/geneB columns, skipping", "WARN")
    return(NULL)
  }

  sl_db$geneA <- standardize_gene_symbols(sl_db$geneA)
  sl_db$geneB <- standardize_gene_symbols(sl_db$geneB)

  pair_keys <- make_pair_key(sl_db$geneA, sl_db$geneB)
  log_message(sprintf("SL pair set: %d known pairs", length(unique(pair_keys))))
  unique(pair_keys)
}

#' Build drug target index
build_drug_target_index <- function(drug_df) {
  if (is.null(drug_df)) return(NULL)

  if (!"gene" %in% names(drug_df) || !"drug_category" %in% names(drug_df)) {
    log_message("Drug data missing required columns, skipping", "WARN")
    return(NULL)
  }

  drug_df$gene <- standardize_gene_symbols(drug_df$gene)
  drug_df <- drug_df[!is.na(drug_df$gene), ]

  # Keep best category per gene
  category_order <- c("FDA_approved", "clinical", "preclinical", "druggable_class", "unknown")
  drug_df$category_rank <- match(drug_df$drug_category, category_order)
  drug_df <- drug_df %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_min(category_rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene, drug_category, drug_name)

  index <- split(drug_df, drug_df$gene)

  log_message(sprintf("Drug index: %d targets", length(index)))
  index
}

#' Build all networks and indices
build_all_networks <- function(data_list) {
  log_message("=== Step 3: Build Networks ===")

  config <- data_list$config

  string_graph    <- build_string_graph(data_list$string_edges, config)
  signaling_graph <- build_signaling_graph(data_list$signaling_edges, config)
  reactome_index  <- build_reactome_index(data_list$reactome_df)
  complex_index   <- build_complex_index(data_list$corum_df)
  sl_pair_set     <- build_sl_pair_set(data_list$sl_db)
  drug_index      <- build_drug_target_index(data_list$drug_df)

  network_objects <- list(
    string_graph    = string_graph,
    signaling_graph = signaling_graph,
    reactome_index  = reactome_index,
    complex_index   = complex_index,
    sl_pair_set     = sl_pair_set,
    drug_index      = drug_index
  )

  log_message("=== Network building complete ===")
  network_objects
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  source("R/01_prepare_data.R")
  data_list <- prepare_data("config/config.yaml")
  network_objects <- build_all_networks(data_list)
  log_message("Network objects:")
  log_message(paste(names(network_objects), collapse = ", "))
}
