# =============================================================================
# 04_compute_network_features.R — STRING network features for each candidate pair
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

#' Convert shortest path distance to a score (from config)
convert_distance_to_score <- function(distance, config) {
  dmap <- config$network_features$distance_score_map
  dmap <- as.list(dmap)
  dkey <- as.character(distance)
  if (dkey %in% names(dmap)) {
    return(as.numeric(dmap[[dkey]]))
  }
  max_len <- config$network_features$max_path_length
  if (distance <= max_len) return(0.1)
  0
}

#' Compute path confidence from edge scores (Bug 5: vanishing product)
compute_path_confidence <- function(edge_confs, method = "product", max_len = 5) {
  if (length(edge_confs) == 0) return(0)
  if (length(edge_confs) > max_len) return(0)  # Bug 5: cap long paths

  switch(method,
    product = prod(edge_confs),
    geometric_mean = prod(edge_confs)^(1 / length(edge_confs)),
    log_sum = sum(log(edge_confs + 1e-10)),
    prod(edge_confs)  # default
  )
}

#' Compute hub penalty for intermediate genes (Bug 3)
#' Formula: 1 / log2(mean_degree + 2)
compute_hub_penalty <- function(graph, intermediate_genes) {
  if (!is.null(graph) && length(intermediate_genes) > 0) {
    degrees <- igraph::degree(graph, v = intermediate_genes)
    mean_deg <- mean(degrees, na.rm = TRUE)
    return(1 / log2(mean_deg + 2))
  }
  1  # No intermediates = no penalty
}

#' Compute common neighbors between two genes
compute_common_neighbors <- function(graph, geneA, geneB) {
  if (is.null(graph)) return(0)
  vnames <- igraph::V(graph)$name
  if (!(geneA %in% vnames) || !(geneB %in% vnames)) return(0)

  nA <- igraph::neighbors(graph, geneA)$name
  nB <- igraph::neighbors(graph, geneB)$name
  length(intersect(nA, nB))
}

#' Compute network features for a single gene pair
compute_network_features_pair <- function(graph, geneA, geneB, config) {
  hub_enabled <- config$network_features$hub_penalty_enabled
  method      <- config$network_features$path_confidence_method
  max_len     <- config$network_features$max_path_length

  # Default disconnected result
  disconnected <- data.frame(
    connected          = FALSE,
    shortest_distance  = NA_integer_,
    path_genes         = NA_character_,
    path_confidence    = 0,
    common_neighbors   = 0L,
    direct_interaction = FALSE,
    hub_penalty        = 1,
    distance_score     = 0,
    network_score      = 0,
    stringsAsFactors   = FALSE
  )

  if (is.null(graph)) return(disconnected)

  sp <- safe_shortest_path(graph, geneA, geneB,
                           weights = igraph::E(graph)$distance_weight)

  if (is.null(sp) || sp$distance > max_len) {
    return(disconnected)
  }

  intermediates <- if (sp$distance > 1) sp$path_genes[2:(length(sp$path_genes) - 1)] else character(0)

  distance_score  <- convert_distance_to_score(sp$distance, config)
  path_conf       <- compute_path_confidence(sp$edge_confs, method, max_len)
  hub_pen         <- if (hub_enabled) compute_hub_penalty(graph, intermediates) else 1
  network_score   <- distance_score * path_conf * hub_pen
  common_n        <- compute_common_neighbors(graph, geneA, geneB)

  data.frame(
    connected          = TRUE,
    shortest_distance  = sp$distance,
    path_genes         = paste(sp$path_genes, collapse = ";"),
    path_confidence    = path_conf,
    common_neighbors   = common_n,
    direct_interaction = (sp$distance == 1),
    hub_penalty        = hub_pen,
    distance_score     = distance_score,
    network_score      = network_score,
    stringsAsFactors   = FALSE
  )
}

#' Compute network features for all unique pairs in candidate_df
compute_network_features_batch <- function(graph, candidate_df, config) {
  log_message("=== Step 4: Compute Network Features ===")

  if (nrow(candidate_df) == 0) {
    log_message("No candidates to compute network features for", "WARN")
    return(candidate_df)
  }

  # Get unique pairs
  unique_pairs <- candidate_df %>%
    dplyr::distinct(mutation_gene, target_gene)

  log_message(sprintf("Computing network features for %d unique pairs", nrow(unique_pairs)))

  # Compute features for each unique pair
  features_list <- vector("list", nrow(unique_pairs))
  for (i in seq_len(nrow(unique_pairs))) {
    mg <- unique_pairs$mutation_gene[i]
    tg <- unique_pairs$target_gene[i]
    features_list[[i]] <- cbind(
      data.frame(mutation_gene = mg, target_gene = tg, stringsAsFactors = FALSE),
      compute_network_features_pair(graph, mg, tg, config)
    )
    if (i %% 100 == 0) {
      log_message(sprintf("  Progress: %d / %d pairs", i, nrow(unique_pairs)))
    }
  }

  network_features <- do.call(rbind, features_list)

  # Join back to candidate_df
  enriched <- candidate_df %>%
    dplyr::left_join(network_features, by = c("mutation_gene", "target_gene"))

  # Fill NA network scores with 0 (Bug 15)
  enriched <- replace_na_scores(enriched, c("network_score", "distance_score",
                                            "path_confidence", "hub_penalty"))
  enriched$connected <- ifelse(is.na(enriched$connected), FALSE, enriched$connected)
  enriched$common_neighbors <- ifelse(is.na(enriched$common_neighbors), 0L, enriched$common_neighbors)

  n_connected <- sum(enriched$connected)
  log_message(sprintf("Network features: %d / %d pairs connected in STRING",
                      n_connected, nrow(enriched)))

  # Save
  ensure_dir(config$paths$processed_dir)
  write.csv(enriched, config$output_files$network_features, row.names = FALSE)

  enriched
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  source("R/01_prepare_data.R")
  source("R/02_association_analysis.R")
  source("R/03_build_networks.R")
  data_list <- prepare_data("config/config.yaml")
  candidate_df <- run_association_pipeline(data_list)
  network_objects <- build_all_networks(data_list)
  candidate_df <- compute_network_features_batch(
    network_objects$string_graph, candidate_df, data_list$config
  )
  print(head(candidate_df))
}
