# =============================================================================
# 07_null_model_permutation.R — Degree-matched empirical network p-values
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

#' Build a target gene feature table for matching
build_target_feature_table <- function(graph, dependency_mat, config) {
  if (is.null(graph)) return(NULL)

  id_col <- config$data_prep$sample_id_column

  vnames <- igraph::V(graph)$name
  degrees <- igraph::degree(graph)

  dep_vars <- NULL
  if (!is.null(dependency_mat)) {
    dep_genes <- setdiff(names(dependency_mat), id_col)
    dep_genes <- intersect(dep_genes, vnames)
    if (length(dep_genes) > 0) {
      dep_vars <- apply(dependency_mat[, dep_genes, drop = FALSE], 2, var, na.rm = TRUE)
    }
  }

  df <- data.frame(
    gene      = vnames,
    degree    = degrees,
    dep_var   = if (!is.null(dep_vars)) dep_vars[match(vnames, names(dep_vars))] else NA_real_,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  df
}

#' Sample a degree-matched target gene
sample_matched_target <- function(target_gene, target_features, match_on = "degree") {
  if (is.null(target_features)) return(NA_character_)

  tf_row <- target_features[target_features$gene == target_gene, ]
  if (nrow(tf_row) == 0) return(NA_character_)

  target_deg <- tf_row$degree[1]

  # Find genes in similar degree bin
  deg_tolerance <- max(2, ceiling(target_deg * 0.3))
  matched <- target_features[
    abs(target_features$degree - target_deg) <= deg_tolerance &
    target_features$gene != target_gene,
  ]

  if (nrow(matched) < 5) {
    # Widen tolerance
    matched <- target_features[
      target_features$gene != target_gene,
    ]
  }

  if (nrow(matched) == 0) return(NA_character_)

  sample(matched$gene, 1)
}

#' Compute empirical network p-value for a single pair (Bug 3 mitigation)
compute_empirical_network_p <- function(graph, geneA, geneB, target_features,
                                        n_perm = 1000, match_on = "degree",
                                        config) {
  if (is.null(graph)) {
    return(data.frame(
      empirical_network_p = NA_real_,
      observed_network_score = NA_real_,
      mean_random_score = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  vnames <- igraph::V(graph)$name
  if (!(geneA %in% vnames) || !(geneB %in% vnames)) {
    return(data.frame(
      empirical_network_p = NA_real_,
      observed_network_score = 0,
      mean_random_score = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  # Observed score
  obs <- compute_network_features_pair(graph, geneA, geneB, config)
  observed <- obs$network_score[1]

  # Permutation
  random_scores <- numeric(n_perm)
  set.seed(config$permutation$random_seed)

  for (i in seq_len(n_perm)) {
    B_rand <- sample_matched_target(geneB, target_features, match_on)
    if (is.na(B_rand)) {
      random_scores[i] <- 0
    } else {
      rand_res <- compute_network_features_pair(graph, geneA, B_rand, config)
      random_scores[i] <- rand_res$network_score[1]
    }
  }

  # Empirical p-value
  emp_p <- (sum(random_scores >= observed, na.rm = TRUE) + 1) / (n_perm + 1)
  mean_rand <- mean(random_scores, na.rm = TRUE)

  data.frame(
    empirical_network_p   = emp_p,
    observed_network_score = observed,
    mean_random_score     = mean_rand,
    stringsAsFactors      = FALSE
  )
}

#' Run permutation tests on top candidate pairs
run_permutation_tests <- function(ranked_df, graph, dependency_mat, config) {
  log_message("=== Step 7: Permutation Tests ===")

  if (is.null(graph)) {
    log_message("No STRING graph available, skipping permutation tests", "WARN")
    ranked_df$empirical_network_p <- NA_real_
    ranked_df$observed_network_score <- NA_real_
    ranked_df$mean_random_score <- NA_real_
    return(ranked_df)
  }

  if (!config$permutation$enabled) {
    log_message("Permutation testing disabled in config, skipping")
    ranked_df$empirical_network_p <- NA_real_
    ranked_df$observed_network_score <- NA_real_
    ranked_df$mean_random_score <- NA_real_
    return(ranked_df)
  }

  # Only test top pairs (expensive)
  top_n <- min(config$permutation$top_n_pairs, nrow(ranked_df))
  log_message(sprintf("Running permutation tests on top %d pairs (%d permutations each)",
                      top_n, config$permutation$n_perm))

  target_features <- build_target_feature_table(graph, dependency_mat, config)

  perm_results <- vector("list", top_n)
  for (i in seq_len(top_n)) {
    perm_results[[i]] <- cbind(
      data.frame(
        mutation_gene = ranked_df$mutation_gene[i],
        target_gene   = ranked_df$target_gene[i],
        stringsAsFactors = FALSE
      ),
      compute_empirical_network_p(
        graph,
        ranked_df$mutation_gene[i],
        ranked_df$target_gene[i],
        target_features,
        config$permutation$n_perm,
        config$permutation$match_on,
        config
      )
    )
    if (i %% 50 == 0) {
      log_message(sprintf("  Progress: %d / %d", i, top_n))
    }
  }

  perm_df <- do.call(rbind, perm_results)

  # Join back
  ranked_df <- ranked_df %>%
    dplyr::left_join(perm_df, by = c("mutation_gene", "target_gene"))

  # Fill NAs for untested rows
  ranked_df$empirical_network_p   <- ifelse(is.na(ranked_df$empirical_network_p), NA_real_, ranked_df$empirical_network_p)
  ranked_df$observed_network_score <- ifelse(is.na(ranked_df$observed_network_score), NA_real_, ranked_df$observed_network_score)
  ranked_df$mean_random_score     <- ifelse(is.na(ranked_df$mean_random_score), NA_real_, ranked_df$mean_random_score)

  # Save
  ensure_dir(config$paths$processed_dir)
  write.csv(ranked_df, config$output_files$permutation_results, row.names = FALSE)

  n_sig <- sum(ranked_df$empirical_network_p < 0.05, na.rm = TRUE)
  log_message(sprintf("Permutation testing complete. %d pairs with empirical p < 0.05", n_sig))

  ranked_df
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  source("R/03_build_networks.R")
  source("R/04_compute_network_features.R")
  source("R/06_score_candidates.R")
  config <- yaml::read_yaml("config/config.yaml")
  ranked_df <- read.csv(config$output_files$final_ranked, stringsAsFactors = FALSE)
  data_list <- list(config = config)
  data_list$string_edges <- read.csv(config$input_files$string_edges, stringsAsFactors = FALSE)
  data_list$dependency_mat <- read.csv(config$input_files$dependency_matrix, stringsAsFactors = FALSE, check.names = FALSE)
  network_objects <- build_all_networks(data_list)
  ranked_df <- run_permutation_tests(ranked_df, network_objects$string_graph, data_list$dependency_mat, config)
}
