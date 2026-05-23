# =============================================================================
# 06_score_candidates.R — Final weighted scoring and evidence classification
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

#' Compute scaled statistical score (Bug 6: safe_log10, Bug 7: scale_0_1 safety)
compute_stat_score <- function(feature_df) {
  stat_score_raw <- safe_log10(feature_df$FDR) * abs(feature_df$beta)
  scale_0_1(stat_score_raw)
}

#' Compute all three final score variants
compute_final_scores <- function(feature_df, config) {
  sc <- config$scoring
  scheme <- sc$default_scheme
  weights <- sc[[scheme]]

  # Ensure all component scores are scaled 0-1
  feature_df$stat_score <- compute_stat_score(feature_df)
  feature_df$network_score_scaled <- scale_0_1(feature_df$network_score)

  # Compute each score variant
  compute_weighted <- function(df, w) {
    w$stat_weight      * df$stat_score +
    w$network_weight   * df$network_score_scaled +
    w$pathway_weight   * df$pathway_score +
    w$signaling_weight * df$signaling_score +
    w$complex_weight   * df$complex_score +
    w$sl_prior_weight  * df$sl_prior_score +
    w$drug_weight      * df$drug_score
  }

  feature_df$final_SL_score         <- compute_weighted(feature_df, sc$balanced)
  feature_df$final_mechanism_score  <- compute_weighted(feature_df, sc$mechanism)
  feature_df$final_therapeutic_score <- compute_weighted(feature_df, sc$therapeutic)

  feature_df
}

#' Classify candidates into Class I-IV (Bug 11: Class IV preserves statistical-only)
classify_candidates <- function(feature_df, config) {
  cl <- config$classification

  feature_df <- feature_df %>%
    dplyr::mutate(
      evidence_class = dplyr::case_when(

        # Class I: high-confidence
        (!is.na(FDR) & FDR < cl$class_I$fdr_max) &
        (!is.na(beta) & beta < cl$class_I$beta_max) &
        (!is.na(shortest_distance) & shortest_distance <= cl$class_I$network_distance_max) &
        (!is.na(pathway_score) & pathway_score >= cl$class_I$pathway_score_min) &
        (!is.na(final_SL_score) & final_SL_score >= cl$class_I$final_score_min) ~
          cl$class_I$label,

        # Class II: network-supported
        (!is.na(FDR) & FDR < cl$class_II$fdr_max) &
        (!is.na(beta) & beta < cl$class_II$beta_max) &
        (!is.na(shortest_distance) & shortest_distance <= cl$class_II$network_distance_max) &
        (!is.na(final_SL_score) & final_SL_score >= cl$class_II$final_score_min) ~
          cl$class_II$label,

        # Class III: SL-prior-supported
        (!is.na(FDR) & FDR < cl$class_III$fdr_max) &
        (!is.na(beta) & beta < cl$class_III$beta_max) &
        (!is.na(sl_prior_score) & sl_prior_score > cl$class_III$sl_prior_score_min) ~
          cl$class_III$label,

        # Class IV: statistical-only (Bug 11)
        (!is.na(FDR) & FDR < cl$class_IV$fdr_max) &
        (!is.na(beta) & beta < cl$class_IV$beta_max) &
        (!is.na(network_score_scaled) & network_score_scaled <= cl$class_IV$network_score_scaled_max) ~
          cl$class_IV$label,

        TRUE ~ cl$low_priority$label
      )
    )

  # Report class distribution
  class_counts <- table(feature_df$evidence_class)
  for (cname in names(class_counts)) {
    log_message(sprintf("  %s: %d", cname, class_counts[cname]))
  }

  feature_df
}

#' Score, rank, and classify candidates
score_and_rank <- function(feature_df, config) {
  log_message("=== Step 6: Score Candidates ===")

  if (nrow(feature_df) == 0) {
    log_message("No candidates to score", "WARN")
    return(feature_df)
  }

  # Compute scores
  feature_df <- compute_final_scores(feature_df, config)

  # Classify
  feature_df <- classify_candidates(feature_df, config)

  # Rank by final score
  feature_df <- feature_df %>%
    dplyr::arrange(dplyr::desc(final_SL_score)) %>%
    dplyr::mutate(rank = dplyr::row_number())

  # Select and order output columns
  output_cols <- c(
    "rank", "mutation_gene", "target_gene",
    "beta", "FDR", "delta_dependency", "n_mut",
    "stat_score", "stat_score_raw",
    "shortest_distance", "path_genes", "path_confidence",
    "hub_penalty", "distance_score", "network_score", "network_score_scaled",
    "pathway_score", "shared_pathways",
    "signaling_score", "signaling_direction", "signaling_path",
    "complex_score", "complex_evidence",
    "sl_prior_score", "sl_prior_type",
    "codep_r", "codependency_score",
    "drug_score", "drug_category",
    "final_SL_score", "final_mechanism_score", "final_therapeutic_score",
    "evidence_class"
  )
  output_cols <- intersect(output_cols, names(feature_df))
  ranked_df <- feature_df[, output_cols, drop = FALSE]

  # Save
  ensure_dir(config$paths$processed_dir)
  write.csv(ranked_df, config$output_files$final_ranked, row.names = FALSE)

  log_message(sprintf("Ranked %d candidates saved", nrow(ranked_df)))
  log_message("=== Scoring complete ===")
  ranked_df
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  # Read previously saved features
  config <- yaml::read_yaml("config/config.yaml")
  feature_df <- read.csv(config$output_files$empirical_features, stringsAsFactors = FALSE)
  ranked_df <- score_and_rank(feature_df, config)
  print(head(ranked_df[, c("rank", "mutation_gene", "target_gene", "final_SL_score", "evidence_class")]))
}
