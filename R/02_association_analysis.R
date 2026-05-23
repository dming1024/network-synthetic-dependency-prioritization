# =============================================================================
# 02_association_analysis.R — Regression-based mutation–dependency association
# =============================================================================

#' Build regression formula string from config
build_regression_formula <- function(config) {
  covariates <- config$association$covariates
  cov_str <- paste(covariates, collapse = " + ")
  if (config$association$use_expression_covariate) {
    formula_str <- sprintf("dependency ~ mutation + target_expression + %s", cov_str)
  } else {
    formula_str <- sprintf("dependency ~ mutation + %s", cov_str)
  }
  formula_str
}

#' Run a single linear regression for one (mutation_gene, dependency_gene) pair
#' Handles singular model matrix, perfect separation, small sample size (Bug 8, Bug 9)
run_single_regression <- function(dep_vec, mut_vec, expr_vec, metadata, formula_str, config) {
  # Count mutant samples
  n_mut <- sum(!is.na(mut_vec) & mut_vec == 1)
  n_wt  <- sum(!is.na(mut_vec) & mut_vec == 0)

  min_mut <- config$association$min_mut_count
  if (n_mut < min_mut) {
    return(data.frame(
      beta = NA_real_, p_value = NA_real_, delta_dependency = NA_real_,
      n_mut = n_mut, n_wt = n_wt, converged = FALSE,
      error_message = sprintf("n_mut (%d) < min_mut_count (%d)", n_mut, min_mut),
      stringsAsFactors = FALSE
    ))
  }

  # Build model data
  df <- data.frame(
    dependency = dep_vec,
    mutation   = mut_vec,
    stringsAsFactors = FALSE
  )

  if (config$association$use_expression_covariate && !is.null(expr_vec)) {
    df$target_expression <- expr_vec
  }

  # Add metadata covariates
  if (!is.null(metadata)) {
    for (cov in config$association$covariates) {
      if (cov %in% names(metadata)) {
        df[[cov]] <- metadata[[cov]]
      }
    }
  }

  # Remove rows with NA in any variable
  df <- df[complete.cases(df), ]

  # Check remaining mutant samples
  n_mut_eff <- sum(df$mutation == 1)
  if (n_mut_eff < min_mut) {
    return(data.frame(
      beta = NA_real_, p_value = NA_real_, delta_dependency = NA_real_,
      n_mut = n_mut, n_wt = n_wt, converged = FALSE,
      error_message = sprintf("Effective n_mut (%d) < min after NA removal", n_mut_eff),
      stringsAsFactors = FALSE
    ))
  }

  # Fit model
  fit <- tryCatch({
    lm(as.formula(formula_str), data = df)
  }, error = function(e) {
    return(list(error = e$message))
  })

  if (!inherits(fit, "lm")) {
    return(data.frame(
      beta = NA_real_, p_value = NA_real_, delta_dependency = NA_real_,
      n_mut = n_mut, n_wt = n_wt, converged = FALSE,
      error_message = fit$error,
      stringsAsFactors = FALSE
    ))
  }

  # Extract coefficients
  coef_summary <- summary(fit)$coefficients
  if (!"mutation" %in% rownames(coef_summary)) {
    return(data.frame(
      beta = NA_real_, p_value = NA_real_, delta_dependency = NA_real_,
      n_mut = n_mut, n_wt = n_wt, converged = FALSE,
      error_message = "mutation coefficient not in model (singular model matrix)",
      stringsAsFactors = FALSE
    ))
  }

  beta    <- coef_summary["mutation", "Estimate"]
  p_value <- coef_summary["mutation", "Pr(>|t|)"]

  # Compute delta dependency
  dep_mut <- df$dependency[df$mutation == 1]
  dep_wt  <- df$dependency[df$mutation == 0]
  delta   <- mean(dep_mut, na.rm = TRUE) - mean(dep_wt, na.rm = TRUE)

  data.frame(
    beta             = beta,
    p_value          = p_value,
    delta_dependency = delta,
    n_mut            = n_mut,
    n_wt             = n_wt,
    converged        = TRUE,
    error_message    = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Run mutation-dependency association for all gene pairs (Bug 20: chunked processing)
run_mutation_dependency_association <- function(data_list) {
  log_message("=== Step 2: Association Analysis ===")

  config <- data_list$config
  id_col <- config$data_prep$sample_id_column

  mut_mat <- data_list$mutation_mat
  dep_mat <- data_list$dependency_mat
  expr_mat <- data_list$expression_mat
  meta    <- data_list$metadata

  if (is.null(mut_mat) || is.null(dep_mat)) {
    stop("Mutation and dependency matrices are required for association analysis")
  }

  # Get gene lists
  mut_genes <- setdiff(names(mut_mat), id_col)
  dep_genes <- setdiff(names(dep_mat), id_col)
  log_message(sprintf("Testing %d mutation genes x %d dependency genes = %d total pairs",
                      length(mut_genes), length(dep_genes),
                      length(mut_genes) * length(dep_genes)))

  # Align samples
  common_samples <- intersect(mut_mat[[id_col]], dep_mat[[id_col]])
  if (!is.null(expr_mat)) {
    common_samples <- intersect(common_samples, expr_mat[[id_col]])
  }
  if (!is.null(meta)) {
    common_samples <- intersect(common_samples, meta[[id_col]])
  }
  log_message(sprintf("Common samples: %d", length(common_samples)))

  # Subset to common samples
  mut_idx <- match(common_samples, mut_mat[[id_col]])
  dep_idx <- match(common_samples, dep_mat[[id_col]])
  meta_idx <- if (!is.null(meta)) match(common_samples, meta[[id_col]]) else NULL

  # Build formula
  formula_str <- build_regression_formula(config)

  # Chunked processing to manage memory (Bug 20)
  chunk_size <- config$association$chunk_size
  n_chunks <- ceiling(length(dep_genes) / chunk_size)

  results_list <- vector("list", length(mut_genes) * length(dep_genes))
  result_idx <- 0

  for (i in seq_along(mut_genes)) {
    mg <- mut_genes[i]
    mut_vec <- mut_mat[mut_idx, mg]

    for (chunk in seq_len(n_chunks)) {
      start_g <- (chunk - 1) * chunk_size + 1
      end_g   <- min(chunk * chunk_size, length(dep_genes))
      chunk_genes <- dep_genes[start_g:end_g]

      for (tg in chunk_genes) {
        dep_vec <- dep_mat[dep_idx, tg]
        expr_vec <- if (!is.null(expr_mat) && tg %in% names(expr_mat)) {
          expr_mat[match(common_samples, expr_mat[[id_col]]), tg]
        } else NULL

        result_idx <- result_idx + 1
        res <- run_single_regression(dep_vec, mut_vec, expr_vec, meta[meta_idx, , drop = FALSE],
                                     formula_str, config)
        results_list[[result_idx]] <- data.frame(
          mutation_gene = mg,
          target_gene   = tg,
          res,
          stringsAsFactors = FALSE
        )
      }
    }
    if (i %% 20 == 0) {
      log_message(sprintf("  Progress: %d / %d mutation genes", i, length(mut_genes)))
    }
  }

  assoc_df <- do.call(rbind, results_list)
  log_message(sprintf("Association results: %d pairs tested", nrow(assoc_df)))
  assoc_df
}

#' Apply FDR correction (Bug 6: safe handling of NA p-values)
apply_fdr_correction <- function(assoc_df, method = "BH") {
  valid_idx <- !is.na(assoc_df$p_value)
  assoc_df$FDR <- NA_real_
  if (sum(valid_idx) > 0) {
    assoc_df$FDR[valid_idx] <- p.adjust(assoc_df$p_value[valid_idx], method = method)
  }
  n_sig <- sum(assoc_df$FDR < 0.05, na.rm = TRUE)
  log_message(sprintf("FDR correction applied (%s). %d pairs at FDR < 0.05", method, n_sig))
  assoc_df
}

#' Filter candidates by statistical thresholds
filter_candidates <- function(assoc_df, config) {
  cf <- config$candidate_filtering
  fdr_thresh <- cf$fdr_threshold

  candidates <- assoc_df

  if (cf$beta_direction == "negative") {
    candidates <- candidates[candidates$beta < 0 & !is.na(candidates$beta), ]
  }
  candidates <- candidates[candidates$FDR < fdr_thresh & !is.na(candidates$FDR), ]
  candidates <- candidates[candidates$n_mut >= cf$min_n_mut, ]
  candidates <- candidates[abs(candidates$delta_dependency) >= cf$min_abs_delta_dependency, ]
  candidates <- candidates[!is.na(candidates$delta_dependency), ]

  log_message(sprintf("Candidates after filtering: %d pairs", nrow(candidates)))
  candidates
}

#' Run the full association pipeline
run_association_pipeline <- function(data_list) {
  config <- data_list$config

  assoc_df <- run_mutation_dependency_association(data_list)
  assoc_df <- apply_fdr_correction(assoc_df, config$association$fdr_method)

  # Filter candidates
  candidate_df <- filter_candidates(assoc_df, config)

  # Add raw stat score
  candidate_df$stat_score_raw <- with(candidate_df,
    safe_log10(FDR) * abs(beta)
  )

  # Save
  ensure_dir(config$paths$processed_dir)
  write.csv(assoc_df, config$output_files$association_results, row.names = FALSE)
  write.csv(candidate_df, config$output_files$candidate_pairs, row.names = FALSE)
  log_message(sprintf("Saved association results and %d candidates", nrow(candidate_df)))

  candidate_df
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  source("R/01_prepare_data.R")
  data_list <- prepare_data("config/config.yaml")
  candidate_df <- run_association_pipeline(data_list)
  print(head(candidate_df))
}
