# =============================================================================
# 01_prepare_data.R — Data loading, validation, and gene symbol standardization
# =============================================================================

#' Load config from YAML file
load_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop(sprintf("Config file not found: %s", config_path))
  }
  config <- yaml::read_yaml(config_path)
  log_message(sprintf("Config loaded from %s", config_path))
  config
}

#' Read CSV with error handling
read_input_csv <- function(filepath) {
  if (!file.exists(filepath)) {
    log_message(sprintf("Input file not found: %s (will be skipped)", filepath), "WARN")
    return(NULL)
  }
  df <- read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE)
  log_message(sprintf("Read %s: %d rows x %d columns", basename(filepath), nrow(df), ncol(df)))
  df
}

#' Standardize sample IDs: trim, check duplicates
standardize_sample_ids <- function(df, id_col = "sample_id") {
  if (!id_col %in% names(df)) {
    stop(sprintf("ID column '%s' not found in data.frame", id_col))
  }
  df[[id_col]] <- trimws(as.character(df[[id_col]]))
  dup_ids <- df[[id_col]][duplicated(df[[id_col]])]
  if (length(dup_ids) > 0) {
    stop(sprintf("Duplicate sample IDs found: %s", paste(unique(dup_ids), collapse = ", ")))
  }
  df
}

#' Standardize mutation matrix: ensure binary, filter dead columns (Bug 8 pre-filter)
standardize_mutation_matrix <- function(mat, id_col = "sample_id") {
  if (is.null(mat)) return(NULL)

  ids <- mat[[id_col]]
  gene_cols <- setdiff(names(mat), id_col)

  for (col in gene_cols) {
    vals <- mat[[col]]
    if (!all(vals %in% c(0, 1, NA))) {
      log_message(sprintf("Non-binary values in mutation column '%s', coercing to 0/1", col), "WARN")
      mat[[col]] <- ifelse(!is.na(vals) & vals != 0, 1, 0)
    }
  }

  # Remove genes with zero mutations
  mut_counts <- colSums(mat[, gene_cols, drop = FALSE], na.rm = TRUE)
  dead_genes <- names(mut_counts[mut_counts == 0])
  if (length(dead_genes) > 0) {
    log_message(sprintf("Removing %d genes with zero mutations", length(dead_genes)), "WARN")
    mat <- mat[, !names(mat) %in% dead_genes, drop = FALSE]
  }

  log_message(sprintf("Mutation matrix: %d samples x %d genes", nrow(mat), ncol(mat) - 1))
  mat
}

#' Standardize dependency matrix: remove near-zero-variance genes (Bug 18)
standardize_dependency_matrix <- function(mat, id_col = "sample_id", min_var = 0.01) {
  if (is.null(mat)) return(NULL)

  ids <- mat[[id_col]]
  gene_cols <- setdiff(names(mat), id_col)

  gene_vars <- apply(mat[, gene_cols, drop = FALSE], 2, var, na.rm = TRUE)
  low_var <- names(gene_vars[gene_vars < min_var])
  if (length(low_var) > 0) {
    log_message(sprintf("Removing %d low-variance dependency genes", length(low_var)), "WARN")
    mat <- mat[, !names(mat) %in% low_var, drop = FALSE]
  }

  # Flag pan-essential genes (mean dependency < -0.5) (Bug 18)
  gene_cols <- setdiff(names(mat), id_col)
  gene_means <- colMeans(mat[, gene_cols, drop = FALSE], na.rm = TRUE)
  pan_essential <- names(gene_means[gene_means < -0.5])
  if (length(pan_essential) > 0) {
    log_message(sprintf(
      "Note: %d genes flagged as potentially pan-essential (mean dependency < -0.5)",
      length(pan_essential)), "INFO")
  }

  log_message(sprintf("Dependency matrix: %d samples x %d genes", nrow(mat), ncol(mat) - 1))
  mat
}

#' Harmonize gene symbols across all data sources (Bug 1)
harmonize_gene_symbols_across_data <- function(data_list) {
  log_message("Harmonizing gene symbols...")

  # Collect all gene sets
  gene_sets <- list()

  if (!is.null(data_list$mutation_mat)) {
    gene_sets$mutation <- setdiff(names(data_list$mutation_mat),
                                  data_list$config$data_prep$sample_id_column)
  }
  if (!is.null(data_list$dependency_mat)) {
    gene_sets$dependency <- setdiff(names(data_list$dependency_mat),
                                    data_list$config$data_prep$sample_id_column)
  }
  if (!is.null(data_list$string_edges)) {
    gene_sets$STRING <- unique(c(data_list$string_edges$geneA, data_list$string_edges$geneB))
  }
  if (!is.null(data_list$signaling_edges)) {
    gene_sets$signaling <- unique(c(data_list$signaling_edges$source_gene,
                                    data_list$signaling_edges$target_gene))
  }

  # Standardize each set
  for (nm in names(gene_sets)) {
    gene_sets[[nm]] <- standardize_gene_symbols(gene_sets[[nm]])
    gene_sets[[nm]] <- gene_sets[[nm]][!is.na(gene_sets[[nm]])]
  }

  # Report coverage
  if (length(gene_sets$STRING) > 0) {
    check_gene_overlap(gene_sets$mutation, gene_sets$STRING, "mutation genes in STRING")
    check_gene_overlap(gene_sets$dependency, gene_sets$STRING, "dependency genes in STRING")
  }

  log_message("Gene symbol harmonization complete")
  invisible(gene_sets)
}

#' Prepare data: the main orchestrator
prepare_data <- function(config_path) {
  log_message("=== Step 1: Prepare Data ===")

  config <- load_config(config_path)
  pf      <- config$input_files
  id_col  <- config$data_prep$sample_id_column

  # Load raw data
  log_message("Loading input files...")
  data_list <- list(
    config          = config,
    mutation_mat    = NULL,
    dependency_mat  = NULL,
    expression_mat  = NULL,
    metadata        = NULL,
    string_edges    = NULL,
    reactome_df     = NULL,
    signaling_edges = NULL,
    corum_df        = NULL,
    sl_db           = NULL,
    drug_df         = NULL
  )

  data_list$metadata        <- read_input_csv(pf$sample_metadata)
  data_list$mutation_mat    <- read_input_csv(pf$mutation_matrix)
  data_list$dependency_mat  <- read_input_csv(pf$dependency_matrix)
  data_list$expression_mat  <- read_input_csv(pf$expression_matrix)
  data_list$string_edges    <- read_input_csv(pf$string_edges)
  data_list$reactome_df     <- read_input_csv(pf$reactome_pathway)
  data_list$signaling_edges <- read_input_csv(pf$omnipath_edges)
  data_list$corum_df        <- read_input_csv(pf$corum_complexes)
  data_list$sl_db           <- read_input_csv(pf$synlethdb_pairs)
  data_list$drug_df         <- read_input_csv(pf$druggable_targets)

  # Standardize sample IDs
  if (!is.null(data_list$metadata)) {
    data_list$metadata <- standardize_sample_ids(data_list$metadata, id_col)
  }
  if (!is.null(data_list$mutation_mat)) {
    data_list$mutation_mat <- standardize_sample_ids(data_list$mutation_mat, id_col)
    data_list$mutation_mat <- standardize_mutation_matrix(data_list$mutation_mat, id_col)
  }
  if (!is.null(data_list$dependency_mat)) {
    data_list$dependency_mat <- standardize_sample_ids(data_list$dependency_mat, id_col)
    data_list$dependency_mat <- standardize_dependency_matrix(
      data_list$dependency_mat, id_col,
      config$data_prep$min_dependency_variance
    )
  }
  if (!is.null(data_list$expression_mat)) {
    data_list$expression_mat <- standardize_sample_ids(data_list$expression_mat, id_col)
  }

  # Standardize gene symbols in network files
  if (config$data_prep$standardize_genes) {
    if (!is.null(data_list$string_edges)) {
      data_list$string_edges$geneA <- standardize_gene_symbols(data_list$string_edges$geneA)
      data_list$string_edges$geneB <- standardize_gene_symbols(data_list$string_edges$geneB)
    }
    if (!is.null(data_list$reactome_df)) {
      data_list$reactome_df$gene <- standardize_gene_symbols(data_list$reactome_df$gene)
    }
    if (!is.null(data_list$signaling_edges)) {
      data_list$signaling_edges$source_gene <- standardize_gene_symbols(data_list$signaling_edges$source_gene)
      data_list$signaling_edges$target_gene <- standardize_gene_symbols(data_list$signaling_edges$target_gene)
    }
    if (!is.null(data_list$corum_df)) {
      data_list$corum_df$gene <- standardize_gene_symbols(data_list$corum_df$gene)
    }
    if (!is.null(data_list$sl_db)) {
      data_list$sl_db$geneA <- standardize_gene_symbols(data_list$sl_db$geneA)
      data_list$sl_db$geneB <- standardize_gene_symbols(data_list$sl_db$geneB)
    }
    if (!is.null(data_list$drug_df)) {
      data_list$drug_df$gene <- standardize_gene_symbols(data_list$drug_df$gene)
    }
  }

  # Report gene coverage across networks
  harmonize_gene_symbols_across_data(data_list)

  # Create output directories
  ensure_dir(config$paths$processed_dir)

  log_message("=== Data preparation complete ===")
  data_list
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  data_list <- prepare_data("config/config.yaml")
  log_message("Data list names:")
  log_message(paste(names(data_list), collapse = ", "))
}
