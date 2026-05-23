# =============================================================================
# 00_utils.R — Shared utility functions
# =============================================================================

# -- Scaling ----------------------------------------------------------------

#' Scale vector to [0, 1] range
#' Handles all-NA, all-identical, and empty vectors safely (Bug 7)
scale_0_1 <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# -- Pair keys --------------------------------------------------------------

#' Create an unordered pair key: sorts genes alphabetically, joins with "__"
make_pair_key <- function(geneA, geneB) {
  mapply(function(a, b) {
    paste(sort(c(a, b)), collapse = "__")
  }, geneA, geneB, USE.NAMES = FALSE)
}

#' Parse a pair key back into geneA, geneB
parse_pair_key <- function(key) {
  parts <- strsplit(key, "__", fixed = TRUE)
  do.call(rbind, lapply(parts, function(x) {
    data.frame(geneA = x[1], geneB = x[2], stringsAsFactors = FALSE)
  }))
}

# -- Safe shortest paths ----------------------------------------------------

#' Safe shortest path on undirected graph (Bug 1, Bug 4)
#' Returns list with path info or NULL-compatible disconnected result
safe_shortest_path <- function(graph, from, to, weights = NULL) {
  if (is.null(graph)) return(NULL)
  vnames <- igraph::V(graph)$name
  if (!(from %in% vnames) || !(to %in% vnames)) return(NULL)

  sp <- tryCatch({
    igraph::shortest_paths(graph, from = from, to = to,
                           weights = weights, output = "both")
  }, error = function(e) NULL)

  if (is.null(sp) || length(sp$vpath[[1]]) == 0) return(NULL)

  path_vertices <- sp$vpath[[1]]$name
  if (length(path_vertices) < 2) return(NULL)

  path_edges <- sp$epath[[1]]
  list(
    path_genes   = path_vertices,
    distance     = length(path_vertices) - 1,
    edge_confs   = if (is.null(path_edges)) numeric(0) else path_edges$edge_confidence,
    edge_weights = if (is.null(path_edges)) numeric(0) else path_edges$distance_weight
  )
}

#' Safe directed shortest path
safe_shortest_path_directed <- function(graph, from, to, weights = NULL) {
  if (is.null(graph)) return(NULL)
  vnames <- igraph::V(graph)$name
  if (!(from %in% vnames) || !(to %in% vnames)) return(NULL)

  sp <- tryCatch({
    igraph::shortest_paths(graph, from = from, to = to,
                           weights = weights, output = "both", mode = "out")
  }, error = function(e) NULL)

  if (is.null(sp) || length(sp$vpath[[1]]) == 0) return(NULL)

  path_vertices <- sp$vpath[[1]]$name
  if (length(path_vertices) < 2) return(NULL)

  path_edges <- sp$epath[[1]]
  list(
    path_genes   = path_vertices,
    distance     = length(path_vertices) - 1,
    edge_confs   = if (is.null(path_edges)) numeric(0) else path_edges$edge_confidence,
    edge_weights = if (is.null(path_edges)) numeric(0) else path_edges$distance_weight
  )
}

# -- Gene symbol standardization (Bug 1, Bug 2) -----------------------------

#' Standardize gene symbols: trim, uppercase, remove empty
standardize_gene_symbols <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  x[x == ""] <- NA_character_
  x
}

# -- Logging ----------------------------------------------------------------

#' Timestamped log message
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
}

# -- Gene overlap checking (Bug 1) ------------------------------------------

#' Report gene overlap between query and network
check_gene_overlap <- function(query_genes, network_genes, label = "") {
  query_genes <- unique(query_genes)
  total <- length(query_genes)
  found <- sum(query_genes %in% network_genes)
  coverage <- if (total > 0) round(100 * found / total, 1) else 0
  log_message(sprintf("Gene coverage [%s]: %d / %d (%.1f%%)",
                      label, found, total, coverage))
  if (coverage < 70 && total > 10) {
    log_message(sprintf(
      "WARNING: Low gene coverage (%.1f%%) for [%s]. Check gene symbol standardization.",
      coverage, label), "WARN")
  }
  invisible(coverage)
}

# -- Safe log10 for FDR (Bug 6) ---------------------------------------------

safe_log10 <- function(x, pseudocount = 1e-300) {
  -log10(x + pseudocount)
}

# -- NA replacement (Bug 15) ------------------------------------------------

replace_na_scores <- function(df, score_cols) {
  for (col in score_cols) {
    if (col %in% names(df)) {
      df[[col]] <- ifelse(is.na(df[[col]]), 0, df[[col]])
    }
  }
  df
}

# -- Directory helper -------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    log_message(paste("Created directory:", path))
  }
}

# -- Detect STRING protein IDs (Bug 2) --------------------------------------

#' Check if a vector looks like STRING protein IDs (e.g., 9606.ENSP...)
is_protein_id <- function(x) {
  grepl("^9606\\.ENSP", x)
}

#' Map STRING protein IDs to gene symbols using an alias file
#' Expects alias_df with columns: protein_id, gene_symbol
map_protein_to_gene <- function(ids, alias_df) {
  if (!is.data.frame(alias_df)) return(ids)
  lookup <- setNames(alias_df[[2]], alias_df[[1]])
  mapped <- lookup[ids]
  ifelse(is.na(mapped), ids, mapped)
}

# -- Regression helper (Bug 8) ----------------------------------------------

#' Safely count non-NA mutation samples
count_mut_samples <- function(mut_vec) {
  sum(!is.na(mut_vec) & mut_vec == 1)
}

# =============================================================================
# Log startup
# =============================================================================
log_message("Utilities (00_utils.R) loaded successfully")
