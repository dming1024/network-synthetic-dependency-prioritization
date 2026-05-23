# =============================================================================
# 05_compute_empirical_evidence.R — Pathway, signaling, complex, SL prior, drug scores
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

#' Compute pathway concordance score
compute_pathway_score <- function(geneA, geneB, reactome_index, config) {
  default <- data.frame(
    pathway_score = 0, shared_pathways = NA_character_, stringsAsFactors = FALSE
  )
  if (is.null(reactome_index)) return(default)

  g2p <- reactome_index$gene_to_pathways
  pw_A <- g2p[[geneA]]
  pw_B <- g2p[[geneB]]

  if (is.null(pw_A) || is.null(pw_B)) return(default)

  shared <- intersect(pw_A, pw_B)
  if (length(shared) > 0) {
    return(data.frame(
      pathway_score   = config$empirical_evidence$pathway$exact_match_score,
      shared_pathways = paste(shared, collapse = "; "),
      stringsAsFactors = FALSE
    ))
  }

  # Parent pathway heuristic: check if any pathway shares a common prefix (parent)
  # Simple approach: extract parent from pathway name by splitting on common delimiters
  parents_A <- unique(gsub("_[0-9]+$|\\.[0-9]+$", "", pw_A))
  parents_B <- unique(gsub("_[0-9]+$|\\.[0-9]+$", "", pw_B))
  shared_parents <- intersect(parents_A, parents_B)

  if (length(shared_parents) > 0) {
    return(data.frame(
      pathway_score   = config$empirical_evidence$pathway$parent_match_score,
      shared_pathways = paste("Parent:", paste(shared_parents, collapse = "; ")),
      stringsAsFactors = FALSE
    ))
  }

  default
}

#' Compute directional signaling evidence score (Bug 13: cautious interpretation)
compute_signaling_score <- function(signaling_graph, geneA, geneB, config) {
  default <- data.frame(
    signaling_score     = 0,
    signaling_direction = "none",
    signaling_path      = NA_character_,
    stringsAsFactors    = FALSE
  )
  if (is.null(signaling_graph)) return(default)

  ec <- config$empirical_evidence$signaling
  max_len <- ec$max_directed_length

  # Try A -> B
  sp_AB <- safe_shortest_path_directed(signaling_graph, geneA, geneB,
                                       weights = igraph::E(signaling_graph)$distance_weight)
  # Try B -> A
  sp_BA <- safe_shortest_path_directed(signaling_graph, geneB, geneA,
                                       weights = igraph::E(signaling_graph)$distance_weight)

  # A -> B exists and is short
  if (!is.null(sp_AB) && sp_AB$distance <= ec$direct_short_max_len) {
    return(data.frame(
      signaling_score     = ec$direct_short_score,
      signaling_direction = "mutation_gene_to_target_gene",
      signaling_path      = paste(sp_AB$path_genes, collapse = " -> "),
      stringsAsFactors    = FALSE
    ))
  }

  # B -> A exists and is short
  if (!is.null(sp_BA) && sp_BA$distance <= ec$direct_short_max_len) {
    return(data.frame(
      signaling_score     = ec$reverse_short_score,
      signaling_direction = "target_gene_to_mutation_gene",
      signaling_path      = paste(sp_BA$path_genes, collapse = " -> "),
      stringsAsFactors    = FALSE
    ))
  }

  # Any directed path within max length
  if (!is.null(sp_AB) && sp_AB$distance <= max_len) {
    return(data.frame(
      signaling_score     = ec$any_path_score,
      signaling_direction = "weak_directional_support",
      signaling_path      = paste(sp_AB$path_genes, collapse = " -> "),
      stringsAsFactors    = FALSE
    ))
  }
  if (!is.null(sp_BA) && sp_BA$distance <= max_len) {
    return(data.frame(
      signaling_score     = ec$any_path_score,
      signaling_direction = "weak_directional_support",
      signaling_path      = paste(sp_BA$path_genes, collapse = " -> "),
      stringsAsFactors    = FALSE
    ))
  }

  default
}

#' Compute protein complex evidence score
compute_complex_score <- function(geneA, geneB, complex_index, config) {
  default <- data.frame(
    complex_score    = 0,
    complex_evidence = NA_character_,
    stringsAsFactors = FALSE
  )
  if (is.null(complex_index)) return(default)

  g2c <- complex_index$gene_to_complexes
  cA <- g2c[[geneA]]
  cB <- g2c[[geneB]]

  if (is.null(cA) || is.null(cB)) return(default)

  shared <- intersect(cA, cB)
  if (length(shared) > 0) {
    return(data.frame(
      complex_score    = config$empirical_evidence$complex$same_complex_score,
      complex_evidence = paste(shared, collapse = "; "),
      stringsAsFactors = FALSE
    ))
  }

  default
}

#' Compute synthetic lethality prior evidence score (Bug 12: separate from external validation)
compute_sl_prior_score <- function(geneA, geneB, sl_pair_set, string_graph, reactome_index, config) {
  default <- data.frame(
    sl_prior_score = 0,
    sl_prior_type  = "none",
    stringsAsFactors = FALSE
  )
  if (is.null(sl_pair_set)) return(default)

  ec <- config$empirical_evidence$sl_prior

  # Tier 1: known pair
  pair_key <- make_pair_key(geneA, geneB)
  if (pair_key %in% sl_pair_set) {
    return(data.frame(
      sl_prior_score = ec$known_pair_score,
      sl_prior_type  = "known_SL_pair",
      stringsAsFactors = FALSE
    ))
  }

  # Tier 2: neighbor-supported (one STRING neighbor is SL with the other gene)
  if (!is.null(string_graph)) {
    vnames <- igraph::V(string_graph)$name
    depth <- ec$neighbor_search_depth

    if (geneA %in% vnames) {
      nA <- igraph::ego(string_graph, order = depth, nodes = geneA)[[1]]$name
      nA <- setdiff(nA, geneA)
      for (n in nA) {
        if (make_pair_key(n, geneB) %in% sl_pair_set) {
          return(data.frame(
            sl_prior_score = ec$neighbor_supported_score,
            sl_prior_type  = "neighbor_supported_SL",
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    if (geneB %in% vnames) {
      nB <- igraph::ego(string_graph, order = depth, nodes = geneB)[[1]]$name
      nB <- setdiff(nB, geneB)
      for (n in nB) {
        if (make_pair_key(geneA, n) %in% sl_pair_set) {
          return(data.frame(
            sl_prior_score = ec$neighbor_supported_score,
            sl_prior_type  = "neighbor_supported_SL",
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  # Tier 3: pathway-level support
  if (!is.null(reactome_index)) {
    g2p <- reactome_index$gene_to_pathways
    pw_A <- g2p[[geneA]]
    pw_B <- g2p[[geneB]]
    if (!is.null(pw_A) && !is.null(pw_B) && length(intersect(pw_A, pw_B)) > 0) {
      return(data.frame(
        sl_prior_score = ec$pathway_level_score,
        sl_prior_type  = "pathway_level_SL_support",
        stringsAsFactors = FALSE
      ))
    }
  }

  default
}

#' Compute DepMap co-dependency score (Bug 18 note: co-dependency != SL)
compute_codependency_score <- function(geneA, geneB, dependency_mat, config) {
  default <- data.frame(
    codep_r           = NA_real_,
    codependency_score = 0,
    stringsAsFactors   = FALSE
  )

  if (is.null(dependency_mat)) return(default)

  id_col <- config$data_prep$sample_id_column
  dep_genes <- setdiff(names(dependency_mat), id_col)

  if (!(geneA %in% dep_genes) || !(geneB %in% dep_genes)) return(default)

  r <- tryCatch({
    cor(dependency_mat[[geneA]], dependency_mat[[geneB]],
        use = "pairwise.complete.obs")
  }, error = function(e) NA_real_)

  data.frame(
    codep_r            = r,
    codependency_score = ifelse(is.na(r), 0, abs(r)),
    stringsAsFactors   = FALSE
  )
}

#' Compute druggability score for a target gene (Bug 14: small weight in final score)
compute_druggability_score <- function(target_gene, drug_index, config) {
  default <- data.frame(
    drug_score    = 0,
    drug_category = "unknown",
    stringsAsFactors = FALSE
  )
  if (is.null(drug_index)) return(default)

  entry <- drug_index[[target_gene]]
  if (is.null(entry)) return(default)

  ec <- config$empirical_evidence$drug
  cat <- entry$drug_category[1]

  score <- switch(cat,
    FDA_approved     = ec$fda_approved_score,
    clinical         = ec$clinical_target_score,
    preclinical      = ec$preclinical_target_score,
    druggable_class  = ec$druggable_class_score,
    unknown          = ec$unknown_score,
    ec$unknown_score
  )

  data.frame(
    drug_score    = score,
    drug_category = cat,
    stringsAsFactors = FALSE
  )
}

#' Compute all empirical evidence features for all candidate pairs
compute_all_empirical_features <- function(candidate_df, network_objects, data_list, config) {
  log_message("=== Step 5: Compute Empirical Evidence ===")

  if (nrow(candidate_df) == 0) {
    log_message("No candidates to compute empirical features for", "WARN")
    return(candidate_df)
  }

  # Get unique pairs
  unique_pairs <- candidate_df %>%
    dplyr::distinct(mutation_gene, target_gene)

  log_message(sprintf("Computing empirical features for %d unique pairs", nrow(unique_pairs)))

  # Pre-compute drug features for all unique target genes (more efficient)
  unique_targets <- unique(unique_pairs$target_gene)
  drug_features <- do.call(rbind, lapply(unique_targets, function(tg) {
    cbind(data.frame(target_gene = tg, stringsAsFactors = FALSE),
          compute_druggability_score(tg, network_objects$drug_index, config))
  }))

  # Compute pair-level features
  features_list <- vector("list", nrow(unique_pairs))
  for (i in seq_len(nrow(unique_pairs))) {
    mg <- unique_pairs$mutation_gene[i]
    tg <- unique_pairs$target_gene[i]

    pw  <- compute_pathway_score(mg, tg, network_objects$reactome_index, config)
    sig <- compute_signaling_score(network_objects$signaling_graph, mg, tg, config)
    cx  <- compute_complex_score(mg, tg, network_objects$complex_index, config)
    sl  <- compute_sl_prior_score(mg, tg, network_objects$sl_pair_set,
                                  network_objects$string_graph,
                                  network_objects$reactome_index, config)
    cod <- compute_codependency_score(mg, tg, data_list$dependency_mat, config)

    features_list[[i]] <- data.frame(
      mutation_gene = mg, target_gene = tg,
      pw, sig, cx, sl, cod,
      stringsAsFactors = FALSE
    )

    if (i %% 100 == 0) {
      log_message(sprintf("  Progress: %d / %d pairs", i, nrow(unique_pairs)))
    }
  }

  empirical_features <- do.call(rbind, features_list)

  # Join drug features per target gene
  empirical_features <- empirical_features %>%
    dplyr::left_join(drug_features, by = "target_gene")

  # Merge back to candidate_df
  enriched <- candidate_df %>%
    dplyr::left_join(empirical_features, by = c("mutation_gene", "target_gene"))

  # Replace NA scores with 0 (Bug 15)
  score_cols <- c("pathway_score", "signaling_score", "complex_score",
                  "sl_prior_score", "codependency_score", "drug_score")
  enriched <- replace_na_scores(enriched, score_cols)
  enriched$shared_pathways  <- ifelse(is.na(enriched$shared_pathways), NA_character_, enriched$shared_pathways)
  enriched$signaling_path   <- ifelse(is.na(enriched$signaling_path), NA_character_, enriched$signaling_path)
  enriched$complex_evidence <- ifelse(is.na(enriched$complex_evidence), NA_character_, enriched$complex_evidence)
  enriched$drug_category    <- ifelse(is.na(enriched$drug_category), "unknown", enriched$drug_category)
  enriched$sl_prior_type    <- ifelse(is.na(enriched$sl_prior_type), "none", enriched$sl_prior_type)

  # Save
  ensure_dir(config$paths$processed_dir)
  write.csv(enriched, config$output_files$empirical_features, row.names = FALSE)

  log_message("=== Empirical evidence computation complete ===")
  enriched
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  source("R/01_prepare_data.R")
  source("R/02_association_analysis.R")
  source("R/03_build_networks.R")
  source("R/04_compute_network_features.R")
  data_list <- prepare_data("config/config.yaml")
  candidate_df <- run_association_pipeline(data_list)
  network_objects <- build_all_networks(data_list)
  candidate_df <- compute_network_features_batch(network_objects$string_graph, candidate_df, data_list$config)
  candidate_df <- compute_all_empirical_features(candidate_df, network_objects, data_list, data_list$config)
  print(head(candidate_df))
}
