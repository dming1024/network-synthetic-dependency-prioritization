#!/usr/bin/env Rscript
# =============================================================================
# simulate_data.R — Generate realistic synthetic test data
#
# Creates all 10 input data files with injected true-positive SL signals
# so the pipeline can be developed and validated without real genomic data.
# =============================================================================

# -- Helpers ----------------------------------------------------------------

#' Generate random gene names
random_genes <- function(n, prefix) {
  sprintf("%s_%d", prefix, seq_len(n))
}

#' Draw from Beta-like distribution for STRING scores
rstring_score <- function(n, shape1 = 5, shape2 = 2) {
  round(rbeta(n, shape1, shape2) * 1000)
}

# -- Generator functions ----------------------------------------------------

generate_sample_metadata <- function(n_samples, seed) {
  set.seed(seed)
  lineages <- c("LUAD", "BRCA", "COAD", "GBM", "SKCM", "OV", "PRAD", "LAML", "HNSC", "KIRC")
  data.frame(
    sample_id   = sprintf("S%d", seq_len(n_samples)),
    lineage     = sample(lineages, n_samples, replace = TRUE),
    TMB         = round(10^rnorm(n_samples, mean = 0.5, sd = 0.4), 1),
    MSI_status  = rbinom(n_samples, 1, 0.05),
    CNV_burden  = round(rnorm(n_samples, mean = 0.15, sd = 0.08), 3),
    batch       = sample(1:3, n_samples, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

generate_mutation_matrix <- function(samples, n_genes, mutation_rate_mean, seed) {
  set.seed(seed)
  n <- length(samples)
  gene_names <- random_genes(n_genes, "MUT")
  # Some genes have higher mutation rates (simulate cancer genes)
  rates <- rbeta(n_genes, 2, 30) * mutation_rate_mean * 20
  rates <- pmax(rates, 0.005)

  mat <- matrix(0, nrow = n, ncol = n_genes)
  for (j in seq_len(n_genes)) {
    mat[, j] <- rbinom(n, 1, rates[j])
  }
  # Ensure no column is all-zero
  for (j in seq_len(n_genes)) {
    if (sum(mat[, j]) == 0) {
      mat[sample(n, 1), j] <- 1
    }
  }
  colnames(mat) <- gene_names
  cbind(data.frame(sample_id = samples, stringsAsFactors = FALSE), as.data.frame(mat))
}

generate_dependency_matrix <- function(samples, n_genes, dep_mean, dep_sd, seed) {
  set.seed(seed)
  n <- length(samples)
  gene_names <- random_genes(n_genes, "DEP")

  mat <- matrix(0, nrow = n, ncol = n_genes)
  for (j in seq_len(n_genes)) {
    gmean <- rnorm(1, mean = dep_mean, sd = 0.3)
    gsd   <- dep_sd * runif(1, 0.5, 2)
    # Some genes are pan-essential (more negative mean)
    if (runif(1) < 0.1) {
      gmean <- gmean - 0.6
    }
    mat[, j] <- rnorm(n, mean = gmean, sd = gsd)
  }
  colnames(mat) <- gene_names
  cbind(data.frame(sample_id = samples, stringsAsFactors = FALSE), as.data.frame(mat))
}

inject_true_sl_signals <- function(mut_mat, dep_mat, meta, n_true_pairs, effect_size, seed) {
  set.seed(seed)

  mut_genes  <- setdiff(names(mut_mat), "sample_id")
  dep_genes  <- setdiff(names(dep_mat), "sample_id")
  n_true     <- min(n_true_pairs, length(mut_genes) * length(dep_genes))
  true_pairs <- data.frame(
    mutation_gene = sample(mut_genes, n_true, replace = TRUE),
    target_gene   = sample(dep_genes, n_true, replace = TRUE),
    stringsAsFactors = FALSE
  )
  true_pairs <- true_pairs[true_pairs$mutation_gene != true_pairs$target_gene, ]
  true_pairs <- true_pairs[!duplicated(true_pairs), ]

  for (i in seq_len(nrow(true_pairs))) {
    mg <- true_pairs$mutation_gene[i]
    tg <- true_pairs$target_gene[i]
    mut_samples <- which(mut_mat[[mg]] == 1)
    if (length(mut_samples) > 0) {
      dep_mat[mut_samples, tg] <- dep_mat[mut_samples, tg] + effect_size
    }
  }

  list(dependency_mat = dep_mat, true_pairs = true_pairs)
}

generate_expression_matrix <- function(samples, mut_genes, dep_genes, seed) {
  set.seed(seed)
  n <- length(samples)
  all_genes <- unique(c(mut_genes, dep_genes))
  mat <- matrix(0, nrow = n, ncol = length(all_genes))
  for (j in seq_along(all_genes)) {
    mat[, j] <- 2^rnorm(n, mean = 5, sd = 1.5)
  }
  colnames(mat) <- all_genes
  cbind(data.frame(sample_id = samples, stringsAsFactors = FALSE), as.data.frame(mat))
}

generate_string_edges <- function(mut_genes, dep_genes, true_pairs, coverage, seed) {
  set.seed(seed)
  all_genes <- unique(c(mut_genes, dep_genes))
  n_total   <- length(all_genes)
  n_covered <- round(n_total * coverage)
  covered   <- sample(all_genes, n_covered)

  edges <- data.frame(
    geneA = character(0),
    geneB = character(0),
    combined_score = integer(0),
    stringsAsFactors = FALSE
  )

  # For true SL pairs, ensure a 1-hop or 2-hop path exists
  if (nrow(true_pairs) > 0) {
    for (i in seq_len(nrow(true_pairs))) {
      mg <- true_pairs$mutation_gene[i]
      tg <- true_pairs$target_gene[i]
      if (!(mg %in% covered) || !(tg %in% covered)) next
      if (runif(1) < 0.3 && mg != tg) {
        # Direct edge
        edges <- rbind(edges, data.frame(
          geneA = mg, geneB = tg,
          combined_score = rstring_score(1, 8, 2),
          stringsAsFactors = FALSE
        ))
      } else {
        # 2-hop: mg -> bridge -> tg
        bridge <- paste0("BRIDGE_", mg, "_", tg)
        if (bridge %in% covered) next  # avoid collision
        covered <- c(covered, bridge)
        edges <- rbind(edges, data.frame(
          geneA = c(mg, bridge),
          geneB = c(bridge, tg),
          combined_score = rstring_score(2, 7, 2),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Background edges
  n_bg <- round(n_covered * 3)
  for (k in seq_len(n_bg)) {
    pair <- sample(covered, 2, replace = FALSE)
    if (pair[1] == pair[2]) next
    edges <- rbind(edges, data.frame(
      geneA = pair[1],
      geneB = pair[2],
      combined_score = rstring_score(1),
      stringsAsFactors = FALSE
    ))
  }

  edges <- edges[!duplicated(edges), ]
  edges <- edges[edges$geneA != edges$geneB, ]
  edges
}

generate_reactome_pathways <- function(mut_genes, dep_genes, n_pathways, true_pairs, seed) {
  set.seed(seed)
  all_genes <- unique(c(mut_genes, dep_genes))
  pathways  <- sprintf("REACTOME_PATHWAY_%d", seq_len(n_pathways))

  mapping <- data.frame(
    gene    = character(0),
    pathway = character(0),
    stringsAsFactors = FALSE
  )

  # Assign random genes to pathways
  for (pw in pathways) {
    n_genes <- sample(10:50, 1)
    genes   <- sample(all_genes, min(n_genes, length(all_genes)))
    mapping <- rbind(mapping, data.frame(
      gene = genes, pathway = pw, stringsAsFactors = FALSE
    ))
  }

  # For true SL pairs, ensure shared pathway
  if (nrow(true_pairs) > 0) {
    for (i in seq_len(nrow(true_pairs))) {
      mg <- true_pairs$mutation_gene[i]
      tg <- true_pairs$target_gene[i]
      if (runif(1) < 0.6) {
        pw <- sample(pathways, 1)
        mapping <- rbind(mapping, data.frame(
          gene = c(mg, tg), pathway = pw, stringsAsFactors = FALSE
        ))
      }
    }
  }

  mapping[!duplicated(mapping), ]
}

generate_omnipath_edges <- function(mut_genes, dep_genes, true_pairs, seed) {
  set.seed(seed)
  all_genes <- unique(c(mut_genes, dep_genes))

  edges <- data.frame(
    source_gene = character(0),
    target_gene = character(0),
    effect      = character(0),
    confidence  = numeric(0),
    stringsAsFactors = FALSE
  )

  # For some true pairs, create a directed path
  effects <- c("activation", "inhibition", "unknown")
  if (nrow(true_pairs) > 0) {
    for (i in seq_len(nrow(true_pairs))) {
      mg <- true_pairs$mutation_gene[i]
      tg <- true_pairs$target_gene[i]
      if (runif(1) < 0.25 && mg != tg) {
        # Direct directed edge
        edges <- rbind(edges, data.frame(
          source_gene = mg, target_gene = tg,
          effect = sample(effects, 1),
          confidence = round(runif(1, 0.5, 1), 2),
          stringsAsFactors = FALSE
        ))
      } else if (runif(1) < 0.25) {
        # 2-hop directed
        bridge <- paste0("SIGBRIDGE_", mg, "_", tg)
        edges <- rbind(edges, data.frame(
          source_gene = c(mg, bridge),
          target_gene = c(bridge, tg),
          effect = sample(effects, 2, replace = TRUE),
          confidence = round(runif(2, 0.5, 1), 2),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Background directed edges
  n_bg <- min(500, length(all_genes) * 2)
  for (k in seq_len(n_bg)) {
    pair <- sample(all_genes, 2, replace = FALSE)
    if (pair[1] == pair[2]) next
    edges <- rbind(edges, data.frame(
      source_gene = pair[1], target_gene = pair[2],
      effect = sample(effects, 1),
      confidence = round(runif(1, 0.3, 1), 2),
      stringsAsFactors = FALSE
    ))
  }

  edges[!duplicated(edges), ]
}

generate_corum_complexes <- function(mut_genes, dep_genes, n_complexes, true_pairs, seed) {
  set.seed(seed)
  all_genes <- unique(c(mut_genes, dep_genes))

  mapping <- data.frame(
    gene       = character(0),
    complex_id = character(0),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_complexes)) {
    cid   <- sprintf("COMPLEX_%d", i)
    n_g   <- sample(3:15, 1)
    genes <- sample(all_genes, min(n_g, length(all_genes)))
    mapping <- rbind(mapping, data.frame(
      gene = genes, complex_id = cid, stringsAsFactors = FALSE
    ))
  }

  # Some true pairs share a complex
  if (nrow(true_pairs) > 0) {
    for (i in seq_len(nrow(true_pairs))) {
      mg <- true_pairs$mutation_gene[i]
      tg <- true_pairs$target_gene[i]
      if (runif(1) < 0.15) {
        cid <- sprintf("COMPLEX_SL_%d", i)
        mapping <- rbind(mapping, data.frame(
          gene = c(mg, tg), complex_id = cid, stringsAsFactors = FALSE
        ))
      }
    }
  }

  mapping[!duplicated(mapping), ]
}

generate_synlethdb_pairs <- function(mut_genes, dep_genes, true_pairs, seed) {
  set.seed(seed)
  all_genes <- unique(c(mut_genes, dep_genes))

  pairs <- data.frame(
    geneA        = character(0),
    geneB        = character(0),
    evidence_type = character(0),
    stringsAsFactors = FALSE
  )

  # Include some true pairs directly in SL DB
  if (nrow(true_pairs) > 0) {
    n_known <- min(round(nrow(true_pairs) * 0.4), nrow(true_pairs))
    known   <- true_pairs[sample(nrow(true_pairs), n_known), ]
    pairs <- rbind(pairs, data.frame(
      geneA = known$mutation_gene,
      geneB = known$target_gene,
      evidence_type = "literature",
      stringsAsFactors = FALSE
    ))
  }

  # Background SL pairs
  n_bg <- 200
  for (k in seq_len(n_bg)) {
    pair <- sample(all_genes, 2, replace = FALSE)
    if (pair[1] == pair[2]) next
    pairs <- rbind(pairs, data.frame(
      geneA = pair[1], geneB = pair[2],
      evidence_type = "database",
      stringsAsFactors = FALSE
    ))
  }

  pairs[!duplicated(pairs), ]
}

generate_druggable_targets <- function(dep_genes, seed) {
  set.seed(seed)
  n <- length(dep_genes)
  categories <- c("FDA_approved", "clinical", "preclinical", "druggable_class", "unknown")
  probs      <- c(0.05, 0.08, 0.12, 0.25, 0.50)

  data.frame(
    gene          = dep_genes,
    drug_category = sample(categories, n, replace = TRUE, prob = probs),
    drug_name     = ifelse(runif(n) < 0.5,
                           sprintf("Drug_%s", dep_genes),
                           NA_character_),
    stringsAsFactors = FALSE
  )
}

# -- Main -------------------------------------------------------------------

simulate_main <- function(config) {
  log_message("=== Generating simulated data ===")

  sim <- config$simulation
  raw_dir <- config$paths$raw_dir
  ensure_dir(raw_dir)

  # 1. Sample metadata
  log_message("Generating sample metadata...")
  meta <- generate_sample_metadata(sim$n_samples, sim$random_seed)
  write.csv(meta, file.path(raw_dir, "sample_metadata.csv"), row.names = FALSE)

  # 2. Mutation matrix
  log_message("Generating mutation matrix...")
  mut <- generate_mutation_matrix(meta$sample_id, sim$n_mutation_genes,
                                  sim$mutation_rate_mean, sim$random_seed + 1)
  write.csv(mut, file.path(raw_dir, "mutation_matrix.csv"), row.names = FALSE)

  # 3. Dependency matrix
  log_message("Generating dependency matrix...")
  dep <- generate_dependency_matrix(meta$sample_id, sim$n_dependency_genes,
                                    sim$dependency_mean, sim$dependency_sd,
                                    sim$random_seed + 2)

  # 4. Extract gene lists
  mut_genes <- setdiff(names(mut), "sample_id")
  dep_genes <- setdiff(names(dep), "sample_id")

  # 5. Inject true SL signals
  log_message("Injecting true SL signals...")
  res <- inject_true_sl_signals(mut, dep, meta, sim$n_true_sl_pairs,
                                sim$true_sl_effect_size, sim$random_seed + 3)
  dep <- res$dependency_mat
  true_pairs <- res$true_pairs
  log_message(sprintf("  Injected %d true SL pairs", nrow(true_pairs)))
  write.csv(dep, file.path(raw_dir, "dependency_matrix.csv"), row.names = FALSE)

  # 6. Expression matrix
  log_message("Generating expression matrix...")
  expr <- generate_expression_matrix(meta$sample_id, mut_genes, dep_genes,
                                     sim$random_seed + 4)
  write.csv(expr, file.path(raw_dir, "expression_matrix.csv"), row.names = FALSE)

  # 7. STRING edges
  log_message("Generating STRING edges...")
  string_edges <- generate_string_edges(mut_genes, dep_genes, true_pairs,
                                        sim$string_coverage, sim$random_seed + 5)
  write.csv(string_edges, file.path(raw_dir, "STRING_edges.csv"), row.names = FALSE)
  log_message(sprintf("  %d STRING edges", nrow(string_edges)))

  # 8. Reactome pathways
  log_message("Generating Reactome pathway mapping...")
  reactome <- generate_reactome_pathways(mut_genes, dep_genes, sim$n_pathway_terms,
                                         true_pairs, sim$random_seed + 6)
  write.csv(reactome, file.path(raw_dir, "Reactome_gene_pathway.csv"), row.names = FALSE)

  # 9. OmniPath edges
  log_message("Generating OmniPath signaling edges...")
  omnipath <- generate_omnipath_edges(mut_genes, dep_genes, true_pairs,
                                      sim$random_seed + 7)
  write.csv(omnipath, file.path(raw_dir, "OmniPath_edges.csv"), row.names = FALSE)

  # 10. CORUM complexes
  log_message("Generating CORUM complex data...")
  corum <- generate_corum_complexes(mut_genes, dep_genes, sim$n_complexes,
                                    true_pairs, sim$random_seed + 8)
  write.csv(corum, file.path(raw_dir, "CORUM_complexes.csv"), row.names = FALSE)

  # 11. SynLethDB pairs
  log_message("Generating SynLethDB pairs...")
  sl_db <- generate_synlethdb_pairs(mut_genes, dep_genes, true_pairs,
                                    sim$random_seed + 9)
  write.csv(sl_db, file.path(raw_dir, "SynLethDB_pairs.csv"), row.names = FALSE)

  # 12. Druggable targets
  log_message("Generating druggable target list...")
  drug <- generate_druggable_targets(dep_genes, sim$random_seed + 10)
  write.csv(drug, file.path(raw_dir, "druggable_targets.csv"), row.names = FALSE)

  # Save true pairs for validation
  write.csv(true_pairs, file.path(raw_dir, "true_sl_pairs.csv"), row.names = FALSE)

  log_message("=== Simulation complete ===")
  log_message(sprintf("Files written to %s/", raw_dir))
  invisible(true_pairs)
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  config <- yaml::read_yaml("config/config.yaml")
  simulate_main(config)
}
