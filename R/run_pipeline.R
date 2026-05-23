#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R — Main CLI orchestrator
#
# Usage:
#   Rscript R/run_pipeline.R --simulate --verbose
#   Rscript R/run_pipeline.R --simulate --stop-at 1
#   Rscript R/run_pipeline.R --start-at 2 --skip-permutation
#   Rscript R/run_pipeline.R --config config/config.yaml
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# -- CLI Option Definition --------------------------------------------------

option_list <- list(
  make_option("--config", default = "config/config.yaml",
              help = "Path to YAML configuration file [default: %default]"),
  make_option("--simulate", action = "store_true", default = FALSE,
              help = "Generate simulated data before running"),
  make_option("--skip-association", action = "store_true", default = FALSE,
              help = "Skip association analysis (use existing results)"),
  make_option("--skip-permutation", action = "store_true", default = FALSE,
              help = "Skip permutation testing (slow)"),
  make_option("--skip-visualization", action = "store_true", default = FALSE,
              help = "Skip visualization generation"),
  make_option("--skip-report", action = "store_true", default = FALSE,
              help = "Skip report generation"),
  make_option("--start-at", type = "integer", default = 1,
              help = "Start pipeline at step N (1-9) [default: %default]"),
  make_option("--stop-at", type = "integer", default = 9,
              help = "Stop pipeline after step N (1-9) [default: %default]"),
  make_option("--verbose", action = "store_true", default = FALSE,
              help = "Verbose logging")
)

# -- Main -------------------------------------------------------------------

main <- function() {
  # Parse CLI
  opt <- parse_args(OptionParser(option_list = option_list))

  # Resolve config path relative to project root
  if (!grepl("^/", opt$config)) {
    opt$config <- file.path(getwd(), opt$config)
  }

  # Load utilities
  source("R/00_utils.R")
  config <- yaml::read_yaml(opt$config)

  log_message("============================================")
  log_message("Network-Constrained Synthetic Dependency Prioritization Pipeline")
  log_message("============================================")

  # Set up logging
  ensure_dir(config$paths$log_dir)

  # -- Step 0: Simulate data (optional) ------------------------------------
  if (opt$simulate) {
    log_message("=== Step 0: Simulate Data ===")
    source("R/simulate_data.R")
    simulate_main(config)
  }

  # -- Pipeline state ------------------------------------------------------
  data_list       <- NULL
  candidate_df    <- NULL
  network_objects <- NULL
  ranked_df       <- NULL

  for (step in opt$`start-at`:opt$`stop-at`) {

    # -- Step 1: Prepare Data ----------------------------------------------
    if (step == 1) {
      log_message("=== Step 1: Prepare Data ===")
      source("R/01_prepare_data.R")
      data_list <- prepare_data(opt$config)
    }

    # -- Step 2: Association Analysis --------------------------------------
    else if (step == 2 && !opt$`skip-association`) {
      if (is.null(data_list)) {
        log_message("Loading data for association step...")
        source("R/01_prepare_data.R")
        data_list <- prepare_data(opt$config)
      }
      log_message("=== Step 2: Association Analysis ===")
      source("R/02_association_analysis.R")
      candidate_df <- run_association_pipeline(data_list)
    }

    # -- Step 3: Build Networks --------------------------------------------
    else if (step == 3) {
      if (is.null(data_list)) {
        source("R/01_prepare_data.R")
        data_list <- prepare_data(opt$config)
      }
      log_message("=== Step 3: Build Networks ===")
      source("R/03_build_networks.R")
      network_objects <- build_all_networks(data_list)
    }

    # -- Step 4: Network Features ------------------------------------------
    else if (step == 4) {
      if (is.null(candidate_df)) {
        cf <- config$output_files$candidate_pairs
        if (file.exists(cf)) {
          candidate_df <- read.csv(cf, stringsAsFactors = FALSE)
          log_message(sprintf("Loaded %d candidates from %s", nrow(candidate_df), cf))
        } else {
          stop("No candidate pairs found. Run steps 1-2 first.")
        }
      }
      if (is.null(network_objects)) {
        if (is.null(data_list)) {
          source("R/01_prepare_data.R")
          data_list <- prepare_data(opt$config)
        }
        source("R/03_build_networks.R")
        network_objects <- build_all_networks(data_list)
      }
      log_message("=== Step 4: Compute Network Features ===")
      source("R/04_compute_network_features.R")
      candidate_df <- compute_network_features_batch(
        network_objects$string_graph, candidate_df, config
      )
    }

    # -- Step 5: Empirical Evidence ----------------------------------------
    else if (step == 5) {
      if (is.null(candidate_df)) {
        nf <- config$output_files$network_features
        if (file.exists(nf)) {
          candidate_df <- read.csv(nf, stringsAsFactors = FALSE)
          log_message(sprintf("Loaded enriched candidates from %s", nf))
        } else {
          stop("No network features found. Run steps 1-4 first.")
        }
      }
      if (is.null(network_objects)) {
        if (is.null(data_list)) {
          source("R/01_prepare_data.R")
          data_list <- prepare_data(opt$config)
        }
        source("R/03_build_networks.R")
        network_objects <- build_all_networks(data_list)
      }
      log_message("=== Step 5: Compute Empirical Evidence ===")
      source("R/05_compute_empirical_evidence.R")
      candidate_df <- compute_all_empirical_features(
        candidate_df, network_objects, data_list, config
      )
    }

    # -- Step 6: Score Candidates ------------------------------------------
    else if (step == 6) {
      if (is.null(candidate_df)) {
        ef <- config$output_files$empirical_features
        if (file.exists(ef)) {
          candidate_df <- read.csv(ef, stringsAsFactors = FALSE)
          log_message(sprintf("Loaded features from %s", ef))
        } else {
          stop("No empirical features found. Run steps 1-5 first.")
        }
      }
      log_message("=== Step 6: Score Candidates ===")
      source("R/06_score_candidates.R")
      ranked_df <- score_and_rank(candidate_df, config)
    }

    # -- Step 7: Permutation Tests -----------------------------------------
    else if (step == 7 && !opt$`skip-permutation`) {
      if (is.null(ranked_df)) {
        fr <- config$output_files$final_ranked
        if (file.exists(fr)) {
          ranked_df <- read.csv(fr, stringsAsFactors = FALSE)
          log_message(sprintf("Loaded ranked candidates from %s", fr))
        } else {
          stop("No ranked candidates found. Run steps 1-6 first.")
        }
      }
      if (is.null(network_objects)) {
        if (is.null(data_list)) {
          source("R/01_prepare_data.R")
          data_list <- prepare_data(opt$config)
        }
        source("R/03_build_networks.R")
        network_objects <- build_all_networks(data_list)
      }
      if (is.null(data_list)) {
        source("R/01_prepare_data.R")
        data_list <- prepare_data(opt$config)
      }
      log_message("=== Step 7: Permutation Tests ===")
      source("R/07_null_model_permutation.R")
      ranked_df <- run_permutation_tests(
        ranked_df, network_objects$string_graph, data_list$dependency_mat, config
      )
    }

    # -- Step 8: Visualization ---------------------------------------------
    else if (step == 8 && !opt$`skip-visualization`) {
      if (is.null(ranked_df)) {
        fr <- config$output_files$final_ranked
        if (file.exists(fr)) {
          ranked_df <- read.csv(fr, stringsAsFactors = FALSE)
          log_message(sprintf("Loaded ranked candidates from %s", fr))
        } else {
          stop("No ranked candidates found. Run steps 1-6 first.")
        }
      }
      if (is.null(network_objects)) {
        if (is.null(data_list)) {
          source("R/01_prepare_data.R")
          data_list <- prepare_data(opt$config)
        }
        source("R/03_build_networks.R")
        network_objects <- build_all_networks(data_list)
      }
      log_message("=== Step 8: Visualization ===")
      source("R/08_visualization.R")
      generate_all_figures(ranked_df, network_objects$string_graph, config)
    }

    # -- Step 9: Report ----------------------------------------------------
    else if (step == 9 && !opt$`skip-report`) {
      if (is.null(ranked_df)) {
        fr <- config$output_files$final_ranked
        if (file.exists(fr)) {
          ranked_df <- read.csv(fr, stringsAsFactors = FALSE)
          log_message(sprintf("Loaded ranked candidates from %s", fr))
        } else {
          stop("No ranked candidates found. Run steps 1-6 first.")
        }
      }
      if (is.null(network_objects)) {
        if (is.null(data_list)) {
          source("R/01_prepare_data.R")
          data_list <- prepare_data(opt$config)
        }
        source("R/03_build_networks.R")
        network_objects <- build_all_networks(data_list)
      }
      log_message("=== Step 9: Generate Report ===")
      source("R/09_generate_report.R")
      generate_report(ranked_df, network_objects$string_graph, config)
    }
  }

  log_message("============================================")
  log_message("Pipeline Complete")
  log_message("============================================")
}

# -- Run --------------------------------------------------------------------
main()
