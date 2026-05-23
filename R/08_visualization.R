# =============================================================================
# 08_visualization.R — Publication-quality figures
# =============================================================================

# Color palette for evidence classes
evidence_class_colors <- c(
  "Class I: high-confidence synthetic dependency candidate" = "#E41A1C",
  "Class II: network-supported candidate"                   = "#377EB8",
  "Class III: SL-prior-supported candidate"                 = "#4DAF4A",
  "Class IV: statistical-only candidate"                    = "#FF7F00",
  "Low priority"                                            = "#999999"
)

#' Ranking dot plot
plot_ranking_dotplot <- function(ranked_df, config) {
  top_n <- min(config$visualization$top_n_candidates, nrow(ranked_df))
  if (top_n == 0) {
    log_message("No candidates for dot plot", "WARN")
    return(invisible(NULL))
  }

  plot_df <- ranked_df[1:top_n, ]
  plot_df$label <- paste0(plot_df$mutation_gene, " -> ", plot_df$target_gene)
  plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))
  plot_df$nlogFDR <- safe_log10(plot_df$FDR)

  pal <- evidence_class_colors
  available_classes <- intersect(names(pal), unique(plot_df$evidence_class))
  pal <- pal[available_classes]

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x = final_SL_score, y = label,
    color = evidence_class, size = nlogFDR
  )) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_manual(values = pal, name = "Evidence Class") +
    ggplot2::scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +
    ggplot2::labs(
      title = "Synthetic Dependency Candidate Ranking",
      x = "Final SL Score",
      y = ""
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      legend.box = "vertical"
    )

  ensure_dir(config$paths$figures_dir)
  w <- config$visualization$dot_plot$width
  h <- config$visualization$dot_plot$height
  ggplot2::ggsave(file.path(config$paths$figures_dir, "ranking_dotplot.pdf"),
                  p, width = w, height = h, device = "pdf")
  ggplot2::ggsave(file.path(config$paths$figures_dir, "ranking_dotplot.png"),
                  p, width = w, height = h, dpi = 150)

  log_message("Ranking dot plot saved")
  invisible(p)
}

#' Network path plot for a single candidate pair
plot_network_path <- function(ranked_df, string_graph, pair_index, config) {
  if (is.null(string_graph)) {
    log_message("No STRING graph available for path plot", "WARN")
    return(invisible(NULL))
  }

  row <- ranked_df[pair_index, ]
  mg <- row$mutation_gene
  tg <- row$target_gene
  vnames <- igraph::V(string_graph)$name

  if (!(mg %in% vnames) || !(tg %in% vnames)) {
    log_message(sprintf("Genes not in graph: %s, %s", mg, tg), "WARN")
    return(invisible(NULL))
  }

  sp <- safe_shortest_path(string_graph, mg, tg,
                           weights = igraph::E(string_graph)$distance_weight)
  if (is.null(sp)) {
    log_message(sprintf("No path between %s and %s", mg, tg), "WARN")
    return(invisible(NULL))
  }

  # Extract subgraph
  path_nodes <- sp$path_genes
  subg <- igraph::induced_subgraph(string_graph, path_nodes)

  # Node colors
  node_colors <- rep("gray60", length(path_nodes))
  node_colors[path_nodes == mg] <- "#E41A1C"   # mutation gene - red
  node_colors[path_nodes == tg] <- "#377EB8"   # target gene - blue

  # Node sizes
  node_sizes <- rep(5, length(path_nodes))
  node_sizes[path_nodes %in% c(mg, tg)] <- 8

  # Layout
  if (length(path_nodes) <= 5) {
    lay <- igraph::layout_with_fr(subg)
  } else {
    lay <- igraph::layout_as_tree(subg)
  }

  ensure_dir(config$paths$figures_dir)
  fname <- sprintf("network_path_%s_%s.pdf", mg, tg)
  w <- config$visualization$network_path_plot$width
  h <- config$visualization$network_path_plot$height

  pdf(file.path(config$paths$figures_dir, fname), width = w, height = h)
  igraph::plot.igraph(subg,
    layout      = lay,
    vertex.color = node_colors,
    vertex.size  = node_sizes,
    vertex.label = path_nodes,
    vertex.label.cex = 0.7,
    edge.width   = igraph::E(subg)$edge_confidence * 3,
    edge.color   = "gray50",
    main         = sprintf("%s -> %s (distance=%d, confidence=%.3f)",
                           mg, tg, sp$distance,
                           compute_path_confidence(sp$edge_confs, "product", 5))
  )
  dev.off()

  invisible(subg)
}

#' Plot all network paths for top examples
plot_all_network_paths <- function(ranked_df, string_graph, config) {
  max_examples <- config$visualization$network_path_plot$max_examples
  n_to_plot <- min(max_examples, nrow(ranked_df))

  if (n_to_plot == 0) {
    log_message("No candidates for network path plots", "WARN")
    return(invisible(NULL))
  }

  n_plotted <- 0
  for (i in seq_len(n_to_plot)) {
    conn <- as.logical(ranked_df$connected[i])
    if (isTRUE(conn)) {
      tryCatch({
        plot_network_path(ranked_df, string_graph, i, config)
        n_plotted <- n_plotted + 1
      }, error = function(e) {
        log_message(sprintf("Failed to plot path for row %d: %s", i, e$message), "WARN")
      })
    }
  }

  log_message(sprintf("Network path plots: %d generated", n_plotted))
  invisible(NULL)
}

#' Evidence heatmap
plot_evidence_heatmap <- function(ranked_df, config) {
  top_n <- min(config$visualization$top_n_candidates, nrow(ranked_df))
  if (top_n == 0) {
    log_message("No candidates for heatmap", "WARN")
    return(invisible(NULL))
  }

  plot_df <- ranked_df[1:top_n, ]

  score_cols <- c("stat_score", "network_score_scaled", "pathway_score",
                  "signaling_score", "complex_score", "sl_prior_score", "drug_score")
  score_cols <- intersect(score_cols, names(plot_df))

  if (length(score_cols) == 0) {
    log_message("No score columns for heatmap", "WARN")
    return(invisible(NULL))
  }

  # Build matrix
  heat_mat <- as.matrix(plot_df[, score_cols, drop = FALSE])
  rownames(heat_mat) <- paste0(plot_df$mutation_gene, " -> ", plot_df$target_gene)

  # Row annotations
  ann_colors <- evidence_class_colors
  available_classes <- intersect(names(ann_colors), unique(plot_df$evidence_class))

  ann_row <- data.frame(
    Class = plot_df$evidence_class,
    row.names = rownames(heat_mat)
  )
  ann_col <- list(Class = ann_colors[available_classes])

  ensure_dir(config$paths$figures_dir)
  w <- config$visualization$heatmap$width
  h <- config$visualization$heatmap$height

  pdf(file.path(config$paths$figures_dir, "evidence_heatmap.pdf"), width = w, height = h)

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(heat_mat,
      color         = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
      annotation_row = ann_row,
      annotation_colors = ann_col,
      cluster_rows   = TRUE,
      cluster_cols   = FALSE,
      fontsize_row   = 6,
      fontsize_col   = 10,
      main           = "Evidence Score Heatmap"
    )
  } else {
    # Fallback using ggplot2
    heat_long <- as.data.frame(heat_mat)
    heat_long$pair <- rownames(heat_long)
    heat_long <- tidyr::pivot_longer(heat_long, -pair, names_to = "evidence", values_to = "score")
    heat_long$pair <- factor(heat_long$pair, levels = rev(unique(heat_long$pair)))

    p <- ggplot2::ggplot(heat_long, ggplot2::aes(x = evidence, y = pair, fill = score)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "white", high = "darkblue", na.value = "grey90") +
      ggplot2::labs(title = "Evidence Score Heatmap", x = "", y = "") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    print(p)
  }

  dev.off()
  log_message("Evidence heatmap saved")
  invisible(NULL)
}

#' Score distribution plot
plot_score_distributions <- function(ranked_df) {
  score_cols <- c("stat_score", "network_score_scaled", "pathway_score",
                  "signaling_score", "complex_score", "sl_prior_score", "drug_score",
                  "final_SL_score")
  score_cols <- intersect(score_cols, names(ranked_df))

  if (length(score_cols) == 0) return(invisible(NULL))

  plot_long <- ranked_df[, score_cols, drop = FALSE]
  plot_long <- tidyr::pivot_longer(plot_long, dplyr::everything(),
                                   names_to = "score_type", values_to = "value")

  p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = value, fill = score_type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ score_type, scales = "free", ncol = 3) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = "Score Distributions", x = "Score", y = "Density")

  ensure_dir("results/figures")
  ggplot2::ggsave("results/figures/score_distributions.pdf", p, width = 12, height = 8)
  log_message("Score distribution plot saved")
  invisible(p)
}

#' Generate all figures
generate_all_figures <- function(ranked_df, string_graph, config) {
  log_message("=== Step 8: Visualization ===")

  if (nrow(ranked_df) == 0) {
    log_message("No candidates for visualization", "WARN")
    return(invisible(NULL))
  }

  plot_ranking_dotplot(ranked_df, config)
  plot_evidence_heatmap(ranked_df, config)
  plot_all_network_paths(ranked_df, string_graph, config)
  plot_score_distributions(ranked_df)

  log_message("=== Visualization complete ===")
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  config <- yaml::read_yaml("config/config.yaml")
  ranked_df <- read.csv(config$output_files$final_ranked, stringsAsFactors = FALSE)
  generate_all_figures(ranked_df, NULL, config)
}
