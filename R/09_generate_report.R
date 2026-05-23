# =============================================================================
# 09_generate_report.R — Mechanism cards and HTML/Markdown report
# =============================================================================

#' Generate a single mechanism interpretation card (spec section 24)
generate_mechanism_card <- function(row, string_graph) {
  mg <- row$mutation_gene
  tg <- row$target_gene

  # Network path
  conn <- isTRUE(as.logical(row$connected))
  network_path_text <- if (!is.na(row$shortest_distance) && conn) {
    path_genes <- strsplit(row$path_genes, ";")[[1]]
    paste(path_genes, collapse = " -> ")
  } else {
    "No STRING path found"
  }

  # SL prior text
  sl_text <- switch(row$sl_prior_type,
    known_SL_pair            = "present (known SL pair)",
    neighbor_supported_SL    = "neighbor-supported",
    pathway_level_SL_support = "pathway-level support",
    "absent"
  )

  # Working hypothesis
  hypothesis <- sprintf(
    "Mutation of %s may rewire %s, increasing cellular reliance on %s. This supports %s as a candidate synthetic dependency in %s-altered contexts.",
    mg,
    if (!is.na(row$shared_pathways) && nchar(row$shared_pathways) > 0) row$shared_pathways else "cellular processes",
    tg, tg, mg
  )

  # Validation suggestions
  validation <- c(
    sprintf("Compare dependency of %s in %s-mutant vs wild-type models", tg, mg),
    sprintf("Knockout or inhibit %s in %s-mutant and wild-type cells", tg, mg),
    "Rescue mutation phenotype or measure pathway biomarkers"
  )
  if (!is.na(row$drug_category) && row$drug_category != "unknown") {
    validation <- c(validation,
      sprintf("Evaluate %s inhibitor sensitivity in %s-mutant vs wild-type cells", tg, mg))
  }

  list(
    candidate          = sprintf("%s mutation -> %s dependency", mg, tg),
    statistical        = sprintf("beta = %.3f, FDR = %.2e, delta_dependency = %.3f (n_mut = %d)",
                                 row$beta, row$FDR, row$delta_dependency, row$n_mut),
    network            = network_path_text,
    network_confidence = if (conn) sprintf("%.3f", row$path_confidence) else "N/A",
    hub_penalty        = if (conn) sprintf("%.3f", row$hub_penalty) else "N/A",
    pathway            = if (!is.na(row$shared_pathways) && nchar(row$shared_pathways) > 0) row$shared_pathways else "No shared pathway",
    signaling          = if (!is.na(row$signaling_direction) && row$signaling_direction != "none") {
      sprintf("%s (path: %s)", row$signaling_direction, row$signaling_path)
    } else "No directional signaling evidence",
    complex            = if (!is.na(row$complex_evidence)) row$complex_evidence else "No shared complex",
    sl_prior           = sl_text,
    drug               = sprintf("%s (category: %s)", row$drug_score, row$drug_category),
    final_score        = sprintf("%.3f", row$final_SL_score),
    evidence_class     = row$evidence_class,
    hypothesis         = hypothesis,
    validation         = validation
  )
}

#' Generate mechanism cards for top N candidates
generate_all_mechanism_cards <- function(ranked_df, string_graph, config) {
  top_n <- min(config$report$top_n_mechanism_cards, nrow(ranked_df))
  if (top_n == 0) {
    log_message("No candidates for mechanism cards", "WARN")
    return(list())
  }

  cards <- vector("list", top_n)
  for (i in seq_len(top_n)) {
    cards[[i]] <- generate_mechanism_card(ranked_df[i, ], string_graph)
  }
  cards
}

#' Render a single card as markdown
render_card_markdown <- function(card, i) {
  c(sprintf("## %d. %s\n", i, card$candidate),
    sprintf("**Statistical evidence:** %s\n", card$statistical),
    sprintf("**Network path:** %s\n", card$network),
    sprintf("**Path confidence:** %s | **Hub penalty:** %s\n",
            card$network_confidence, card$hub_penalty),
    sprintf("**Pathway evidence:** %s\n", card$pathway),
    sprintf("**Signaling evidence:** %s\n", card$signaling),
    sprintf("**Complex evidence:** %s\n", card$complex),
    sprintf("**SL prior:** %s\n", card$sl_prior),
    sprintf("**Druggability:** %s\n", card$drug),
    sprintf("**Final SL Score:** %s | **Class:** %s\n", card$final_score, card$evidence_class),
    sprintf("\n**Working hypothesis:** %s\n", card$hypothesis),
    "\n**Suggested validation:**\n",
    paste0("- ", card$validation, "\n"),
    "\n---\n\n"
  )
}

#' Build markdown report
render_markdown_report <- function(ranked_df, cards, config) {
  top_n <- min(config$report$top_n_mechanism_cards, nrow(ranked_df))

  lines <- c(
    "# Synthetic Dependency Prioritization Report\n",
    sprintf("**Date:** %s\n", Sys.Date()),
    sprintf("**Total candidates:** %d\n", nrow(ranked_df)),
    sprintf("**Top candidates shown:** %d\n\n", top_n),
    "## Summary Statistics\n\n"
  )

  # Class distribution
  class_counts <- table(ranked_df$evidence_class)
  lines <- c(lines, "### Evidence Class Distribution\n\n")
  for (cname in names(class_counts)) {
    lines <- c(lines, sprintf("- %s: %d\n", cname, class_counts[cname]))
  }

  # Top candidates table
  lines <- c(lines, "\n## Top Candidates\n\n")
  table_cols <- c("rank", "mutation_gene", "target_gene", "beta", "FDR",
                  "shortest_distance", "final_SL_score", "evidence_class")
  table_cols <- intersect(table_cols, names(ranked_df))
  top_rows <- min(30, nrow(ranked_df))

  # Simple pipe-separated table
  lines <- c(lines, paste(sprintf("| %s ", table_cols), collapse = ""), "|\n")
  lines <- c(lines, paste(rep("|---", length(table_cols)), collapse = ""), "|\n")
  for (i in seq_len(top_rows)) {
    vals <- sapply(table_cols, function(cn) as.character(ranked_df[i, cn]))
    lines <- c(lines, paste(sprintf("| %s ", vals), collapse = ""), "|\n")
  }

  # Mechanism cards
  lines <- c(lines, "\n## Candidate Mechanism Cards\n\n")
  for (i in seq_along(cards)) {
    lines <- c(lines, render_card_markdown(cards[[i]], i))
  }

  # Figures
  lines <- c(lines, "\n## Figures\n\n")
  lines <- c(lines,
    "![Ranking Dot Plot](../figures/ranking_dotplot.png)\n\n",
    "![Evidence Heatmap](../figures/evidence_heatmap.pdf)\n\n"
  )

  ensure_dir(config$paths$reports_dir)
  report_path <- file.path(config$paths$reports_dir, "synthetic_dependency_report.md")
  writeLines(unlist(lines), report_path)
  log_message(sprintf("Markdown report saved: %s", report_path))
}

#' Build simple HTML report
render_html_report <- function(ranked_df, cards, config) {
  top_n <- min(config$report$top_n_mechanism_cards, nrow(ranked_df))

  html <- c(
    "<!DOCTYPE html><html><head><meta charset='utf-8'>",
    "<title>Synthetic Dependency Prioritization Report</title>",
    "<style>",
    "body { font-family: Arial, sans-serif; max-width: 1000px; margin: 0 auto; padding: 20px; }",
    "h1 { color: #333; border-bottom: 2px solid #333; }",
    "h2 { color: #555; border-bottom: 1px solid #ccc; }",
    "table { border-collapse: collapse; width: 100%; font-size: 12px; }",
    "th, td { border: 1px solid #ddd; padding: 6px; text-align: left; }",
    "th { background-color: #f2f2f2; }",
    ".class-I { color: #E41A1C; font-weight: bold; }",
    ".class-II { color: #377EB8; font-weight: bold; }",
    ".class-III { color: #4DAF4A; }",
    ".class-IV { color: #FF7F00; }",
    ".card { background: #f9f9f9; border: 1px solid #ddd; padding: 15px; margin: 15px 0; border-radius: 5px; }",
    ".card h3 { margin-top: 0; }",
    ".validation { background: #e8f5e9; padding: 10px; border-left: 4px solid #4CAF50; }",
    "</style></head><body>",
    "<h1>Synthetic Dependency Prioritization Report</h1>",
    sprintf("<p><strong>Date:</strong> %s | <strong>Total candidates:</strong> %d | <strong>Top shown:</strong> %d</p>",
            Sys.Date(), nrow(ranked_df), top_n),
    "<h2>Summary Statistics</h2><ul>"
  )

  class_counts <- table(ranked_df$evidence_class)
  for (cname in names(class_counts)) {
    css_class <- if (grepl("Class I:", cname)) "class-I" else
                 if (grepl("Class II:", cname)) "class-II" else
                 if (grepl("Class III:", cname)) "class-III" else
                 if (grepl("Class IV:", cname)) "class-IV" else ""
    html <- c(html, sprintf("<li class='%s'>%s: %d</li>", css_class, cname, class_counts[cname]))
  }
  html <- c(html, "</ul>")

  # Top candidates table
  html <- c(html, "<h2>Top Candidates</h2><table><tr>")
  table_cols <- c("rank", "mutation_gene", "target_gene", "beta", "FDR",
                  "shortest_distance", "final_SL_score", "evidence_class")
  table_cols <- intersect(table_cols, names(ranked_df))
  for (cn in table_cols) {
    html <- c(html, sprintf("<th>%s</th>", cn))
  }
  html <- c(html, "</tr>")
  for (i in seq_len(min(30, nrow(ranked_df)))) {
    html <- c(html, "<tr>")
    for (cn in table_cols) {
      val <- if (is.numeric(ranked_df[i, cn])) sprintf("%.3e", ranked_df[i, cn]) else as.character(ranked_df[i, cn])
      html <- c(html, sprintf("<td>%s</td>", val))
    }
    html <- c(html, "</tr>")
  }
  html <- c(html, "</table>")

  # Mechanism cards
  html <- c(html, "<h2>Candidate Mechanism Cards</h2>")
  for (i in seq_along(cards)) {
    card <- cards[[i]]
    html <- c(html, sprintf("<div class='card'><h3>%d. %s</h3>", i, card$candidate))
    html <- c(html, sprintf("<p><strong>Statistical:</strong> %s</p>", card$statistical))
    html <- c(html, sprintf("<p><strong>Network path:</strong> %s</p>", card$network))
    html <- c(html, sprintf("<p><strong>Pathway:</strong> %s</p>", card$pathway))
    html <- c(html, sprintf("<p><strong>Signaling:</strong> %s</p>", card$signaling))
    html <- c(html, sprintf("<p><strong>SL prior:</strong> %s</p>", card$sl_prior))
    html <- c(html, sprintf("<p><strong>Druggability:</strong> %s</p>", card$drug))
    html <- c(html, sprintf("<p><strong>Final Score:</strong> %s | <strong>Class:</strong> %s</p>",
                           card$final_score, card$evidence_class))
    html <- c(html, sprintf("<p><strong>Hypothesis:</strong> %s</p>", card$hypothesis))
    html <- c(html, "<div class='validation'><strong>Suggested validation:</strong><ul>")
    for (v in card$validation) {
      html <- c(html, sprintf("<li>%s</li>", v))
    }
    html <- c(html, "</ul></div></div>")
  }

  html <- c(html, "<h2>Figures</h2>",
    "<p><img src='../figures/ranking_dotplot.png' alt='Ranking Dot Plot' style='max-width:100%'></p>",
    "</body></html>"
  )

  ensure_dir(config$paths$reports_dir)
  report_path <- file.path(config$paths$reports_dir, "synthetic_dependency_report.html")
  writeLines(unlist(html), report_path)
  log_message(sprintf("HTML report saved: %s", report_path))
}

#' Generate full report
generate_report <- function(ranked_df, string_graph, config) {
  log_message("=== Step 9: Generate Report ===")

  if (nrow(ranked_df) == 0) {
    log_message("No candidates for report generation", "WARN")
    return(invisible(NULL))
  }

  cards <- generate_all_mechanism_cards(ranked_df, string_graph, config)

  fmt <- config$report$output_format
  if (fmt %in% c("markdown", "both")) {
    render_markdown_report(ranked_df, cards, config)
  }
  if (fmt %in% c("html", "both")) {
    render_html_report(ranked_df, cards, config)
  }

  log_message("=== Report generation complete ===")
}

# -- Standalone execution ---------------------------------------------------
if (sys.nframe() == 0) {
  source("R/00_utils.R")
  config <- yaml::read_yaml("config/config.yaml")
  ranked_df <- read.csv(config$output_files$final_ranked, stringsAsFactors = FALSE)
  generate_report(ranked_df, NULL, config)
}
