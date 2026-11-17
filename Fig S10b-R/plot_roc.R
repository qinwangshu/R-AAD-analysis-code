# Draw the three-panel covariate-adjusted ROC figure from aggregate data.
# Run with Rscript plot_roc.R, or open this file in RStudio and click Source.

script_directory <- function() {
  source_files <- vapply(sys.frames(), function(frame) {
    if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else ""
  }, character(1))
  source_files <- source_files[nzchar(source_files)]
  if (length(source_files)) return(dirname(normalizePath(tail(source_files, 1L))))

  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg)) return(dirname(normalizePath(sub("^--file=", "", script_arg[1L]))))
  getwd()
}

xml_text <- function(value) {
  value <- gsub("&", "&amp;", value, fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  gsub(">", "&gt;", value, fixed = TRUE)
}

number <- function(value) formatC(value, digits = 5, format = "f", decimal.mark = ".")

svg_text <- function(x, y, label, size = 11, anchor = "start", weight = "normal",
                     extra = "") {
  sprintf(
    '<text x="%s" y="%s" font-family="Arial" font-size="%s" text-anchor="%s" fill="#111111" font-weight="%s" %s>%s</text>',
    number(x), number(y), size, anchor, weight, extra, xml_text(label)
  )
}

svg_line <- function(x1, y1, x2, y2) {
  sprintf(
    '<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="#111111" stroke-width="1.1"/>',
    number(x1), number(y1), number(x2), number(y2)
  )
}

svg_curve <- function(x, y) {
  points <- paste(paste(number(x), number(y), sep = ","), collapse = " ")
  sprintf(
    '<polyline points="%s" fill="none" stroke="#F02222" stroke-width="1.55" stroke-linejoin="miter"/>',
    points
  )
}

read_inputs <- function(directory) {
  summary <- read.csv(file.path(directory, "ROC_summary.csv"))
  coordinates <- read.csv(file.path(directory, "ROC_coordinates.csv"))
  groups <- c("CTRL", "HP", "AMI")
  comparisons <- paste("AAD vs", groups)
  required_summary <- c(
    "panel", "comparison", "n_AAD", "n_comparator", "AUC", "CI95_lower",
    "CI95_upper", "sensitivity", "specificity"
  )
  required_coordinates <- c("comparison", "false_positive_rate", "sensitivity")

  stopifnot(
    identical(names(summary), required_summary),
    identical(names(coordinates), required_coordinates),
    nrow(summary) == length(groups),
    setequal(summary$comparison, comparisons),
    setequal(coordinates$comparison, comparisons),
    length(unique(summary$panel)) == 1L,
    summary$panel[1L] %in% c("a", "b")
  )
  list(summary = summary, coordinates = coordinates, groups = groups)
}

draw_panel <- function(index, group, metrics, curve) {
  left <- 61 + (index - 1) * 257
  right <- left + 192
  top <- 40
  bottom <- top + 192
  midpoint <- (left + right) / 2
  x <- function(value) left + value * (right - left)
  y <- function(value) bottom - value * (bottom - top)

  fpr <- curve$false_positive_rate
  sensitivity <- curve$sensitivity
  auc_from_curve <- sum(diff(fpr) * (head(sensitivity, -1L) + tail(sensitivity, -1L)) / 2)
  stopifnot(
    length(fpr) >= 2L, all(is.finite(fpr)), all(is.finite(sensitivity)),
    all(diff(fpr) >= 0), all(diff(sensitivity) >= 0),
    fpr[1L] == 0, sensitivity[1L] == 0, tail(fpr, 1L) == 1, tail(sensitivity, 1L) == 1,
    abs(auc_from_curve - metrics$AUC) < 1e-10
  )

  content <- c(
    sprintf('<g id="ROC_AAD_vs_%s">', group),
    svg_text(midpoint, top - 12, paste("AAD VS.", group), 13, "middle", "bold"),
    svg_line(left, top, left, bottom), svg_line(left, bottom, right, bottom),
    svg_line(right, bottom, right, top), svg_line(right, top, left, top)
  )
  for (tick in seq(0, 1, 0.2)) {
    content <- c(
      content, svg_line(x(tick), bottom, x(tick), bottom + 5),
      svg_text(x(tick), bottom + 18, sprintf("%.1f", tick), 10, "middle"),
      svg_line(left - 5, y(tick), left, y(tick)),
      svg_text(left - 9, y(tick) + 3.5, sprintf("%.1f", tick), 10, "end")
    )
  }

  axis_y <- (top + bottom) / 2
  rotation <- sprintf('transform="rotate(-90 %s %s)"', number(left - 39), number(axis_y))
  content <- c(
    content, svg_curve(x(fpr), y(sensitivity)),
    svg_text(midpoint, bottom + 42, "1-Specificity", 12, "middle", "bold"),
    svg_text(left - 39, axis_y, "Sensitivity", 12, "middle", "bold", extra = rotation)
  )
  labels <- c(
    sprintf("AUC = %.3f", metrics$AUC),
    sprintf("(95%% CI: %.3f-%.3f)", metrics$CI95_lower, metrics$CI95_upper),
    sprintf("Sensitivity = %.1f%%", 100 * metrics$sensitivity),
    sprintf("Specificity = %.1f%%", 100 * metrics$specificity),
    "Covariate-specific cutoff"
  )
  for (i in seq_along(labels)) {
    content <- c(content, svg_text(midpoint, top + 101 + (i - 1L) * 16, labels[i], 11, "middle"))
  }
  c(content, "</g>")
}

plot_roc <- function(directory) {
  inputs <- read_inputs(directory)
  panel <- inputs$summary$panel[1L]
  svg <- c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<svg xmlns="http://www.w3.org/2000/svg" width="790pt" height="288pt" viewBox="0 0 790 288">',
    '<rect width="100%" height="100%" fill="white"/>',
    svg_text(8, 22, panel, 20)
  )

  for (index in seq_along(inputs$groups)) {
    group <- inputs$groups[index]
    comparison <- paste("AAD vs", group)
    metrics <- inputs$summary[inputs$summary$comparison == comparison, , drop = FALSE]
    curve <- inputs$coordinates[inputs$coordinates$comparison == comparison, , drop = FALSE]
    svg <- c(svg, draw_panel(index, group, metrics, curve))
  }

  output <- file.path(directory, paste0("Fig4", panel, "_Covariate_Adjusted_ROC.svg"))
  writeLines(c(svg, "</svg>"), output, useBytes = TRUE)
  message("Saved ", basename(output))
  invisible(output)
}

plot_roc(script_directory())
