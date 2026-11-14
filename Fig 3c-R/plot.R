args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  args <- c(Sys.getenv("FIG3C_SAMPLES", "plot_data.csv"),
            Sys.getenv("FIG3C_RESULTS", "comparison_results.csv"),
            Sys.getenv("FIG3C_SVG", "Fig3c_boxplots_editable.svg"))
}
if (length(args) != 3L || any(!nzchar(args))) {
  stop("Supply the sample CSV, original Fig. 3c results CSV, and output SVG.")
}

samples <- read.csv(args[1], fileEncoding = "UTF-8-BOM")
results <- read.csv(args[2], fileEncoding = "UTF-8-BOM")
groups <- c("CTRL", "HP", "AMI", "AAD")
features <- unique(results$feature_ID)
comparisons <- c("AAD vs CTRL", "HP vs CTRL", "AMI vs CTRL", "AAD vs HP", "AAD vs AMI")
aad_comparisons <- paste("AAD vs", rev(groups[groups != "AAD"]))
group_sizes <- as.integer(table(factor(samples$group, levels = groups)))
stopifnot(length(features) == 7L, nrow(results) == 35L, nrow(samples) == 152L,
          !anyDuplicated(samples$sample_ID),
          identical(group_sizes, c(29L, 29L, 29L, 65L)),
          all(table(results$feature_ID, results$comparison) == 1),
          setequal(results$comparison, comparisons),
          all(is.finite(as.matrix(samples[features]))),
          all(as.matrix(samples[features]) >= 0),
          all(results$q_BH_all_2330_tests >= 0 & results$q_BH_all_2330_tests <= 1))

values <- function(feature, group) samples[samples$group == group, feature]
for (i in seq_len(nrow(results))) {
  row <- results[i, ]
  comparison_groups <- strsplit(row$comparison, " vs ", fixed = TRUE)[[1]]
  a <- values(row$feature_ID, comparison_groups[1])
  b <- values(row$feature_ID, comparison_groups[2])
  stopifnot(length(a) == row$n_numerator, length(b) == row$n_denominator,
            sum(a == 0) == row$zeros_numerator, sum(b == 0) == row$zeros_denominator,
            abs(mean(a) - row$mean_numerator) < 1e-10,
            abs(mean(b) - row$mean_denominator) < 1e-10)
}

ink <- "#202831"
muted <- "#4D5964"
blue <- "#1E6486"
guide <- "#E5EAEE"
axis_color <- "#A3ADB5"
panel_width <- 248
width <- 1042
height <- 550

label <- function(x, y, content, size = 14, anchor = "start", fill = ink, bold = FALSE) {
  content <- gsub("<", "&lt;", content, fixed = TRUE)
  sprintf('<text x="%g" y="%g" font-size="%g" text-anchor="%s" fill="%s" font-weight="%s">%s</text>',
          x, y, size, anchor, fill, if (bold) "bold" else "normal", content)
}
line <- function(x1, y1, x2, y2, color, weight = 1) {
  sprintf('<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="%s" stroke-width="%g"/>',
          x1, y1, x2, y2, color, weight)
}
rect <- function(x, y, w, h, fill, stroke = "none") {
  sprintf('<rect x="%g" y="%g" width="%g" height="%g" fill="%s" stroke="%s" stroke-width="1.4"/>',
          x, y, w, h, fill, stroke)
}
dot <- function(x, y, fill) {
  sprintf('<circle cx="%g" cy="%g" r="2.25" fill="%s" fill-opacity="0.34"/>', x, y, fill)
}

draw_group <- function(abundance, x, pos_y, is_aad) {
  stats <- boxplot.stats(abundance)$stats
  y <- pos_y(stats)
  edge <- if (is_aad) blue else "#596673"
  box_fill <- if (is_aad) "#D9EBF4" else "#EDF0F2"
  x_jitter <- runif(length(abundance), -12, 12)
  c(line(x, y[1], x, y[2], edge, 1.3),
    line(x, y[4], x, y[5], edge, 1.3),
    line(x - 6, y[1], x + 6, y[1], edge, 1.3),
    line(x - 6, y[5], x + 6, y[5], edge, 1.3),
    rect(x - 13, y[4], 26, y[2] - y[4], box_fill, edge),
    dot(x + x_jitter, pos_y(abundance), edge),
    line(x - 13, y[3], x + 13, y[3], edge, 2))
}

draw_panel <- function(feature, index) {
  x0 <- 28 + (index - 1) %% 4 * panel_width
  y0 <- 13 + (index - 1) %/% 4 * 264
  panel <- results[results$feature_ID == feature, ]
  species <- unique(panel$metabolite)
  abundance <- lapply(groups, function(group) values(feature, group))
  upper <- max(unlist(abundance))
  ticks <- pretty(c(0, upper), n = 4)
  ticks <- ticks[ticks >= 0 & ticks <= upper * 1.1]
  top <- y0 + 94
  bottom <- y0 + 224
  limit <- max(ticks, upper) * 1.03
  pos_y <- function(value) bottom - (bottom - top) * value / limit
  centers <- x0 + c(65, 112, 159, 206)
  art <- c(label(x0 + 128, y0 + 20, species, 17, "middle", bold = TRUE),
           sprintf('<g transform="translate(%g %g) rotate(-90)">', x0 - 1, (top + bottom)/2),
           label(0, 0, "Relative abundance", 14, "middle"), '</g>',
           line(x0 + 40, pos_y(ticks), x0 + 229, pos_y(ticks), guide),
           line(x0 + 40, top, x0 + 40, bottom, axis_color),
           line(x0 + 40, bottom, x0 + 229, bottom, axis_color),
           label(x0 + 34, pos_y(ticks) + 4, format(ticks, trim = TRUE), 12.5, "end", muted))
  for (j in seq_along(groups)) {
    art <- c(art, draw_group(abundance[[j]], centers[j], pos_y, groups[j] == "AAD"),
             label(centers[j], bottom + 19, groups[j], 13, "middle", muted))
  }

  q <- panel$q_BH_all_2330_tests[match(aad_comparisons, panel$comparison)]
  for (k in seq_along(aad_comparisons)) {
    x <- centers[c(length(groups) - k, length(groups))]
    y <- top - 8 - (k - 1) * 18
    mark <- if (q[k] < 0.05) "*" else "ns"
    art <- c(art, line(x[1], y + 4, x[1], y, muted),
             line(x[1], y, x[2], y, muted),
             line(x[2], y, x[2], y + 4, muted),
             label(mean(x), y - 2, mark, 14, "middle", muted, bold = q[k] < 0.05))
  }
  art
}

set.seed(31)
svg <- c('<?xml version="1.0" encoding="UTF-8"?>',
         sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="183mm" height="%gmm" viewBox="0 0 %g %g">',
                 183 * height / width, width, height),
         '<title>Fig. 3c. Relative abundances of seven LPE species</title>',
         rect(0, 0, width, height, "#FFFFFF"),
         '<g font-family="Arial, Helvetica, sans-serif">',
         label(8, 23, "c", 22, bold = TRUE))
for (i in seq_along(features)) {
  svg <- c(svg, draw_panel(features[i], i))
}
svg <- c(svg, '</g>', '</svg>')
writeLines(enc2utf8(svg), args[3], useBytes = TRUE)
