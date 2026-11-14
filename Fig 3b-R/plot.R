plot_fig3b <- function(input = "plot_data.csv", output = "Fig3b_editable.svg") {
  data <- read.csv(input, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  classes <- c("LPE", "LPC", "PE", "PC", "PI", "PG")
  effects <- c("log2FC_AAD_vs_HP", "log2FC_AAD_vs_AMI")
  q_values <- c("q_AAD_vs_HP", "q_AAD_vs_AMI")
  required <- c("feature_ID", "lipid_class", effects, q_values)
  stopifnot(all(required %in% names(data)), !anyNA(data[required]),
            !anyDuplicated(data$feature_ID), setequal(data$lipid_class, classes),
            all(is.finite(as.matrix(data[effects]))),
            all(as.matrix(data[q_values]) >= 0 & as.matrix(data[q_values]) <= 1))

  lower_in_both <- data$log2FC_AAD_vs_HP < 0 & data$q_AAD_vs_HP < 0.05 &
                   data$log2FC_AAD_vs_AMI < 0 & data$q_AAD_vs_AMI < 0.05
  groups <- factor(data$lipid_class, levels = classes)
  total <- as.integer(table(groups))
  lower <- as.integer(table(groups[lower_in_both]))
  proportions <- lower / total

  text_svg <- function(x, y, label, size, anchor, fill) {
    sprintf('<text x="%g" y="%g" font-size="%g" text-anchor="%s" fill="%s">%s</text>',
            x, y, size, anchor, fill, label)
  }
  line_svg <- function(x1, y1, x2, y2, stroke) {
    sprintf('<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="%s"/>',
            x1, y1, x2, y2, stroke)
  }

  left <- 151
  right <- 1361
  top <- 40
  baseline <- 498
  rows <- seq(82.5, 465.5, length.out = length(classes))
  ticks <- seq(0, 100, by = 25)
  tick_x <- left + (right - left) * ticks / 100
  bar_widths <- (right - left) * proportions
  shown <- lower > 0

  svg <- c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<svg xmlns="http://www.w3.org/2000/svg" width="183mm" height="80.372mm" viewBox="0 0 1480 650">',
    '<title>Fig. 3b. Lipid species significantly lower in AAD than in both HP and AMI</title>',
    '<rect width="1480" height="650" fill="#FFFFFF"/>',
    '<g font-family="Arial, Helvetica, sans-serif">',
    line_svg(tick_x, top, tick_x, baseline, "#E4E9ED"),
    sprintf('<rect x="%g" y="%g" width="%g" height="35" fill="#246D8E"/>',
            left, rows[shown] - 17.5, bar_widths[shown]),
    text_svg(129, rows + 12, classes, 33, "end", "#202831"),
    text_svg(left + bar_widths + 18.75, rows + 12, paste0(lower, "/", total),
             33, "start", "#253644"),
    line_svg(left, baseline, right, baseline, "#8D9AA5"),
    line_svg(tick_x, baseline, tick_x, baseline + 9, "#8D9AA5"),
    text_svg(tick_x, 546, paste0(ticks, "%"), 28, "middle", "#485761"),
    text_svg(786, 606, "Proportion of measured species", 30, "middle", "#202831"),
    '</g>',
    '</svg>'
  )
  writeLines(svg, output, useBytes = TRUE)
  invisible(data.frame(lipid_class = classes, lower_in_both = lower, total = total))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    plot_fig3b()
  } else if (length(args) == 2L) {
    plot_fig3b(args[1], args[2])
  } else {
    stop("Usage: Rscript plot.R [input.csv output.svg]")
  }
}
