# Fig. 3A: exact proportional markers on a white-background lipid map.
# Base R only; no external packages needed for drawing.
# Run from this folder: source("plot_map.R", encoding = "UTF-8")
# Or: Rscript plot_map.R path/to/plot_data.csv
# The SVG contains editable text, lines, circles and mathematically exact sectors.

draw_map <- function(data_file, output_file = file.path(dirname(data_file), "Fig3A_R.svg")) {
  d <- read.csv(data_file, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  stopifnot(nrow(d) == 21L, !anyDuplicated(d$lipid_class), sum(d$n_features) == 466L)
  stopifnot(all(d$n_FDR_lower >= 0), all(d$n_FDR_higher >= 0))
  stopifnot(all(d$n_FDR_lower + d$n_FDR_higher <= d$n_features))
  # Recalculate proportions from INTEGER COUNTS; do not trust rounded percentages.
  nsig <- d$n_FDR_lower + d$n_FDR_higher
  red <- ifelse(nsig > 0, d$n_FDR_higher / nsig, NA_real_)
  blue <- ifelse(nsig > 0, d$n_FDR_lower / nsig, NA_real_)
  stopifnot(all(abs(red[nsig > 0] + blue[nsig > 0] - 1) < 1e-12))
  elements <- character()
  add <- function(s) elements <<- c(elements, s)
  esc <- function(s) {
    s <- gsub("&", "&amp;", s, fixed=TRUE)
    s <- gsub("<", "&lt;", s, fixed=TRUE)
    gsub(">", "&gt;", s, fixed=TRUE)
  }
  txt <- function(x, y, text, size=28, bold=FALSE, anchor="middle") {
    add(sprintf('<text x="%s" y="%s" font-family="Times New Roman, Liberation Serif, serif" font-size="%s" font-weight="%s" text-anchor="%s" fill="#111111">%s</text>',
      x,y,size,if(bold) "bold" else "normal",anchor,esc(text)))
  }
  pie <- function(x,y,r,pred) {
    fill <- if(is.na(pred)) "#A3A8AD" else if(pred==1) "#DF5648" else "#079AD4"
    add(sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s"/>',x,y,r,fill))
    if(!is.na(pred) && pred>0 && pred<1) {
      # Red starts at 12 o'clock and sweeps CLOCKWISE.
      theta <- -pi/2 + 2*pi*pred
      xe <- x+r*cos(theta); ye <- y+r*sin(theta)
      large_arc <- as.integer(pred>0.5)
      add(sprintf('<path d="M %.6f %.6f L %.6f %.6f A %.6f %.6f 0 %d 1 %.6f %.6f Z" fill="#DF5648"/>',
        x,y,x,y-r,r,r,large_arc,xe,ye))
    }
    add(sprintf('<circle cx="%s" cy="%s" r="%s" fill="none" stroke="#111111" stroke-width="1.7"/>',x,y,r))
  }
  labels <- list(
  list(20, 55, "a", 46, TRUE, "start"),
  list(887, 59, "Lipid alterations in AAD vs CTRL", 50, TRUE, "middle"),
  list(887, 102, "Discovery lipidomics: 466 measured features", 30, FALSE, "middle"),
  list(55, 168, "Fatty-acyl-related lipids", 31, TRUE, "start"),
  list(635, 173, "Glycerophospholipids", 31, TRUE, "start"),
  list(1444, 180, "Sphingolipids", 31, TRUE, "start"),
  list(286, 432, "Acylglycerols", 31, TRUE, "start"),
  list(95, 681, "Other annotated classes", 31, TRUE, "start"),
  list(285, 225, "Acyl-CoA", 26, FALSE, "middle"),
  list(732, 262, "Glycerol-3-P", 26, FALSE, "middle"),
  list(904, 262, "LPA", 26, FALSE, "middle"),
  list(1027, 262, "PA", 26, FALSE, "middle"),
  list(1213, 262, "CDP-DAG", 26, FALSE, "middle")
  )
  paths <- list(
  list("M128 216 H210", FALSE),
  list("M366 216 H438", FALSE),
  list("M815 252 H863", FALSE),
  list("M943 252 H993", FALSE),
  list("M1060 252 H1138", FALSE),
  list("M1200 274 L1164 329", FALSE),
  list("M1232 274 L1271 329", FALSE),
  list("M1279 412 V451", FALSE),
  list("M1244 476 H1300", TRUE),
  list("M1014 274 C918 390 757 463 580 528", FALSE),
  list("M581 541 L907 514", FALSE),
  list("M581 545 L914 598", FALSE),
  list("M405 538 H508", TRUE),
  list("M986 510 H1078", TRUE),
  list("M986 596 H1078", TRUE),
  list("M1542 293 H1645", TRUE),
  list("M1495 349 V415", FALSE),
  list("M1495 502 V563", FALSE)
  )
  add('<svg xmlns="http://www.w3.org/2000/svg" width="300mm" height="151.35mm" viewBox="0 0 1774 895">')
  add('<title>Fig. 3A: lipid alterations in AAD vs CTRL</title>')
  add('<defs><marker id="end" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto-start-reverse"><path d="M0 0 L7 3.5 L0 7 Z" fill="#111111"/></marker></defs>')
  add('<rect width="1774" height="895" fill="#FFFFFF"/>')
  for(p in paths) add(sprintf('<path d="%s" fill="none" stroke="#111111" stroke-width="2" marker-end="url(#end)"%s/>',
    p[[1]],if(p[[2]]) ' marker-start="url(#end)"' else ""))
  for(z in labels) do.call(txt,z)
  for(i in seq_len(nrow(d))) {
    txt(d$x[i],d$y[i]-24,d$lipid_class[i])
    pie(d$x[i],d$y[i],d$r[i],red[i])
  }
  add('<rect x="170" y="789" width="1440" height="53" rx="7" fill="#FFFFFF" stroke="#111111" stroke-width="1.5"/>')
  pie(210,815,17,0);txt(240,823,"Decreased species",21,FALSE,"start")
  pie(470,815,17,1);txt(500,823,"Increased species",21,FALSE,"start")
  pie(739,815,17,.25);txt(769,823,"Proportions among significant species",21,FALSE,"start")
  pie(1317,815,17,NA_real_);txt(1347,823,"No significant species",21,FALSE,"start")
  txt(887,866,"Sector areas: proportions of significantly altered species (global FDR < 0.05); nonsignificant species excluded.",20)
  txt(887,889,"Arrows indicate biochemical relationships, not measured flux.",20)
  add("</svg>")
  writeLines(enc2utf8(elements),output_file,useBytes=TRUE)
  message("Saved SVG: ",normalizePath(output_file))
  invisible(d)
}
args <- commandArgs(trailingOnly=TRUE)
data_file <- if(length(args)) args[1] else "plot_data.csv"
if(!file.exists(data_file)) data_file <- file.choose()
draw_map(data_file)

