# Base R only. No external packages are required.
# Set the working directory to this folder, then source("plot.R").
# Alternatively run: Rscript plot.R path/to/plot_data.csv
# Values are read from the supplied, verified plotting table.
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0L) args[1L] else "plot_data.csv"
if (!file.exists(input_file)) input_file <- file.choose()
d <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE,
              fileEncoding = "UTF-8-BOM")
output_file <- file.path(dirname(normalizePath(input_file)), "figure_R.svg")
W <- 1100
INK <- "#202A33"; MUTED <- "#53626E"; GRID <- "#DCE3E8"
BLUE <- "#0072B2"; ORANGE <- "#D55E00"
elements <- character()
esc <- function(s) {
  s <- gsub("&", "&amp;", as.character(s), fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  s <- gsub(">", "&gt;", s, fixed = TRUE)
  gsub('"', "&quot;", s, fixed = TRUE)
}
add <- function(s) elements <<- c(elements, s)
txt <- function(x, y, s, size = 16, color = INK, bold = FALSE, anchor = "start") {
  add(sprintf('<text x="%.8g" y="%.8g" font-family="Arial, sans-serif" font-size="%.8g" font-weight="%s" fill="%s" text-anchor="%s">%s</text>',
              x, y, size, if (bold) "bold" else "normal", color, anchor, esc(s)))
}
rect <- function(x, y, w, h, fill = "#FFFFFF", stroke = "none", sw = 1) {
  add(sprintf('<rect x="%.8g" y="%.8g" width="%.8g" height="%.8g" fill="%s" stroke="%s" stroke-width="%.8g"/>',x,y,w,h,fill,stroke,sw))
}
line <- function(x1,y1,x2,y2,color=GRID,width=1,dashed=FALSE) {
  add(sprintf('<line x1="%.8g" y1="%.8g" x2="%.8g" y2="%.8g" stroke="%s" stroke-width="%.8g"%s/>',x1,y1,x2,y2,color,width,if(dashed) ' stroke-dasharray="6 5"' else ""))
}
dot <- function(x,y,color,square=FALSE,r=4.5) {
  if(square) rect(x-r,y-r,2*r,2*r,color) else
    add(sprintf('<circle cx="%.8g" cy="%.8g" r="%.8g" fill="%s"/>',x,y,r,color))
}
heat <- function(v,lim=1.5) {
  t <- min(abs(v)/lim,1)
  target <- if(v<0) c(87,146,186) else c(207,127,108)
  z <- round(c(250,251,252)*(1-t)+target*t)
  sprintf("#%02x%02x%02x",as.integer(z[1]),as.integer(z[2]),as.integer(z[3]))
}
start <- function(height,tag,title) {
  H <<- height
  rect(0,0,W,H)
  txt(24,40,tag,27,bold=TRUE);txt(67,39,title,23,bold=TRUE)
}
finish <- function() {
  svg <- c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="%.3fmm" viewBox="0 0 %d %d">',180*H/W,W,H),elements,'</svg>')
  writeLines(enc2utf8(svg),output_file,useBytes=TRUE)
  message("Saved: ",output_file)
}
colorbar <- function(y) {
  x <- 752; w <- 298
  for(i in 0:149) rect(x+i*w/150,y,w/150+.15,12,heat(-1.5+3*i/149))
  for(v in c(-1.5,0,1.5)) txt(x+(v+1.5)/3*w,y+31,as.character(v),13,color=MUTED,anchor="middle")
  txt(x-16,y+11,"log2FC",14,color=MUTED,anchor="end")
}
comparisons <- c("AAD vs CTRL","HP vs CTRL","AMI vs CTRL","AAD vs HP","AAD vs AMI")
matrix_header <- function() {
  txt(28,73,"Discovery cohort: AAD 65; CTRL, HP and AMI 29 each",15,color=MUTED)
  txt(441,113,"Disease vs healthy control",15,color=MUTED,anchor="middle")
  txt(897,113,"Direct disease comparisons",15,color=MUTED,anchor="middle")
  line(218,124,675,124);line(733,124,1060,124)
  centers <- c(278,438,598,818,978)
  for(j in seq_along(comparisons)) txt(centers[j],151,comparisons[j],17,bold=TRUE,anchor="middle")
  centers
}

# All seven LPE features in the discovery lipidomics dataset.
# Cell labels are FC and global q; fill color represents log2FC.
# Do NOT redo BH adjustment over only these 35 displayed comparisons.
species <- c("LPE(15:0)","LPE(16:0)","LPE(18:0)","LPE(18:1)","LPE(18:2)","LPE(20:4)","LPE(22:6)")
stopifnot(nrow(d)==35,!anyDuplicated(paste(d$metabolite,d$comparison)))
stopifnot(all(abs(log2(d$FC)-d$log2FC)<1e-10))
start(680,"b","LPE species in the discovery lipidomics dataset")
centers <- matrix_header()
for(i in seq_along(species)) {
  n <- species[i];y <- 174+(i-1)*56
  txt(28,y+29,n,20,bold=n=="LPE(18:1)")
  for(j in seq_along(comparisons)) {
    r <- d[d$metabolite==n & d$comparison==comparisons[j],];stopifnot(nrow(r)==1)
    x <- centers[j];rect(x-74,y,148,49,heat(r$log2FC))
    txt(x,y+21,sprintf("%.2f",r$FC),19,bold=TRUE,anchor="middle")
    txt(x,y+40,sprintf("q = %.2g",r$q_BH_all_2330_tests),13,anchor="middle")
  }
}
colorbar(584)
txt(28,636,"Cell values: mean abundance ratio and global FDR q value; color: log2FC.",15,color=MUTED)
txt(28,660,"Lipidomics annotations are shown as recorded, without assigning sn positions.",15,color=MUTED)
finish()
