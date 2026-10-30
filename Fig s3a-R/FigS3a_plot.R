# Supplementary Figure 3a: minimal reproducible Venn diagram.
# Keep this script and FigS3a_data.csv in one folder. In RStudio click Source.
# Base R only; no package installation or internet access is required.
# Default output: FigS3a_editable.svg. PDF/PNG are optional switches below.
#
# Audit, 2026-10-20:
# Latest source workbook: Source data/FIG S3/A/脂质组 代谢组合并   2026.10.16.xlsx.
# Source worksheet 特征合并 contains the current 175-feature set (A3:A177).
# 合并标志物 A3:A83 contains the current 81 candidates. Its later 28-row list
# is separate and must not be included in the 81-feature count.
# 4-vinylguaiacol sulfate is already at 合并标志物 A10 and Perillyl Acetate at A46.
# LPC(24:0) sn-1 and LPE(18:0) sn-1 are already absent from this 81-feature set.
# All 12,312 abundance values of these 81 records match the current source data.
# Therefore, do not repeat the two-out/two-in replacement on this workbook.
# Workbook SHA256: 2cdf82401a2739bc7f827fc9f509c1ea21e859d53a69ff3439ab171da2c453bf
#
# Set definitions:
# A = the current 175 AAD-specific candidate features, before VIP filtering.
# B = AAD-vs-HP differential features: 316 metabolomics + 280 lipidomics
#     records, integrated by exact feature name to 585 entries.
# C = AAD-vs-AMI differential features: 265 metabolomics + 185 lipidomics
#     records, integrated by exact feature name to 441 entries.
# Raw two-sided equal-variance Student P < 0.05, no additional VIP or FDR filter.
# Quantification source: 数据/human-全组分数据分析重  2026.10.16.xlsx, 原始数据去重;
#                       数据/human-脂质组数据分析新  2026.10.16.xlsx, 原始数据.
# All 152 discovery observations retained; H031 and H048 were not removed.
# A reproduces the existing cross-platform integration/exclusion rules used
# in Fig.2b. Do not additionally remove the four Fig.2b exclusions from B/C.
#
# Current set totals: 175 / 585 / 441.
# Region order: A only, B only, C only, AB only, AC only, BC only, ABC.
# Current region counts: 32 / 244 / 126 / 44 / 18 / 216 / 81.
# If the old figure remains: B total 575 -> 585; B only 242 -> 244;
# C only 134 -> 126; BC only 208 -> 216; all other displayed counts unchanged.
# Plot counts below are calculated from the membership CSV, not hard-coded.
# The CSV has one row per feature in the union (761 rows), with 0/1 set flags
# and source-platform IDs. Triple-positive rows are the 81 candidates.
# This script plots those already audited sets; it does not rerun statistical
# tests or VIP models. Circle areas are schematic, not proportional to counts.
# Native SVG text is editable. Tested with R 4.6.1; SVG font Arial, PDF Helvetica.
EXPORT_PDF <- FALSE
EXPORT_PNG <- FALSE

args <- commandArgs(trailingOnly=TRUE)
script_file <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
if(is.null(script_file)) {
  f <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
  if(length(f)) script_file <- sub('^--file=','',f[1])
}
out <- if(length(args))args[1] else if(!is.null(script_file))dirname(normalizePath(script_file,winslash='/',mustWork=TRUE)) else getwd()
d <- read.csv(file.path(out,'FigS3a_data.csv'),fileEncoding='UTF-8-BOM',stringsAsFactors=FALSE,check.names=FALSE)
columns <- c('AAD_specific','AAD_vs_HP','AAD_vs_AMI')
stopifnot(all(c('feature',columns)%in%names(d)),!anyDuplicated(d$feature),all(nzchar(d$feature)))
flags <- as.matrix(d[,columns]);storage.mode(flags) <- 'double'
stopifnot(!anyNA(flags),all(flags%in%c(0L,1L)),all(rowSums(flags)>0L))
bits <- rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1))
counts <- apply(bits,1,function(b)sum(rowSums(sweep(flags,2,b,'=='))==3L))
totals <- colSums(flags)
stopifnot(sum(counts)==nrow(d),all(colSums(sweep(bits,1,counts,'*'))==totals))

W <- 4.2
xlim <- c(-1.70,1.70);ylim <- c(-1.80,1.90)
H <- W*diff(ylim)/diff(xlim)
centers <- rbind(c(0,.57),c(-.53,-.34),c(.53,-.34))
radius <- .95
fills <- c('#CE4858','#DDD17B','#8FC56B')
pos <- rbind(c(0,1.00),c(-.91,-.57),c(.91,-.57),c(-.47,.23),c(.47,.23),c(0,-.83),c(0,-.15))
labelpos <- rbind(c(0,1.79),c(-.86,-1.46),c(.86,-1.46))
labels <- c('AAD-specific','AAD vs HP','AAD vs AMI')

draw <- function() {
  par(mar=c(0,0,0,0),xaxs='i',yaxs='i',family='sans',bg='white')
  plot.new();plot.window(xlim=xlim,ylim=ylim,asp=1)
  theta <- seq(0,2*pi,length.out=1001)
  for(i in 1:3) polygon(centers[i,1]+radius*cos(theta),centers[i,2]+radius*sin(theta),col=adjustcolor(fills[i],alpha.f=.65),border=NA)
  for(i in 1:3) lines(centers[i,1]+radius*cos(theta),centers[i,2]+radius*sin(theta),col='#333333',lwd=.65,lty='dashed')
  text(pos[,1],pos[,2],counts,cex=1.08,font=2,col='#101010')
  for(i in 1:3) {
    text(labelpos[i,1],labelpos[i,2],labels[i],cex=1.00,font=2)
    text(labelpos[i,1],labelpos[i,2]-.19,paste(totals[i],'features'),cex=.87)
  }
  text(-1.59,1.75,'a',cex=1.55,font=2)
}

xml <- function(x) {
  x <- gsub('&','&amp;',x,fixed=TRUE)
  x <- gsub('<','&lt;',x,fixed=TRUE)
  gsub('>','&gt;',x,fixed=TRUE)
}
sx <- function(x)(x-xlim[1])*W*72/diff(xlim)
sy <- function(y)(ylim[2]-y)*H*72/diff(ylim)
text_svg <- function(x,y,label,size=10,bold=FALSE,anchor='middle')sprintf('<text x="%.5f" y="%.5f" font-family="Arial" font-size="%.2f" font-weight="%s" fill="#101010" text-anchor="%s">%s</text>',sx(x),sy(y)+size*.34,size,if(bold)'bold'else'normal',anchor,xml(label))
svg <- c(sprintf('<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%.5fpt" height="%.5fpt" viewBox="0 0 %.5f %.5f">',W*72,H*72,W*72,H*72),sprintf('<rect width="%.5f" height="%.5f" fill="#FFFFFF"/>',W*72,H*72),'<g id="sets">')
for(i in 1:3)svg <- c(svg,sprintf('<circle cx="%.5f" cy="%.5f" r="%.5f" fill="%s" fill-opacity="0.65"/>',sx(centers[i,1]),sy(centers[i,2]),radius*W*72/diff(xlim),fills[i]))
for(i in 1:3)svg <- c(svg,sprintf('<circle cx="%.5f" cy="%.5f" r="%.5f" fill="none" stroke="#333333" stroke-width="0.5" stroke-dasharray="2.5,2.5"/>',sx(centers[i,1]),sy(centers[i,2]),radius*W*72/diff(xlim)))
svg <- c(svg,'</g>','<g id="region_counts">')
for(i in 1:7)svg <- c(svg,text_svg(pos[i,1],pos[i,2],counts[i],10.8,TRUE))
svg <- c(svg,'</g>','<g id="set_labels">')
for(i in 1:3)svg <- c(svg,text_svg(labelpos[i,1],labelpos[i,2],labels[i],10,TRUE),text_svg(labelpos[i,1],labelpos[i,2]-.19,paste(totals[i],'features'),8.7))
svg <- c(svg,'</g>',text_svg(-1.59,1.75,'a',15.5,TRUE),'</svg>')
writeLines(enc2utf8(svg),file.path(out,'FigS3a_editable.svg'),useBytes=TRUE)
if(EXPORT_PDF) {
  pdf(file.path(out,'FigS3a.pdf'),width=W,height=H,pointsize=10,family='Helvetica',useDingbats=FALSE,bg='white')
  draw();dev.off()
}
if(EXPORT_PNG) {
  png(file.path(out,'FigS3a_600dpi.png'),width=W,height=H,units='in',res=600,pointsize=10,bg='white',type='cairo')
  draw();dev.off()
}
cat('Set totals:',paste(totals,collapse=' / '),'\n')
cat('Seven regions:',paste(counts,collapse=' / '),'\n')
cat('Union:',nrow(d),'features; triple intersection:',counts[7],'features.\n')
