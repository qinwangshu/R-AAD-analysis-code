# Fig. 2c: standalone base-R script. No additional R packages required.
# RStudio: open this script and click Source. Rscript: Rscript Fig2c_plot.R
# Optional first command-line argument: directory containing Fig2c_all_pathways.csv.
args <- commandArgs(trailingOnly=TRUE)
script_file <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if(is.null(script_file)) {
 file_arg <- grep('^--file=', commandArgs(trailingOnly=FALSE), value=TRUE)
 if(length(file_arg)) script_file <- sub('^--file=', '', file_arg[1])
}
out <- if(length(args)) args[1] else if(!is.null(script_file)) dirname(normalizePath(script_file, winslash='/', mustWork=TRUE)) else getwd()
if(!file.exists(file.path(out,'Fig2c_all_pathways.csv'))) stop('Place Fig2c_all_pathways.csv alongside this R script, or supply its directory as the first argument.')
d <- read.csv(file.path(out,'Fig2c_all_pathways.csv'),fileEncoding='UTF-8-BOM',stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(d)==129, all(d$Hits==as.integer(d$Hits)),
          all(is.finite(d$Raw_p)), all(d$Raw_p>0 & d$Raw_p<=1),
          !anyDuplicated(paste(d$pathway, d$comparison)))
alpha <- 0.05
minimum_p <- tapply(d$Raw_p, d$pathway, min)
selected <- names(minimum_p)[minimum_p < alpha]
selected <- selected[order(minimum_p[selected], selected)]
d <- d[d$pathway %in% selected, , drop=FALSE]
d$row_order <- match(d$pathway, selected)
d$significant_Raw_p <- d$Raw_p < alpha
group_order <- c('AAD VS CTRL','HP VS CTRL','AMI VS CTRL')
d <- d[order(d$row_order, match(d$comparison, group_order)), ]
d <- d[, c('pathway','comparison','row_order','Hits','Raw_p','log2_P','FDR','significant_Raw_p')]
write.csv(d, file.path(out, 'Fig2c_plot_data.csv'), row.names=FALSE, fileEncoding='UTF-8')
stopifnot(length(selected)==22, nrow(d)==59)

groups <- c('AAD VS CTRL','HP VS CTRL','AMI VS CTRL')
paths <- unique(d$pathway[order(d$row_order)])
n <- length(paths);stopifnot(n==22)
palette <- colorRampPalette(c('#4D5B89','#7FA3BC','#C3C6B4','#F2DB94','#E7875B','#BD513F','#902D31'))(512)
lower <- -10;upper <- 0;neutral <- log2(.05)
stopifnot(all(abs(d$log2_P-log2(d$Raw_p))<1e-12))
colour <- function(p){
 value <- pmax(lower,pmin(upper,log2(p)))
 scaled <- ifelse(value>=neutral,.5*(upper-value)/(upper-neutral),.5+.5*(neutral-value)/(neutral-lower))
 index <- round(scaled*511)+1
 palette[pmax(1,pmin(512,index))]
}
size <- function(h).86*sqrt(h)
draw <- function(){
 par(family='sans',bg='white',xaxs='i',yaxs='i')
 par(fig=c(0,5.71/6.85,0,1),mai=c(.38,3.87,1.18,.09))
 plot.new();plot.window(xlim=c(.48,3.52),ylim=c(.35,n+.65))
 rect(.5,.5,3.5,n+.5,col='#FAFAFA',border=NA)
 abline(h=seq(.5,n+.5,by=1),col='#E8EBEE',lwd=.42)
 abline(v=1:3,col='#D9DEE4',lwd=.5)
 axis(2,at=n:1,labels=paths,las=1,tick=FALSE,line=-.34,cex.axis=.95,col.axis='#17212B',hadj=1)
 axis(3,at=1:3,labels=c('AAD vs CTRL','HP vs CTRL','AMI vs CTRL'),las=2,tick=FALSE,line=-.35,cex.axis=.9,col.axis='#17212B',font=2)
 points(match(d$comparison,groups),n+1-d$row_order,pch=21,bg=colour(d$Raw_p),col='#5A6066',lwd=.30,cex=size(d$Hits))
 box(col='#6E747A',lwd=.55)
 par(fig=c(5.786/6.85,1,.02,.97),mai=c(0,0,0,.07),new=TRUE)
 plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
 text(.06,.975,'Hits',adj=c(0,.5),font=2,cex=.95)
 yy <- seq(.935,.725,length.out=7)
 points(rep(.21,7),yy,pch=21,bg='#8A969F',col='#5A6066',lwd=.3,cex=size(1:7))
 text(.5,yy,labels=1:7,adj=c(0,.5),cex=.88)
 text(.06,.655,'log2(P)',adj=c(0,.5),font=2,cex=.95)
 text(.06,.625,'(raw P)',adj=c(0,.5),cex=.70)
 y0 <- .33;y1 <- .585;cuts <- seq(y0,y1,length.out=513)
 barvalues <- seq(upper,lower,length.out=512)
 barcolours <- colour(2^barvalues)
 for(i in 1:512)rect(.09,cuts[i],.31,cuts[i+1],col=barcolours[i],border=NA)
 rect(.09,y0,.31,y1,border='#6E747A',lwd=.4)
 breaks <- c(0,-2,neutral,-6,-8,lower)
 y <- y0+(-breaks)/(-lower)*(y1-y0)
 segments(.31,y,.36,y,col='#555555',lwd=.4)
 text(.42,y,labels=c('0','-2','-4.32','-6','-8','<= -10'),adj=c(0,.5),cex=.74)
 text(.06,.295,'Yellow: P = 0.05',adj=c(0,.5),cex=.62)
 par(fig=c(0,1,0,1),mai=c(0,0,0,0),new=TRUE)
 plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
 text(.023,.982,'c',adj=c(0,1),font=2,cex=1.42)
}
svg(file.path(out,'Fig2c_22pathways.svg'),width=6.85,height=5.6,pointsize=10,family='sans',bg='white')
draw();dev.off()
png(file.path(out,'Fig2c_22pathways_600dpi.png'),width=6.85,height=5.6,units='in',res=600,pointsize=10,bg='white',type='cairo')
draw();dev.off()


# Native SVG text export: retain vector geometry, replace outlined glyph runs
# with real text elements using each original character's position.
# The known label order is checked against all figure text runs.
svg_file <- file.path(out,'Fig2c_22pathways.svg')
z <- readLines(svg_file, encoding='UTF-8', warn=FALSE)
labels <- c(rev(paths), 'AAD vs CTRL','HP vs CTRL','AMI vs CTRL',
            'Hits', as.character(1:7), 'log2(P)', '(raw P)',
            '0','-2','-4.32','-6','-8','<= -10','Yellow: P = 0.05','c')
font_sizes <- c(rep(9.5,n),rep(9,3),9.5,rep(8.8,7),9.5,7,rep(7.4,6),6.2,14.2)
bold <- c(rep(FALSE,n),rep(TRUE,4),rep(FALSE,7),TRUE,FALSE,rep(FALSE,7),TRUE)
stopifnot(length(labels)==n+21,length(font_sizes)==n+21,length(bold)==n+21)
xml_escape <- function(x) {
 x <- gsub('&','&amp;',x,fixed=TRUE)
 x <- gsub('<','&lt;',x,fixed=TRUE)
 gsub('>','&gt;',x,fixed=TRUE)
}
get_attr <- function(lines,name) {
 sub(paste0('.* ',name,'="([^"]+)".*'),'\\1',lines)
}
result <- character();i <- 1L;k <- 0L
while(i<=length(z)) {
 if(grepl('^<g fill=',z[i]) && i<length(z) && grepl('^<use ',z[i+1L])) {
  j <- i+1L
  while(j<=length(z) && z[j]!='</g>') j <- j+1L
  use <- z[(i+1L):(j-1L)]
  k <- k+1L
  stopifnot(k<=length(labels),length(use)==nchar(labels[k],type='chars'))
  x <- as.numeric(get_attr(use,'x'));y <- as.numeric(get_attr(use,'y'))
  vertical <- k %in% (n+1):(n+3)
  if(vertical) {
   attrs <- sprintf('transform="translate(%s,%s) rotate(-90)" x="%s" y="0"',x[1],y[1],paste(sprintf('%.6f',y[1]-y),collapse=' '))
  } else {
   stopifnot(max(y)-min(y)<0.0001)
   attrs <- sprintf('x="%s" y="%s"',paste(sprintf('%.6f',x),collapse=' '),y[1])
  }
  result <- c(result,z[i],sprintf('<text %s font-family="Arial, Helvetica, sans-serif" font-size="%s" font-weight="%s" xml:space="preserve">%s</text>',attrs,font_sizes[k],if(bold[k])'bold' else 'normal',xml_escape(labels[k])),'</g>')
  i <- j+1L
 } else {result <- c(result,z[i]);i <- i+1L}
}
stopifnot(k==n+21,!any(grepl('^<use ',result)))
writeLines(enc2utf8(result),file.path(out,'Fig2c_22pathways_editable.svg'),useBytes=TRUE)
pdf(file.path(out,'Fig2c_22pathways.pdf'),width=6.85,height=5.6,pointsize=10,family='Helvetica',bg='white',useDingbats=FALSE)
draw();dev.off()
writeLines(capture.output(sessionInfo()), file.path(out,'R_sessionInfo.txt'))
cat(sprintf('Retained %d pathways; %d bubbles; %d blank positions.\n', n, nrow(d), 3*n-nrow(d)))
cat('Created editable SVG, outlined SVG, PDF and 600 dpi PNG.\n')


# Illustrator-friendly SVG export.
# Each label is a single normal text object with a single x coordinate.
# Cairo's 26 grid-line clip paths and unused outlined glyph definitions are removed.
ai_lines <- readLines(file.path(out,'Fig2c_22pathways_editable.svg'), encoding='UTF-8', warn=FALSE)
defs_start <- grep('^<defs>', ai_lines)
defs_end <- grep('^</defs>', ai_lines)
stopifnot(length(defs_start)==1L,length(defs_end)==1L,defs_end>defs_start)
ai_lines <- ai_lines[-seq(defs_start,defs_end)]
ai_lines <- sub('<svg ', '<svg version="1.1" ', ai_lines, fixed=TRUE)
ai_lines <- gsub(' clip-path="[^"]+"', '', ai_lines)
label_count <- 0L
for (ii in seq_along(ai_lines)) {
 if (grepl('<text ',ai_lines[ii],fixed=TRUE)) {
  old_x <- regmatches(ai_lines[ii],regexpr(' x="[^"]+"',ai_lines[ii]))
  first_x <- strsplit(gsub('^ x="|"$', '', old_x),' ',fixed=TRUE)[[1]][1]
  ai_lines[ii] <- sub(' x="[^"]+"',paste0(' x="',first_x,'"'),ai_lines[ii])
  ai_lines[ii] <- sub('font-family="[^"]+"','font-family="Arial"',ai_lines[ii])
  label_count <- label_count+1L
 }
}
stopifnot(label_count==43L,!any(grepl('clip-path|<clipPath|<use |<image',ai_lines)))
background <- grep('^<rect x="-',ai_lines)
if(length(background)) ai_lines[background] <- '<rect x="0" y="0" width="493" height="403" fill="#FFFFFF"/>'
writeLines(enc2utf8(ai_lines),file.path(out,'Fig2c_22pathways_AI_text.svg'),useBytes=TRUE)
file.copy(file.path(out,'Fig2c_22pathways.pdf'),file.path(out,'Fig2c_22pathways_AI_text.pdf'),overwrite=TRUE)
cat('Created Illustrator-friendly SVG: 43 complete text objects, no clipping paths.\n')
