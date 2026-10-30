# Supplementary Figure 3b: all 81 candidate features; no VIP filtering.
# Keep this script and FigS3b_data.csv together. In RStudio click Source.
# Base R only; no extra packages or internet access required.
# Output: FigS3b_editable.svg. No filter/export switches or extra output files.
#
# Source: Source data/FIG S3/b/脂质组 代谢组合并   2026.10.16.xlsx,
# worksheet 合并标志物 A3:EX83 (81 candidate features before VIP filtering).
# Source SHA256: 2cdf82401a2739bc7f827fc9f509c1ea21e859d53a69ff3439ab171da2c453bf
#
# One CSV carries all 81 candidates, their platform, feature ID,
# source row, and raw intensities for all 152 discovery observations.
# Sample headers use GROUP|SAMPLE_ID. All 81 candidates are always plotted.
# Current scope: 51 metabolomics-platform and 30 lipidomics-platform features.
# The separate 64-feature set was selected from 175 features for Fig.2d; it is
# not the Fig.S3b input. Do not use the appended 28-row list as the 81-row input.
#
# Columns: Ctrl 29, HP 29, AMI 29, AAD 65; original order within each group.
# Rows are split by arithmetic mean abundance in AAD versus Ctrl: lower above,
# higher below. Within each block, complete-linkage hierarchical clustering
# uses Euclidean distances between row Z-score profiles over all 152 samples.
# A stable alphabetical input order resolves algorithmic ties reproducibly.
# The specified 3S feature is placed last by rotating the children along its
# path in the higher-block dendrogram. Merge heights, clusters, and distances
# are unchanged; it is not detached from its cluster or excluded from fitting.
# No column clustering. All samples retained, including H031 and H048.
# Recorded zeros remain unchanged nondetection codes; no imputation or log.
# Per-row Z = (abundance - row mean) / sample SD, using all 152 observations.
# Colour display limits are -1 to +1; full raw values are not truncated.
# LPE-prefixed labels are bold dark green. No additional feature filtering.
# Candidate-feature displays are not independent diagnostic validation.
# Tested with R 4.6.1. SVG contains complete native Arial text labels.

args <- commandArgs(trailingOnly=TRUE)
script_file <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
if(is.null(script_file)) {
  f <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
  if(length(f))script_file <- sub('^--file=','',f[1])
}
out <- if(length(args))args[1] else if(!is.null(script_file))dirname(normalizePath(script_file,winslash='/',mustWork=TRUE)) else getwd()
d <- read.csv(file.path(out,'FigS3b_data.csv'),fileEncoding='UTF-8-BOM',stringsAsFactors=FALSE,check.names=FALSE)
metadata <- c('feature_id','platform','feature_name','source_row')
stopifnot(identical(names(d)[seq_along(metadata)],metadata))
headers <- names(d)[-seq_along(metadata)]
stopifnot(length(headers)==152,all(grepl('^[A-Z]+[|].+',headers)))
samples <- data.frame(sample_id=sub('^[^|]+[|]','',headers),group=sub('[|].*$','',headers),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(samples$sample_id),identical(unique(samples$group),c('CTRL','HP','AMI','AAD')))
stopifnot(identical(as.integer(table(factor(samples$group,levels=c('CTRL','HP','AMI','AAD')))),c(29L,29L,29L,65L)))
keys <- paste(d$platform,d$feature_id,sep=':')
stopifnot(nrow(d)==81,!anyDuplicated(keys),!anyNA(d$feature_name),all(nzchar(d$feature_name)))
d <- d[order(tolower(d$feature_name),d$platform,d$feature_id,method='radix'),,drop=FALSE]
X <- as.matrix(d[,headers,drop=FALSE]);storage.mode(X) <- 'double'
stopifnot(nrow(d)>1,all(is.finite(X)),all(X>=0),all(apply(X,1,sd)>0))
Z <- sweep(sweep(X,1,rowMeans(X),'-'),1,apply(X,1,sd),'/')
stopifnot(max(abs(rowMeans(Z)))<1e-12,max(abs(apply(Z,1,sd)-1))<1e-12)
# Split on raw mean abundance; row Z-scores are used only for distance.
mean_aad <- rowMeans(X[,samples$group=='AAD',drop=FALSE])
mean_ctrl <- rowMeans(X[,samples$group=='CTRL',drop=FALSE])
stopifnot(all(mean_aad!=mean_ctrl))
lower <- which(mean_aad<mean_ctrl)
higher <- which(mean_aad>mean_ctrl)
last_name <- '3S_2349_Tetrahydro_1h_pyrido_34b_indole_carboxylic_acid'
last_index <- which(d$feature_name==last_name)
stopifnot(length(last_index)==1L,last_index%in%higher)
cluster_block <- function(ix,anchor_last=integer()) {
  if(length(ix)<2L)return(ix)
  hc <- hclust(dist(Z[ix,,drop=FALSE],method='euclidean'),method='complete')
  if(length(anchor_last)) {
    leaf <- match(anchor_last,ix)
    contains <- logical(nrow(hc$merge))
    for(k in seq_len(nrow(hc$merge))) {
      children <- hc$merge[k,]
      hit <- vapply(children,function(ch)if(ch<0L)(-ch==leaf)else contains[ch],logical(1))
      if(hit[1])hc$merge[k,] <- rev(children)
      contains[k] <- any(hit)
    }
    leaves <- function(node) {
      if(node<0L)return(-node)
      c(leaves(hc$merge[node,1]),leaves(hc$merge[node,2]))
    }
    hc$order <- leaves(nrow(hc$merge))
    stopifnot(tail(ix[hc$order],1)==anchor_last)
  }
  ix[hc$order]
}
ord_lower <- cluster_block(lower)
ord_higher <- cluster_block(higher,last_index)
ord <- c(ord_lower,ord_higher)
stopifnot(length(ord)==nrow(d),!anyDuplicated(ord))
d <- d[ord,,drop=FALSE];X <- X[ord,,drop=FALSE];Z <- Z[ord,,drop=FALSE]
stopifnot(tail(d$feature_name,1)==last_name)
n_lower <- length(lower);n_higher <- length(higher)
n <- nrow(d)
palette <- colorRampPalette(c('#313695','#4575B4','#74ADD1','#ABD9E9','#E0F3F8','#FFFFBF','#FEE090','#FDAE61','#F46D43','#D73027','#A50026'))(513)
fill_z <- function(z)palette[pmax(1,pmin(513,round((pmax(-1,pmin(1,z))+1)*256)+1))]
group_colours <- c(CTRL='#CA8B51',HP='#C9AD4E',AMI='#577CA2',AAD='#567D47')
group_labels <- c(CTRL='Ctrl',HP='HP',AMI='AMI',AAD='AAD')
W <- 9.5
row_height <- .11
direction_gap <- .14
H <- .66+n*row_height+direction_gap+.25
left <- .80;right <- 5.50;gap <- .035
cw <- (right-left-3*gap)/152
xs <- numeric(152);spans <- list();cursor <- left
for(gr in names(group_colours)) {
  ix <- which(samples$group==gr);xs[ix] <- cursor+(seq_along(ix)-1)*cw
  spans[[gr]] <- c(cursor,cursor+length(ix)*cw);cursor <- cursor+length(ix)*cw+gap
}
y <- H-.66-(seq_len(n)-.5)*row_height
y[(n_lower+1L):n] <- y[(n_lower+1L):n]-direction_gap
is_lpe <- startsWith(d$feature_name,'LPE')
feature_title <- paste(n,'candidate features; clustered within AAD/Ctrl direction blocks')
xml <- function(s) {
  s <- gsub('&','&amp;',s,fixed=TRUE);s <- gsub('<','&lt;',s,fixed=TRUE)
  gsub('>','&gt;',s,fixed=TRUE)
}
text_svg <- function(x,yy,label,size=10,bold=FALSE,colour='#17212B',anchor='start')sprintf('<text x="%.5f" y="%.5f" font-family="Arial" font-size="%.3f" font-weight="%s" fill="%s" text-anchor="%s">%s</text>',x*72,(H-yy)*72+size*.34,size,if(bold)'bold'else'normal',colour,anchor,xml(label))
rect_svg <- function(x,y0,w,h,colour)sprintf('<rect x="%.5f" y="%.5f" width="%.5f" height="%.5f" fill="%s"/>',x*72,(H-y0-h)*72,w*72,h*72,colour)
svg <- c(sprintf('<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%.4fpt" height="%.4fpt" viewBox="0 0 %.4f %.4f">',W*72,H*72,W*72,H*72),sprintf('<rect width="%.5f" height="%.5f" fill="#FFFFFF"/>',W*72,H*72),'<g id="heatmap_cells" shape-rendering="crispEdges">')
svg <- c(svg,rect_svg(rep(xs,times=n),rep(y-row_height/2,each=152),cw,row_height,fill_z(as.vector(t(Z)))),'</g>','<g id="sample_groups">')
for(gr in names(spans)) {
  sp <- spans[[gr]];svg <- c(svg,rect_svg(sp[1],H-.58,sp[2]-sp[1],.035,group_colours[gr]),text_svg(mean(sp),H-.34,group_labels[gr],12.2,TRUE,anchor='middle'))
}
svg <- c(svg,text_svg(.08,H-.14,'b',20.5),text_svg(left,H-.10,feature_title,8),'</g>','<g id="feature_labels">')
for(i in seq_len(n))svg <- c(svg,text_svg(right+.07,y[i],d$feature_name[i],7.2,is_lpe[i],if(is_lpe[i])'#3B6F37'else'#17212B'))
svg <- c(svg,'</g>','<g id="direction_blocks" shape-rendering="crispEdges">')
direction_colours <- c(rep('#577CA2',n_lower),rep('#B6443B',n_higher))
svg <- c(svg,rect_svg(rep(left-.095,n),y-row_height/2,.045,row_height,direction_colours),'</g>','<g id="colour_legend" shape-rendering="crispEdges">')
top <- H-1.03;bottom <- top-1.58;edges <- seq(bottom,top,length.out=514)
svg <- c(svg,rect_svg(.13,edges[-514],.14,diff(edges),palette),text_svg(.13,top+.16,'Row Z-score',6.8))
for(v in c(-1,-.5,0,.5,1))svg <- c(svg,text_svg(.32,bottom+(v+1)/2*(top-bottom),as.character(v),7.2))
svg <- c(svg,rect_svg(.10,bottom-.34,.07,.07,'#577CA2'),text_svg(.22,bottom-.305,'AAD < Ctrl',5.9),rect_svg(.10,bottom-.50,.07,.07,'#B6443B'),text_svg(.22,bottom-.465,'AAD > Ctrl',5.9),'</g>','</svg>')
writeLines(enc2utf8(svg),file.path(out,'FigS3b_editable.svg'),useBytes=TRUE)
cat(sprintf('All %d candidates (%d metabolomics; %d lipidomics); no VIP filtering.\n',n,sum(d$platform=='met'),sum(d$platform=='lip')))
cat(sprintf('%d features x 152 samples; %d cells; %d LPE feature labels.\n',n,n*152,sum(is_lpe)))

cat(sprintf('AAD vs Ctrl: %d lower, %d higher. Complete-linkage / Euclidean row clustering within blocks.\n',n_lower,n_higher))
cat('Last feature:',tail(d$feature_name,1),'\n')
