# Figure 2d: minimal reproducible 77-feature heatmap.
# Keep only this script and Fig2d_77_data.csv together; open this file in
# RStudio and click Source. Base R only; no extra packages or internet needed.
# Default output: Fig2d_77_editable.svg. Existing SVG at that name is replaced.
# Set EXPORT_PDF/EXPORT_PNG below to TRUE only when those extra outputs are wanted.
# Data: 43 glycerophospholipid and 34 sphingolipid-related features from the
# existing 175-feature candidate set; no additional VIP or P-value filtering.
# CSV first six columns: IDs, platform, name, subclass, existing VIP, source row.
# Remaining columns: GROUP|SAMPLE_ID; these headers carry group assignments.
# All 152 samples are retained, including H048 and H031. Zeros are unchanged
# nondetection codes. No log or imputation. Row Z uses sample SD (n-1).
# Display limits are -1/+1; underlying values are not truncated.
# Source: Source data/FIG 2/d/脂质组 代谢组合并   2026.10.16.xlsx, sheet 特征合并.
# Source SHA256: 6ade8dad31b69a28d31612a8e7c6c4ab7c515abc1ad7347f26a8a3d2b1864057
# This descriptive class subset is not a reconstructed KEGG pathway-hit list.
# Tested with R 4.6.1. SVG has native Arial text; PDF uses Helvetica.
EXPORT_PDF <- FALSE
EXPORT_PNG <- FALSE
args <- commandArgs(trailingOnly=TRUE)
script_file <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
if(is.null(script_file)) {
  filearg <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
  if(length(filearg)) script_file <- sub('^--file=','',filearg[1])
}
out <- if(length(args)) args[1] else if(!is.null(script_file)) dirname(normalizePath(script_file,winslash='/',mustWork=TRUE)) else getwd()
read_input <- function(name) read.csv(file.path(out,name),fileEncoding='UTF-8-BOM',stringsAsFactors=FALSE,check.names=FALSE)
d <- read_input('Fig2d_77_data.csv')
metadata <- c('feature_id','platform','feature_name','subclass','VIP','source_row')
stopifnot(identical(names(d)[seq_along(metadata)],metadata))
sample_headers <- names(d)[-(seq_along(metadata))]
stopifnot(length(sample_headers)==152,all(grepl('^[A-Z]+[|].+',sample_headers)))
samples <- data.frame(
  sample_id=sub('^[^|]+[|]','',sample_headers),
  group=sub('[|].*$','',sample_headers),stringsAsFactors=FALSE)
ann <- d[,metadata,drop=FALSE]
ann$feature_key <- paste(ann$platform,ann$feature_id,sep=':')
ann$VIP_gt_1 <- ann$VIP > 1
ann$is_LPE <- ann$subclass=='LPE'
classes <- data.frame(
  subclass=c('PC','LPC','PE','LPE','dMePE','Cer','SM','Hex1Cer','Hex2Cer'),
  label=c('PC','LPC','PE','LPE','dMePE','Cer','SM','Hex1Cer','Hex2Cer'),
  colour=c('#789F4C','#97B35D','#5E8C3E','#3B6F37','#A2B97B','#70479A','#6E774D','#489FBF','#70BDD0'),
  lipid_block=c(rep('Glycerophospholipids',5),rep('Sphingolipid-related',4)),
  n_features=c(18,6,8,8,3,9,17,7,1),stringsAsFactors=FALSE)
blocks <- data.frame(
  lipid_block=c('Glycerophospholipids','Sphingolipid-related'),
  label=c('Glycerophospholipids','Sphingolipid-related'),
  colour=c('#789F4C','#587E9D'),n_features=c(43,34),stringsAsFactors=FALSE)
ann$lipid_block <- classes$lipid_block[match(ann$subclass,classes$subclass)]
raw <- data.frame(feature_key=ann$feature_key,d[,sample_headers,drop=FALSE],check.names=FALSE)
names(raw)[-1] <- samples$sample_id
stopifnot(identical(as.integer(table(factor(samples$group,levels=c('CTRL','HP','AMI','AAD')))),c(29L,29L,29L,65L)))
stopifnot(nrow(ann)==77,nrow(samples)==152,!anyDuplicated(ann$feature_key),!anyDuplicated(samples$sample_id))
stopifnot(identical(unique(samples$group),c('CTRL','HP','AMI','AAD')))
palette <- colorRampPalette(c('#313695','#4575B4','#74ADD1','#ABD9E9','#E0F3F8','#FFFFBF','#FEE090','#FDAE61','#F46D43','#D73027','#A50026'))(513)
fill_z <- function(z) palette[pmax(1,pmin(513,round((pmax(-1,pmin(1,z))+1)*256)+1))]
group_colours <- c(CTRL='#CA8B51',HP='#C9AD4E',AMI='#577CA2',AAD='#567D47')
group_labels <- c(CTRL='Ctrl',HP='HP',AMI='AMI',AAD='AAD')
W <- 11.6
row_height <- .09
class_gap <- .055
block_extra_gap <- .095
left <- .80
right <- 5.50
column_gap <- .035
cell_width <- (right-left-3*column_gap)/152
xs <- numeric(152)
cursor <- left
group_spans <- list()
for(g in names(group_colours)) {
  ix <- which(samples$group==g)
  xs[ix] <- cursor+(seq_along(ix)-1)*cell_width
  group_spans[[g]] <- c(cursor,cursor+length(ix)*cell_width)
  cursor <- cursor+length(ix)*cell_width+column_gap
}

layout_rows <- function(a) {
  present <- classes$subclass[classes$subclass %in% a$subclass]
  H <- .66+nrow(a)*row_height+(length(present)-1)*class_gap+block_extra_gap+.25
  yy <- numeric(nrow(a));spans <- list();cursor <- H-.66
  for(cl in present) {
    ix <- which(a$subclass==cl)
    yy[ix] <- cursor-(seq_along(ix)-.5)*row_height
    spans[[cl]] <- c(cursor,cursor-length(ix)*row_height)
    cursor <- cursor-length(ix)*row_height-class_gap
    if(cl=='dMePE') cursor <- cursor-block_extra_gap
  }
  block_spans <- list()
  for(bl in blocks$lipid_block) {
    members <- classes$subclass[classes$lipid_block==bl]
    ss <- spans[intersect(members,names(spans))]
    block_spans[[bl]] <- c(max(vapply(ss,function(x)x[1],numeric(1))),min(vapply(ss,function(x)x[2],numeric(1))))
  }
  list(H=H,y=yy,spans=spans,block_spans=block_spans)
}

draw_plot <- function(a,Z,lay) {
  H <- lay$H;y <- lay$y
  par(mar=c(0,0,0,0),xaxs='i',yaxs='i',family='sans',bg='white',xpd=NA)
  plot.new();plot.window(xlim=c(0,W),ylim=c(0,H))
  colours <- matrix(fill_z(as.vector(Z)),nrow=nrow(Z))
  rect(rep(xs,each=nrow(Z)),rep(y-row_height/2,times=152),rep(xs+cell_width,each=nrow(Z)),rep(y+row_height/2,times=152),col=as.vector(colours),border=NA)
  for(g in names(group_spans)) {
    span <- group_spans[[g]]
    rect(span[1],H-.58,span[2],H-.545,col=group_colours[g],border=NA)
    text(mean(span),H-.34,group_labels[g],cex=1.22,font=2)
  }
  text(.08,H-.14,'d',adj=c(0,.5),cex=2.05)
  text(left,H-.10,paste(nrow(a),'features'),adj=c(0,.5),cex=.80)
  for(i in seq_len(nrow(a))) {
    bold <- isTRUE(a$is_LPE[i])
    text(right+.07,y[i],a$feature_name[i],adj=c(0,.5),cex=.63,font=if(bold)2 else 1,col=if(bold)'#3B6F37' else '#17212B')
  }
  # The annotation strip reports VIP status, not abundance or significance.
  rect(left-.095,y-row_height/2,left-.050,y+row_height/2,col=ifelse(a$VIP_gt_1,'#596B79','#E5E9EC'),border=NA)
  text(left-.072,H-.615,'VIP',cex=.52)
  for(cl in names(lay$spans)) {
    s <- lay$spans[[cl]];mid <- mean(s)
    co <- classes$colour[match(cl,classes$subclass)]
    bx <- 7.65
    segments(bx-.14,s[1]-.012,bx,s[1]-.012,col=co,lwd=.8)
    segments(bx,s[1]-.012,bx,s[2]+.012,col=co,lwd=.8)
    segments(bx-.14,s[2]+.012,bx,s[2]+.012,col=co,lwd=.8)
    segments(bx,mid,bx+.12,mid,col=co,lwd=.8)
    text(bx+.20,mid,paste0(cl,' (',sum(a$subclass==cl),')'),adj=c(0,.5),cex=.80,font=2,col=co)
  }
  for(bl in names(lay$block_spans)) {
    sp <- lay$block_spans[[bl]]; mid <- mean(sp); bx <- 9.10
    co <- blocks$colour[match(bl,blocks$lipid_block)]
    segments(bx-.10,sp[1]-.012,bx,sp[1]-.012,col=co,lwd=.8)
    segments(bx,sp[1]-.012,bx,sp[2]+.012,col=co,lwd=.8)
    segments(bx-.10,sp[2]+.012,bx,sp[2]+.012,col=co,lwd=.8)
    segments(bx,mid,bx+.10,mid,col=co,lwd=.8)
    text(bx+.20,mid,blocks$label[match(bl,blocks$lipid_block)],adj=c(0,.5),cex=.90,font=2,col=co)
  }
  # Fixed display limits; complete Z-scores are preserved in CSV.
  top <- H-1.03;bottom <- top-1.58
  edges <- seq(bottom,top,length.out=514)
  rect(.13,edges[-514],.27,edges[-1],col=palette,border=NA)
  text(.13,top+.16,'Row Z-score',adj=c(0,.5),cex=.68)
  ticks <- c(-1,-.5,0,.5,1)
  text(.32,bottom+(ticks+1)/2*(top-bottom),as.character(ticks),adj=c(0,.5),cex=.72)
  rect(.10,bottom-.34,.17,bottom-.27,col='#596B79',border=NA)
  text(.22,bottom-.305,'VIP > 1',adj=c(0,.5),cex=.59)
  rect(.10,bottom-.50,.17,bottom-.43,col='#E5E9EC',border=NA)
  text(.22,bottom-.465,'VIP <= 1',adj=c(0,.5),cex=.59)
}

xml_escape <- function(s) {
  s <- gsub('&','&amp;',s,fixed=TRUE)
  s <- gsub('<','&lt;',s,fixed=TRUE)
  gsub('>','&gt;',s,fixed=TRUE)
}

write_svg <- function(filename,a,Z,lay) {
  # Ordinary SVG 1.1 text with one x/y position per label; no clips or outlined glyphs.
  H <- lay$H;y <- lay$y
  chunks <- list(sprintf('<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%.4fpt" height="%.4fpt" viewBox="0 0 %.4f %.4f">',W*72,H*72,W*72,H*72))
  add <- function(s) chunks[[length(chunks)+1L]] <<- s
  text_svg <- function(x,yy,label,size=10,bold=FALSE,colour='#17212B',anchor='start') {
    sprintf('<text x="%.5f" y="%.5f" font-family="Arial" font-size="%.3f" font-weight="%s" fill="%s" text-anchor="%s">%s</text>',x*72,(H-yy)*72+size*.34,size,if(bold)'bold' else 'normal',colour,anchor,xml_escape(label))
  }
  rect_svg <- function(x,y0,width,height,colour) sprintf('<rect x="%.5f" y="%.5f" width="%.5f" height="%.5f" fill="%s"/>',x*72,(H-y0-height)*72,width*72,height*72,colour)
  add(sprintf('<rect width="%.5f" height="%.5f" fill="#FFFFFF"/>',W*72,H*72))
  add('<g id="heatmap_cells" shape-rendering="crispEdges">')
  add(rect_svg(rep(xs,times=nrow(Z)),rep(y-row_height/2,each=152),cell_width,row_height,fill_z(as.vector(t(Z)))))
  add('</g>')
  add('<g id="sample_groups">')
  for(g in names(group_spans)) {
    s <- group_spans[[g]]
    add(rect_svg(s[1],H-.58,s[2]-s[1],.035,group_colours[g]))
    add(text_svg(mean(s),H-.34,group_labels[g],12.2,TRUE,anchor='middle'))
  }
  add(text_svg(.08,H-.14,'d',20.5))
  add(text_svg(left,H-.10,paste(nrow(a),'features'),8))
  add('</g>')
  add('<g id="feature_labels">')
  for(i in seq_len(nrow(a))) add(text_svg(right+.07,y[i],a$feature_name[i],6.3,isTRUE(a$is_LPE[i]),if(isTRUE(a$is_LPE[i]))'#3B6F37' else '#17212B'))
  add('</g>')
  add('<g id="VIP_annotation" shape-rendering="crispEdges">')
  add(rect_svg(left-.095,y-row_height/2,.045,row_height,ifelse(a$VIP_gt_1,'#596B79','#E5E9EC')))
  add(text_svg(left-.072,H-.615,'VIP',5.2,anchor='middle'))
  add('</g>')
  add('<g id="chemical_classes">')
  for(cl in names(lay$spans)) {
    s <- lay$spans[[cl]];mid <- mean(s);bx <- 7.65
    co <- classes$colour[match(cl,classes$subclass)]
    add(sprintf('<path d="M %.5f %.5f H %.5f V %.5f H %.5f M %.5f %.5f H %.5f" fill="none" stroke="%s" stroke-width="0.6"/>',(bx-.14)*72,(H-s[1]+.012)*72,bx*72,(H-s[2]-.012)*72,(bx-.14)*72,bx*72,(H-mid)*72,(bx+.12)*72,co))
    add(text_svg(bx+.20,mid,paste0(cl,' (',sum(a$subclass==cl),')'),8.0,TRUE,co))
  }
  add('</g>')
  add('<g id="lipid_blocks">')
  for(bl in names(lay$block_spans)) {
    sp <- lay$block_spans[[bl]]; mid <- mean(sp); bx <- 9.10
    co <- blocks$colour[match(bl,blocks$lipid_block)]
    add(sprintf('<path d="M %.5f %.5f H %.5f V %.5f H %.5f M %.5f %.5f H %.5f" fill="none" stroke="%s" stroke-width="0.6"/>',(bx-.10)*72,(H-sp[1]+.012)*72,bx*72,(H-sp[2]-.012)*72,(bx-.10)*72,bx*72,(H-mid)*72,(bx+.10)*72,co))
    add(text_svg(bx+.20,mid,blocks$label[match(bl,blocks$lipid_block)],9.0,TRUE,co))
  }
  add('</g>')
  top <- H-1.03;bottom <- top-1.58;edges <- seq(bottom,top,length.out=514)
  add('<g id="colour_legend" shape-rendering="crispEdges">')
  add(rect_svg(.13,edges[-514],.14,diff(edges),palette))
  add(text_svg(.13,top+.16,'Row Z-score',6.8))
  for(v in c(-1,-.5,0,.5,1)) add(text_svg(.32,bottom+(v+1)/2*(top-bottom),as.character(v),7.2))
  add(rect_svg(.10,bottom-.34,.07,.07,'#596B79'));add(text_svg(.22,bottom-.305,'VIP > 1',5.9))
  add(rect_svg(.10,bottom-.50,.07,.07,'#E5E9EC'));add(text_svg(.22,bottom-.465,'VIP <= 1',5.9))
  add('</g>\n</svg>')
  writeLines(enc2utf8(unlist(chunks,use.names=FALSE)),filename,useBytes=TRUE)
}

for(n in 77L) {
  stopifnot(nrow(raw)==n,identical(names(raw)[-1],samples$sample_id))
  a <- ann[match(raw$feature_key,ann$feature_key),,drop=FALSE]
  stopifnot(!anyNA(a$feature_key),all(!is.na(match(a$subclass,classes$subclass))))
  stopifnot(identical(as.integer(table(factor(a$subclass,levels=classes$subclass))),as.integer(classes$n_features)))
  ord <- order(match(a$subclass,classes$subclass),a$source_row)
  a <- a[ord,,drop=FALSE]
  X <- as.matrix(raw[ord,-1,drop=FALSE]);storage.mode(X) <- 'double'
  stopifnot(all(is.finite(X)),all(X>=0),all(apply(X,1,sd)>0))
  Z <- sweep(sweep(X,1,rowMeans(X),'-'),1,apply(X,1,sd),'/')
  stopifnot(max(abs(rowMeans(Z)))<1e-12,max(abs(apply(Z,1,sd)-1))<1e-12)
  lay <- layout_rows(a)
  prefix <- file.path(out,sprintf('Fig2d_%d',n))
  if(EXPORT_PDF) {
    pdf(paste0(prefix,'.pdf'),width=W,height=lay$H,pointsize=10,family='Helvetica',bg='white',useDingbats=FALSE)
    draw_plot(a,Z,lay);dev.off()
  }
  if(EXPORT_PNG) {
    png(paste0(prefix,'_600dpi.png'),width=W,height=lay$H,units='in',res=600,pointsize=10,bg='white',type='cairo')
    draw_plot(a,Z,lay);dev.off()
  }
  write_svg(paste0(prefix,'_editable.svg'),a,Z,lay)
  cat(sprintf('%d features: %d samples; figure %.2f x %.3f inches; %d displayed cells.\n',n,ncol(X),W,lay$H,length(X)))
}
