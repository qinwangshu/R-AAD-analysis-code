# Fig. S3c: Mfuzz reanalysis of all 81 current candidates; no VIP filtering.
# Keep this script and FigS3c_data.csv together; Source in RStudio or run:
# Rscript FigS3c_plot.R
# Dependencies (install once):
# install.packages('BiocManager')
# BiocManager::install('Mfuzz', ask=FALSE, update=FALSE)
# Tested with R 4.6.1, Mfuzz 2.72.0, e1071 1.7-17.
#
# Source: FIG S3/b/脂质组 代谢组合并   2026.10.16.xlsx,
# 合并标志物 A3:EX83; SHA256:
# 2cdf82401a2739bc7f827fc9f509c1ea21e859d53a69ff3439ab171da2c453bf
# Data are the same 81 x 152 raw intensities used in current Fig. S3b.
# CTRL29, HP29, AMI29, AAD65; H031 and H048 retained. No log, imputation,
# VIP filtering, or outlier deletion. Recorded nondetection zeros retained.
# Compute arithmetic group means FIRST. Then Mfuzz::standardise each feature
# across the four means using the sample SD (n-1). This matches the verified
# old preprocessing. This differs from averaging individual-sample Z-scores.
#
# Fit all 81 profiles jointly, without pre-splitting into up/down blocks.
# Four clusters retained from original Fig. S3c; c=4 is a fixed choice, not
# a newly demonstrated optimal cluster count. Old m and seed were not saved.
# Re-estimate m using Mfuzz::mestimate. Seed=20261020, 100 random starts,
# max 10000 iterations, relative tolerance 1e-10. Retain lowest objective.
# Re-label solely for reproducible display: centroid peak group (CTRL to AAD),
# then centroid CTRL, HP, AMI, AAD. Cluster numbers do not identify old clusters.
# Highest-membership assignment is exhaustive, even if memberships are low;
# there is no membership cutoff. All four memberships are exported for review.
# The clinical groups are cross-sectional categories, not a time course.
# Curves indicate relative profiles; clusters alone do not prove specificity
# or diagnostic performance. No manuscript file is modified by this script.
# Outputs: FigS3c_results.csv (means, Z, all memberships, cluster),
# FigS3c_editable.svg (native editable text and vector lines).

args <- commandArgs(trailingOnly=TRUE)
script_file <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
if(is.null(script_file)) {
 f <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
 if(length(f))script_file <- sub('^--file=','',f[1])
}
out <- if(length(args))args[1] else if(!is.null(script_file))dirname(normalizePath(script_file,winslash='/',mustWork=TRUE)) else getwd()
if(!requireNamespace('Mfuzz',quietly=TRUE))stop('Install BiocManager and Mfuzz as described at the top of this script.')
suppressPackageStartupMessages(library(Mfuzz))
d <- read.csv(file.path(out,'FigS3c_data.csv'),fileEncoding='UTF-8-BOM',check.names=FALSE,stringsAsFactors=FALSE)
sample_cols <- grep('^(CTRL|HP|AMI|AAD)[|]',names(d),value=TRUE)
groups <- sub('[|].*$','',sample_cols)
stopifnot(nrow(d)==81L,length(sample_cols)==152L,!anyDuplicated(paste(d$platform,d$feature_id)),
 identical(as.integer(table(factor(groups,levels=c('CTRL','HP','AMI','AAD')))),c(29L,29L,29L,65L)))
d <- d[order(tolower(d$feature_name),d$platform,d$feature_id,method='radix'),,drop=FALSE]
X <- as.matrix(d[,sample_cols]);storage.mode(X) <- 'double'
stopifnot(all(is.finite(X)),all(X>=0))
g <- c('CTRL','HP','AMI','AAD')
M <- sapply(g,function(k)rowMeans(X[,groups==k,drop=FALSE]))
rownames(M) <- paste(d$platform,d$feature_id,sep=':')
stopifnot(all(apply(M,1,sd)>0))
es <- Mfuzz::standardise(Biobase::ExpressionSet(assayData=M))
Z <- Biobase::exprs(es)
m <- Mfuzz::mestimate(es)
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(20261020)
fits <- lapply(seq_len(100),function(i)Mfuzz::mfuzz(es,centers=4,m=m,iter.max=10000,control=list(reltol=1e-10)))
objectives <- vapply(fits,function(f)sum(f$withinerror),numeric(1))
fit <- fits[[which.min(objectives)]]
C <- fit$centers
ord <- order(max.col(C,ties.method='first'),C[,1],C[,2],C[,3],C[,4],method='radix')
C <- C[ord,,drop=FALSE];U <- fit$membership[,ord,drop=FALSE]
cl <- max.col(U,ties.method='first')
stopifnot(all(abs(rowSums(U)-1)<1e-10),all(is.finite(U)),all(tabulate(cl,4)>0))
colnames(M) <- paste0('mean_',g);colnames(Z) <- paste0('Z_',g)
colnames(U) <- paste0('membership_cluster_',1:4)
result <- cbind(d[,c('feature_id','platform','feature_name','source_row')],M,Z,
 cluster=cl,max_membership=apply(U,1,max),U,fuzzifier_m=m,
 seed=20261020,starts=100,objective=min(objectives))
result <- result[order(result$cluster,-result$max_membership,result$feature_name,method='radix'),]
write.csv(result,file.path(out,'FigS3c_results.csv'),row.names=FALSE,fileEncoding='UTF-8',na='NA')

# Native SVG authoring avoids text outlining and Illustrator clipping issues.
xml <- function(s) {
 s <- gsub('&','&amp;',as.character(s),fixed=TRUE)
 s <- gsub('<','&lt;',s,fixed=TRUE);gsub('>','&gt;',s,fixed=TRUE)
}
num <- function(v)formatC(v,digits=4,format='f',decimal.mark='.')
text <- function(x,y,label,size=11,anchor='start',fill='#111111',weight='normal',extra='')sprintf('<text x="%s" y="%s" font-family="Arial" font-size="%s" text-anchor="%s" fill="%s" font-weight="%s" %s>%s</text>',num(x),num(y),size,anchor,fill,weight,extra,xml(label))
line <- function(x1,y1,x2,y2,col='#222222',width=1,dash='')sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="%s" %s/>',num(x1),num(y1),num(x2),num(y2),col,width,if(nzchar(dash))paste0('stroke-dasharray="',dash,'"')else'')
poly <- function(xx,yy,col,width=1,opacity=1,label='')sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="%s" opacity="%s" stroke-linejoin="round"><title>%s</title></polyline>',paste(paste(num(xx),num(yy),sep=','),collapse=' '),col,width,opacity,xml(label))
W <- 780;H <- 530
svg <- c('<?xml version="1.0" encoding="UTF-8"?>',sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%spt" height="%spt" viewBox="0 0 %s %s">',W,H,W,H),
 '<rect width="100%" height="100%" fill="white"/>',text(8,25,'c',23),
 text(W/2,21,'81 candidate features',12,'middle'))
is_lpe <- grepl('^LPE',d$feature_name)
for(k in 1:4) {
 ix <- which(cl==k);col <- (k-1)%%2;row <- (k-1)%/%2
 left <- 67+col*386;right <- left+295;top <- 64+row*228;bottom <- top+164
 xx <- seq(left,right,length.out=4);yy <- function(v)bottom-(v+1.65)/3.3*(bottom-top)
 svg <- c(svg,sprintf('<g id="cluster_%d">',k),text((left+right)/2,top-19,paste0('Cluster ',k,' (n = ',length(ix),')'),14,'middle'))
 for(v in c(-1.5,-1,-.5,0,.5,1,1.5)) {
  svg <- c(svg,line(left,yy(v),right,yy(v),if(v==0)'#BBBBBB'else'#EAEAEA',if(v==0).8 else .5),line(left-4,yy(v),left,yy(v)),text(left-8,yy(v)+3.5,format(v,trim=TRUE),9,'end'))
 }
 for(i in ix[!is_lpe[ix]])svg <- c(svg,poly(xx,yy(Z[i,]),'#B73935',.85,.55,d$feature_name[i]))
 for(i in ix[is_lpe[ix]])svg <- c(svg,poly(xx,yy(Z[i,]),'#32754C',1.65,.95,d$feature_name[i]))
 svg <- c(svg,poly(xx,yy(C[k,]),'#111111',1.7,1,'Fuzzy cluster centroid'),line(left,top,left,bottom),line(left,bottom,right,bottom))
 for(j in 1:4)svg <- c(svg,line(xx[j],bottom,xx[j],bottom+4),text(xx[j],bottom+18,if(j==1)'Ctrl'else g[j],11,'middle'))
 svg <- c(svg,text(left-43,(top+bottom)/2,'Standardized group mean',10,'middle',extra=sprintf('transform="rotate(-90 %s %s)"',num(left-43),num((top+bottom)/2))),'</g>')
}
svg <- c(svg,line(175,506,197,506,'#B73935',1.1),text(203,510,'Other candidates',10),
 line(325,506,347,506,'#32754C',1.7),text(353,510,'LPE features',10),line(460,506,482,506,'#111111',1.7),text(488,510,'Fuzzy centroid',10),'</svg>')
writeLines(svg,file.path(out,'FigS3c_editable.svg'),useBytes=TRUE)
cat('Mfuzz ',as.character(utils::packageVersion('Mfuzz')),'; m=',format(m,digits=12),'; objective=',format(min(objectives),digits=12),'\n',sep='')
print(table(cluster=cl));print(result[grepl('^LPE',result$feature_name),c('feature_name','cluster','max_membership')],row.names=FALSE)
cat('Membership range:',range(result$max_membership),'\n')
