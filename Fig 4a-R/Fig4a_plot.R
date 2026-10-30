# Fig. 4a: discovery-cohort ROC curves for the fixed FOUR-metabolite panel.
# Keep this script and Fig4a_data.csv together; click Source in RStudio.
# Install once: install.packages(c('brglm2','pROC'))
# Tested package versions are printed in the console at completion.
#
# Input: the three current FIG 4/A discovery XLSX files, each containing
# 65 AAD and 29 comparator participants. Duplicated AAD records across the
# three inputs were checked for equality and merged into 152 unique records.
# Current H049/HP-21 values were verified against the latest metabolomics
# master data; no values were imputed, logged, deleted, or substituted.
#
# All four predictors enter each model together: Phe-Try, Carnitine C11:0,
# LPE(18:1) sn-1, LPC(17:0) sn-1. No stepwise selection and no LPE18:0.
# Three separate binary models; outcome AAD=1 and comparator=0.
# Because AAD vs CTRL exhibits complete separation, ALL THREE comparisons
# use Firth mean-bias-reduced logistic regression, brglm2 type='AS_mean'.
# For logistic regression this is equivalent to the Jeffreys-prior penalty
# l(beta)+0.5*log(det(X'WX)); see brglm2 official documentation:
# https://cran.r-project.org/package=brglm2
# Predictors are centered/scaled within each comparison for numerical
# stability only. Coefficients are exported in the ORIGINAL abundance units.
#
# These are APPARENT / in-sample ROC curves on the model-fitting data.
# They are not cross-validated or externally validated estimates. Three
# pairwise fitted models are not a single locked diagnostic classifier.
# Existing five-variable SPSS scores, AUCs and thresholds are not reused.
#
# AUC: empirical, fixed direction (larger score predicts AAD), pROC.
# 95% CI: DeLong, conditional on fitted scores; does not account for model
# fitting or prior variable selection. At AUC=1 the empirical DeLong interval
# degenerates to [1,1]; this is not proof of certain future discrimination.
# Sensitivity/specificity: maximum Youden index on the same discovery data.
# If several thresholds tie, choose highest sensitivity, then specificity,
# then lowest threshold. Rule: predict AAD if probability >= threshold.
#
# Outputs, each with a distinct purpose:
# Fig4a_editable.svg: native text, editable vector curves, original layout.
# Fig4a_summary.csv: AUC/CI/operating point and original-scale coefficients.
# Fig4a_predictions.csv: model-specific probabilities for all 282 observations.
# Fig4a_ROC_coordinates.csv: all thresholds and ROC coordinates.
# Input Fig4a_data.csv remains unchanged. No original XLSX/SAV is modified.

args <- commandArgs(trailingOnly=TRUE)
script_file <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
if(is.null(script_file)) {
 f <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
 if(length(f))script_file <- sub('^--file=','',f[1])
}
out <- if(length(args))args[1] else if(!is.null(script_file))dirname(normalizePath(script_file,winslash='/',mustWork=TRUE)) else getwd()
for(p in c('brglm2','pROC'))if(!requireNamespace(p,quietly=TRUE))stop(paste('Install package',p,'before running this script.'))
d <- read.csv(file.path(out,'Fig4a_data.csv'),fileEncoding='UTF-8-BOM',check.names=FALSE,stringsAsFactors=FALSE)
features <- c('Phe-Try','Carnitine C11:0','LPE(18:1) sn-1','LPC(17:0) sn-1')
stopifnot(identical(names(d),c('sample','sample_id','group',features)),nrow(d)==152L,
 !anyDuplicated(d$sample),!anyDuplicated(d$sample_id),
 identical(as.integer(table(factor(d$group,levels=c('AAD','CTRL','HP','AMI')))),c(65L,29L,29L,29L)))
Xall <- as.matrix(d[,features]);storage.mode(Xall) <- 'double'
stopifnot(all(is.finite(Xall)),all(Xall>=0))
summaries <- list();predictions <- list();curves <- list()
for(g in c('CTRL','HP','AMI')) {
 ix <- which(d$group%in%c('AAD',g));dd <- d[ix,,drop=FALSE]
 X <- Xall[ix,,drop=FALSE];mu <- colMeans(X);sigma <- apply(X,2,sd)
 Z <- sweep(sweep(X,2,mu,'-'),2,sigma,'/')
 modeldata <- data.frame(y=as.integer(dd$group=='AAD'),Z,check.names=FALSE)
 names(modeldata) <- c('y','x1','x2','x3','x4')
 fit <- glm(y~x1+x2+x3+x4,data=modeldata,family=binomial('logit'),
  method=brglm2::brglmFit,type='AS_mean',control=brglm2::brglmControl(maxit=1000,epsilon=1e-10))
 stopifnot(isTRUE(fit$converged),all(is.finite(coef(fit))))
 score <- as.numeric(predict(fit,type='response'))
 beta <- unname(coef(fit)[-1])/sigma
 intercept <- unname(coef(fit)[1])-sum(beta*mu)
 stopifnot(max(abs(plogis(intercept+X%*%beta)-score))<1e-12)
 ro <- pROC::roc(modeldata$y,score,levels=c(0,1),direction='<',quiet=TRUE)
 auc <- as.numeric(pROC::auc(ro))
 # pROC warns at empirical AUC=1; retain the documented degeneracy explicitly.
 ci <- if(auc==1)c(1,1,1)else as.numeric(pROC::ci.auc(ro,method='delong',conf.level=.95))
 cr <- data.frame(comparison=paste('AAD vs',g),threshold=ro$thresholds,
  false_positive_rate=1-ro$specificities,sensitivity=ro$sensitivities,
  specificity=ro$specificities)
 J <- cr$sensitivity+cr$specificity-1;best <- which(abs(J-max(J))<1e-12)
 best <- best[order(-cr$sensitivity[best],-cr$specificity[best],cr$threshold[best])][1]
 op <- cr[best,]
 predictions[[g]] <- data.frame(comparison=paste('AAD vs',g),dd[,1:3],AAD_positive=modeldata$y,
  predicted_probability=score,Youden_threshold=op$threshold,
  predicted_AAD=as.integer(score>=op$threshold),check.names=FALSE)
 summaries[[g]] <- data.frame(comparison=paste('AAD vs',g),method='Firth logistic; four fixed predictors',
  n_AAD=65,n_comparator=29,AUC=auc,CI95_lower=ci[1],CI95_upper=ci[3],CI_method='DeLong; conditional on fitted scores',
  threshold=op$threshold,sensitivity=op$sensitivity,specificity=op$specificity,
  Youden_J=max(J),intercept=intercept,coefficient_PheTry=beta[1],
  coefficient_CarnitineC110=beta[2],coefficient_LPE181sn1=beta[3],coefficient_LPC170sn1=beta[4],
  model_converged=fit$converged,performance='Apparent; discovery model-fitting data',row.names=NULL)
 curves[[g]] <- cr
}
summary <- do.call(rbind,summaries);rownames(summary) <- NULL
write.csv(summary,file.path(out,'Fig4a_summary.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(do.call(rbind,predictions),file.path(out,'Fig4a_predictions.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(do.call(rbind,curves),file.path(out,'Fig4a_ROC_coordinates.csv'),row.names=FALSE,fileEncoding='UTF-8')

# Draw a native SVG with complete text objects, not outlined font glyphs.
xml <- function(s) {
 s <- gsub('&','&amp;',as.character(s),fixed=TRUE);s <- gsub('<','&lt;',s,fixed=TRUE);gsub('>','&gt;',s,fixed=TRUE)
}
num <- function(v)formatC(v,digits=5,format='f',decimal.mark='.')
text <- function(x,y,label,size=11,anchor='start',weight='normal',extra='')sprintf('<text x="%s" y="%s" font-family="Arial" font-size="%s" text-anchor="%s" fill="#111111" font-weight="%s" %s>%s</text>',num(x),num(y),size,anchor,weight,extra,xml(label))
line <- function(x1,y1,x2,y2,col='#111111',width=1.1) sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="%s"/>',num(x1),num(y1),num(x2),num(y2),col,width)
poly <- function(xx,yy) sprintf('<polyline points="%s" fill="none" stroke="#F02222" stroke-width="1.55" stroke-linejoin="miter"/>',paste(paste(num(xx),num(yy),sep=','),collapse=' '))
W <- 790;H <- 288
svg <- c('<?xml version="1.0" encoding="UTF-8"?>',sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%spt" height="%spt" viewBox="0 0 %s %s">',W,H,W,H),
 '<rect width="100%" height="100%" fill="white"/>',text(8,22,'a',20))
for(k in 1:3) {
 g <- c('CTRL','HP','AMI')[k];cr <- curves[[g]];ss <- summaries[[g]]
 left <- 61+(k-1)*257;right <- left+192;top <- 40;bottom <- top+192
 xp <- function(v)left+v*(right-left);yp <- function(v)bottom-v*(bottom-top)
 svg <- c(svg,sprintf('<g id="ROC_AAD_vs_%s">',g),text((left+right)/2,top-12,paste('AAD VS.',g),13,'middle','bold'),
  line(left,top,left,bottom),line(left,bottom,right,bottom),line(right,bottom,right,top),line(right,top,left,top))
 for(v in seq(0,1,.2)) {
  svg <- c(svg,line(xp(v),bottom,xp(v),bottom+5),text(xp(v),bottom+18,sprintf('%.1f',v),10,'middle'),
   line(left-5,yp(v),left,yp(v)),text(left-9,yp(v)+3.5,sprintf('%.1f',v),10,'end'))
 }
 # Order ROC points from (0,0) to (1,1); no interpolation or smoothing.
 ord <- order(cr$false_positive_rate,cr$sensitivity,method='radix')
 svg <- c(svg,poly(xp(cr$false_positive_rate[ord]),yp(cr$sensitivity[ord])),
  text((left+right)/2,bottom+42,'1-Specificity',12,'middle','bold'),
  text(left-39,(top+bottom)/2,'Sensitivity',12,'middle','bold',extra=sprintf('transform="rotate(-90 %s %s)"',num(left-39),num((top+bottom)/2))))
 labels <- c(sprintf('AUC = %.3f',ss$AUC),sprintf('(95%% CI: %.3f-%.3f)',ss$CI95_lower,ss$CI95_upper),
  sprintf('Sensitivity = %.1f%%',100*ss$sensitivity),sprintf('Specificity = %.1f%%',100*ss$specificity))
 for(j in 1:4)svg <- c(svg,text((left+right)/2,top+101+(j-1)*16,labels[j],11,'middle'))
 svg <- c(svg,'</g>')
}
svg <- c(svg,text(W/2,H-4,'Firth logistic regression; four metabolites; apparent discovery performance',8.5,'middle'),'</svg>')
writeLines(svg,file.path(out,'Fig4a_editable.svg'),useBytes=TRUE)
print(summary[,c('comparison','AUC','CI95_lower','CI95_upper','sensitivity','specificity','threshold')],row.names=FALSE)
cat('R',as.character(getRversion()),'; brglm2',as.character(packageVersion('brglm2')),'; pROC',as.character(packageVersion('pROC')),'\n')
