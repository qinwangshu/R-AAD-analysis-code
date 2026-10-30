# Fig. 4b: FOUR-metabolite Firth logistic models REFITTED ON VALIDATION DATA.
# Required reviewer inputs: Fig4b_plot.R + Fig4b_data.csv, in the same folder.
# Click Source in RStudio; install once: install.packages(c('brglm2','pROC')).
# Tested with R4.6.1, brglm2 1.1.1, pROC 1.19.1.
#
# Data: current FIG 4/B/验证集AAD、CTRL.xlsx, 验证集AAD、HP.xlsx,
# 验证集AAD、AMI.xlsx, Sheet1. All eight source variables were checked
# against the master 最终数据处理验证集  2026.10.16.xlsx, with zero mismatches.
# Identical repeated AAD records across the pairwise tables were consolidated.
# There are 313 unique observations: AAD68, CTRL94, HP90 and AMI61.
# Only four fixed predictors are used: Phe-Try, Carnitine C11:0,
# LPE(18:1) sn-1, LPC(17:0) sn-1. No stepwise or VIP selection.
#
# Three separate Firth logistic models are newly fitted within this cohort,
# AAD=1 versus comparator=0, using brglm2 type='AS_mean'. In logistic models,
# Firth mean bias reduction is equivalent to the Jeffreys-prior penalty.
# Official documentation: https://cran.r-project.org/package=brglm2
# Original discovery-model coefficients and cutoffs are NOT applied.
# Numerical robustness: deterministic BFGS restarts initialize additional
# brglm2 fits; retain the finite converged fit with greatest penalized
# likelihood. AUC, sensitivity and specificity are not used to choose a fit.
# Predictors are centered/scaled within each pair for numerical stability;
# exported coefficients use the original abundance units. No log, imputation,
# outlier removal or data deletion. Phe-Try zeros (AAD8, AMI2) are confirmed
# nondetections and remain zero. Some extreme probabilities numerically reach
# 0/1; all coefficients must be finite and convergence is explicitly checked.
#
# AUC is empirical with a fixed direction: higher score predicts AAD.
# 95% CI is DeLong conditional on fitted scores; it does not include model-
# fitting or variable-selection uncertainty. AUC=1 gives the degenerate
# empirical DeLong interval [1,1], not certain perfect future performance.
# Each cutoff is independently selected by maximum Youden index within the
# current cohort; no cutoff is fixed to the discovery results. Ties: highest
# sensitivity, then specificity, then lowest threshold. Classify AAD when
# probability >= threshold. No cutoff markers are drawn on the curves.
#
# Both model fitting and cutoff selection use these same observations.
# Report the curves as APPARENT / within-cohort refitted performance, NOT
# independent external validation of the discovery models. This distinction
# must also be retained in the manuscript. Three pairwise models are not a
# single unified locked classifier. No original XLSX or manuscript is modified.
#
# Automatic outputs: editable native-text SVG, summary with original-scale
# coefficients, predicted probabilities (449 pairwise rows), ROC coordinates.
# Input Fig4b_data.csv remains unchanged. No extra input files are required.

args <- commandArgs(trailingOnly=TRUE)
script_file <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
if(is.null(script_file)) {
 f <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
 if(length(f))script_file <- sub('^--file=','',f[1])
}
out <- if(length(args))args[1] else if(!is.null(script_file))dirname(normalizePath(script_file,winslash='/',mustWork=TRUE)) else getwd()
for(p in c('brglm2','pROC'))if(!requireNamespace(p,quietly=TRUE))stop(paste('Install package',p,'before running this script.'))
d <- read.csv(file.path(out,'Fig4b_data.csv'),fileEncoding='UTF-8-BOM',check.names=FALSE,stringsAsFactors=FALSE)
features <- c('Phe-Try','Carnitine C11:0','LPE(18:1) sn-1','LPC(17:0) sn-1')
stopifnot(identical(names(d),c('sample','group',features)),nrow(d)==313L,
 !anyDuplicated(d$sample),
 identical(as.integer(table(factor(d$group,levels=c('AAD','CTRL','HP','AMI')))),c(68L,94L,90L,61L)))
Xall <- as.matrix(d[,features]);storage.mode(Xall) <- 'double'
stopifnot(all(is.finite(Xall)),all(Xall>=0))
fit_firth <- function(md) {
 Q <- model.matrix(y~x1+x2+x3+x4,md);y <- md$y
 objective <- function(b,gradient=FALSE) {
  eta <- as.numeric(Q%*%b);p <- plogis(eta);w <- p*(1-p);I <- crossprod(Q,w*Q)
  detI <- determinant(I,logarithm=TRUE)
  if(detI$sign<=0)return(if(gradient)rep(0,ncol(Q))else 1e100)
  if(gradient) {
   h <- w*rowSums(Q*t(solve(I,t(Q))))
   return(-as.numeric(crossprod(Q,y-p+h*(.5-p))))
  }
  # Stable log(1+exp(eta)).
  -sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))))-.5*as.numeric(detI$modulus)
 }
 fit <- function(start=NULL)glm(y~x1+x2+x3+x4,data=md,family=binomial('logit'),
  method=brglm2::brglmFit,type='AS_mean',start=start,
  control=brglm2::brglmControl(maxit=1000,epsilon=1e-10))
 f1 <- fit();fits <- list(f1)
 starts <- list(rep(0,ncol(Q)),unname(coef(f1))*.25,unname(coef(f1))*.5,
  unname(coef(f1))*2,c(0,-1,-1,-1,-1),c(1,-1,1,-1,1))
 RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(20261020)
 starts <- c(starts,lapply(seq_len(30),function(i)rnorm(ncol(Q),sd=if(i<=15)1 else 3)))
 for(start in starts) {
  ff <- tryCatch({
   opt <- optim(start,objective,function(b)objective(b,TRUE),method='BFGS',control=list(maxit=10000,reltol=1e-13))
   fit(unname(opt$par))
  },error=function(e)NULL)
  if(!is.null(ff))fits[[length(fits)+1L]] <- ff
 }
 ok <- vapply(fits,function(f)isTRUE(f$converged)&&all(is.finite(coef(f))),logical(1))
 stopifnot(any(ok));fits <- fits[ok]
 loss <- vapply(fits,function(f)objective(coef(f)),numeric(1))
 best <- fits[[which.min(loss)]]
 stopifnot(max(abs(objective(coef(best),TRUE)))<1e-6)
 list(fit=best,penalized_loglik=-min(loss))
}

summaries <- list();predictions <- list();curves <- list()
for(g in c('CTRL','HP','AMI')) {
 ix <- which(d$group%in%c('AAD',g));dd <- d[ix,,drop=FALSE]
 X <- Xall[ix,,drop=FALSE];mu <- colMeans(X);sigma <- apply(X,2,sd)
 Z <- sweep(sweep(X,2,mu,'-'),2,sigma,'/')
 modeldata <- data.frame(y=as.integer(dd$group=='AAD'),Z,check.names=FALSE)
 names(modeldata) <- c('y','x1','x2','x3','x4')
 ff <- fit_firth(modeldata);fit <- ff$fit
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
 predictions[[g]] <- data.frame(comparison=paste('AAD vs',g),dd[,1:2],AAD_positive=modeldata$y,
  predicted_probability=score,Youden_threshold=op$threshold,
  predicted_AAD=as.integer(score>=op$threshold),check.names=FALSE)
 summaries[[g]] <- data.frame(comparison=paste('AAD vs',g),method='Firth logistic; four fixed predictors',
  n_AAD=sum(modeldata$y==1),n_comparator=sum(modeldata$y==0),AUC=auc,CI95_lower=ci[1],CI95_upper=ci[3],CI_method='DeLong; conditional on fitted scores',
  threshold=op$threshold,sensitivity=op$sensitivity,specificity=op$specificity,
  Youden_J=max(J),intercept=intercept,coefficient_PheTry=beta[1],
  coefficient_CarnitineC110=beta[2],coefficient_LPE181sn1=beta[3],coefficient_LPC170sn1=beta[4],
  model_converged=fit$converged,penalized_loglik=ff$penalized_loglik,performance='Apparent; refitted on validation cohort; cutoff selected within cohort',row.names=NULL)
 curves[[g]] <- cr
}
summary <- do.call(rbind,summaries);rownames(summary) <- NULL
write.csv(summary,file.path(out,'Fig4b_summary.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(do.call(rbind,predictions),file.path(out,'Fig4b_predictions.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(do.call(rbind,curves),file.path(out,'Fig4b_ROC_coordinates.csv'),row.names=FALSE,fileEncoding='UTF-8')

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
 '<rect width="100%" height="100%" fill="white"/>',text(8,22,'b',20))
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
 for(j in seq_along(labels))svg <- c(svg,text((left+right)/2,top+101+(j-1)*16,labels[j],11,'middle'))
 svg <- c(svg,'</g>')
}
svg <- c(svg,'</svg>')
writeLines(svg,file.path(out,'Fig4b_editable.svg'),useBytes=TRUE)
print(summary[,c('comparison','AUC','CI95_lower','CI95_upper','sensitivity','specificity','threshold')],row.names=FALSE)
cat('R',as.character(getRversion()),'; brglm2',as.character(packageVersion('brglm2')),'; pROC',as.character(packageVersion('pROC')),'\n')
