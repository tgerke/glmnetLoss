# t0 is defined as in 0 < FPR < t0 (Pepe, 2000)
# this estimator follows from Wang2011
# assumes continuity in x (does not account for ties)
pauc <- function(y, x, t0=.2) {
   cases <- which(y==1)
   controls <- which(y==0)
   controlq <- quantile(x[controls], probs=1-t0, type=7)
   
   sum(sapply(x[cases], function(z) {
      intsum <- sum(x[controls]>z)/length(controls)
      if(intsum < t0) {return(t0-intsum)}
      else {
         return(0)
      }
   })) / length(cases)
}

# this will replace the cv.lognet function in glmnet
cvfunc <- function (outlist, lambda, x, y, weights, offset, foldid,
             type.measure, grouped, keep = FALSE)
{
   prob_min = 1e-05
   prob_max = 1 - prob_min
   nc = dim(y)
   if (is.null(nc)) {
      y = as.factor(y)
      ntab = table(y)
      nc = as.integer(length(ntab))
      y = diag(nc)[as.numeric(y),]
   }
   N = nrow(y)
   nfolds = max(foldid)
   if (!is.null(offset)) {
      is.offset = TRUE
      offset = drop(offset)
   }
   else
      is.offset = FALSE   
   mlami = max(sapply(outlist, function(obj)
      min(obj$lambda)))
   which_lam = lambda >= mlami
   predmat = matrix(NA, nrow(y), length(lambda))
   nlams = double(nfolds)
   for (i in seq(nfolds)) {
      which = foldid == i
      fitobj = outlist[[i]]
      if (is.offset)
         off_sub = offset[which]
      preds = predict(fitobj,
                      x[which, , drop = FALSE],
                      s = lambda[which_lam],
                      offset = off_sub,
                      type = "response")
      nlami = sum(which_lam)
      predmat[which, seq(nlami)] = preds
      nlams[i] = nlami
   }
   
   cvraw = matrix(NA, nfolds, length(lambda))
   good = matrix(0, nfolds, length(lambda))
   for (i in seq(nfolds)) {
      good[i, seq(nlams[i])] = 1
      which = foldid == i
      for (j in seq(nlams[i])) {
         cvraw[i, j] = 	pauc(y[which, 1], predmat[which, j]) #takes the default t0
         #cvraw[i, j] = auc.mat(y[which, ], predmat[which, j], weights[which])
      }
   }
   N = apply(good, 2, sum)
   weights = tapply(weights, foldid, sum)
   
   cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
   cvsd = sqrt(apply(
      scale(cvraw, cvm, FALSE) ^ 2,
      2,
      weighted.mean,
      w = weights,
      na.rm = TRUE
   ) / (N - 1))
   out = list(cvm = cvm, cvsd = cvsd, name = "pAUC")
   if (keep) {out$fit.preval = predmat}
   out
}

# this replaces the cv.glmnet function
rose.glmnet <- function (x, y, weights, offset = NULL, lambda = NULL, 
                         type.measure = c("mse", "deviance", "class", "auc", "mae"), nfolds = 10, foldid, 
                         grouped = TRUE, keep = FALSE, parallel = FALSE, 
                         rose.formula=NULL, rose.N, rose.p, ...) 
{
   if (missing(type.measure)) 
      type.measure = "default"
   else type.measure = match.arg(type.measure)
   if (!is.null(lambda) && length(lambda) < 2) 
      stop("Need more than one value of lambda for cv.glmnet")
   N = nrow(x)
   if (missing(weights)) 
      weights = rep(1, N)
   else weights = as.double(weights)
   y = drop(y)
   glmnet.call = match.call(expand.dots = TRUE)
   which = match(c("type.measure", "nfolds", "foldid", "grouped", 
                   "keep"), names(glmnet.call), F)
   if (any(which)) 
      glmnet.call = glmnet.call[-which]
   glmnet.call[[1]] = as.name("glmnet")
   glmnet.object = glmnet(x, y, weights = weights, offset = offset, 
                          lambda = lambda, ...)
   glmnet.object$call = glmnet.call
   is.offset = glmnet.object$offset
   if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
      nz = predict(glmnet.object, type = "nonzero")
      nz = sapply(nz, function(x) sapply(x, length))
      nz = ceiling(apply(nz, 1, median))
   }
   else nz = sapply(predict(glmnet.object, type = "nonzero"), 
                    length)
   if (missing(foldid)) 
      foldid = sample(rep(seq(nfolds), length = N))
   else nfolds = max(foldid)
   if (nfolds < 3) 
      stop("nfolds must be bigger than 3; nfolds=10 recommended")
   outlist = as.list(seq(nfolds))
   # if (parallel) {
   #    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% 
   #    {
   #       which = foldid == i
   #       if (is.matrix(y)) 
   #          y_sub = y[!which, ]
   #       else y_sub = y[!which]
   #       if (is.offset) 
   #          offset_sub = as.matrix(offset)[!which, ]
   #       else offset_sub = NULL
   #       glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda, 
   #              offset = offset_sub, weights = weights[!which], 
   #              ...)
   #    }
   # }
   # else {
   for (i in seq(nfolds)) {
      which = foldid == i
      if (is.matrix(y)) 
         y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset) 
         offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      outlist[[i]] = glmnet(x[!which, , drop = FALSE], 
                            y_sub, lambda = lambda)
      #offset = offset_sub, weights = weights[!which], ...)
   }
   # }
   fun = paste("cv", class(glmnet.object)[[1]], sep = ".")
   lambda = glmnet.object$lambda
   cvstuff = cvfunc(outlist, lambda, x, y, weights, 
                               offset, foldid, type.measure, grouped, keep)
   cvm = cvstuff$cvm
   cvsd = cvstuff$cvsd
   nas = is.na(cvsd)
   if (any(nas)) {
      lambda = lambda[!nas]
      cvm = cvm[!nas]
      cvsd = cvsd[!nas]
      nz = nz[!nas]
   }
   cvname = cvstuff$name
   out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + 
                 cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
   if (keep) 
      out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
   #lamin = if (cvname == "AUC") 
   lamin = getmin(lambda, -cvm, cvsd)
   #else getmin(lambda, cvm, cvsd)
   obj = c(out, as.list(lamin))
   class(obj) = "cv.glmnet"
   obj
}

unlockBinding("cv.glmnet", getNamespace("glmnet"))
assignInNamespace("cv.glmnet", rose.glmnet, ns="glmnet", envir= getNamespace("glmnet"))
assign("cv.glmnet", rose.glmnet, envir= getNamespace("glmnet"))
lockBinding("cv.glmnet", getNamespace("glmnet"))
#cv.glmnet <- glmnet:::cv.glmnet

# # test
# set.seed(1010)
# n=1000;p=1000
# nzc=trunc(p/10)
# x=matrix(rnorm(n*p),n,p)
# beta=rnorm(nzc)
# fx= x[,seq(nzc)] %*% beta
# eps=rnorm(n)*5
# y=drop(fx+eps)
# px=exp(fx)
# px=px/(1+px)
# y <- as.numeric(y>.5)
# #t0 <- 0.2
# set.seed(1011)
# cvob1=glmnet:::cv.glmnet(x,y,family="binomial", type.measure="auc")
# plot(cvob1)
# 
# library(pROC)
# # check the first 10 pAUCs
# for (i in 1:10) {
#    print(auc(roc(response=y, predictor=x[,i]),
#              partial.auc=c(1,1-.2)))
#    print(pauc(y, x[,i], .2))
# }