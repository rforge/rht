msvyratio <-
function(numerator, denominator, design, k, na.rm=FALSE, control=rht.control(...), ...){
   # retrieve the control parameters
   acc <- control$acc
   steps <- control$steps
   quietly <- control$quietly
   asymmetric <- control$asymmetric
   exact <- control$exact
   ddigits <- control$ddigits
   # test the design object
   if(!is(design, "survey.design2")){
   stop("Design must be an object of the 'survey.design' class\n")
   }
   original.design <- design
   method <- "robust ratio estimator"
   # add the one-step estimator as an original estimator
   add.method <- ifelse(steps==1, " (One-step M-estimation)", " (M-estimation)")
   method <- paste(method, add.method, sep="")
   if(design$call[[1]] == "subset"){
      method <- paste(method, " (for Domains)", sep="")
   }
   # check whether the numerator formula object is valid
   if(inherits(numerator, "formula")){
      num <- model.frame(numerator, design$variables, na.action=na.pass)
      if(NCOL(num) > 1){
	 stop("Only one 'y'-variable allowed!\n")
      }
      mat <- matrix(nrow=NROW(num), ncol=3)
      mat[,1] <- as.matrix(num)
      mat[,2] <- as.matrix(weights(design))
   }
   else{
      stop("Numerator must be a formula object!\n")
   }
   # check whether the denominator formula object is valid
   if(inherits(denominator, "formula")){
      den <- model.frame(denominator, design$variables, na.action=na.pass)
      if(NCOL(den) > 1){
	 stop("Only one 'y'-variable allowed!\n")
      }
      mat[,3] <- as.matrix(den)
      colnames(mat) <- c("numerator", "w", "denominator")
   }
   else{
      stop("Denominator must be a formula object!\n")
   }
   # na checks (problem: rsvycheck tests only for one variable...)
   # FIXME
   if(na.rm){
      nas <- rowSums(is.na(mat))
      design <- design[nas == 0]
      mat <- mat[nas == 0, , drop=F]
      na <- length(which(nas!=0))
   }
   if(any(mat[,2] == 0) & !quietly) cat("Some weights are zero!\n")
   if(any(mat[,3] == 0) & !quietly) cat("There are denominators=0\n")
   # compute the quantiles
   x.quant <- .svyquantiles(x=mat[,3], w=mat[,2], probs=c(0.5, 0.8))
   y.quant <- .svyquantiles(x=mat[,1], w=mat[,2], probs=c(0.5, 0.8))
   # set the starting values
   if (x.quant[1] == 0){
      beta.med <- y.quant[1]
      initial <- round(y.quant[1], 3)
   }
   else{
      beta.med <- y.quant[1] / x.quant[1]
      initial <- round(beta.med, 3)
   }
   if (beta.med == 0){
      beta.med <- y.quant[2] / x.quant[2]
   }
   result <- .irls(k, mat, acc, steps, design, beta.med, mer=FALSE, true.mean=NULL, asymmetric=asymmetric)  
   #
   ratio <- result$m
   names(ratio) <- paste(names(num), "/", names(den), sep="")
   # FIXME: variance if wrong: must be linearized variance of the ratio estimator
   attr(ratio, "var") <- result$se^2
   attr(ratio, "nobs") <- NROW(mat)
   attr(ratio, "na") <- na
   attr(ratio, "method") <- method
   attr(ratio, "statistic") <- "ratio"
   attr(ratio, "outliers") <- sum(result$u<1)
   attr(ratio, "robweights") <- result$u
   attr(ratio, "resid") <- as.vector(result$ures)
   optim <- cbind(result$niter, round(initial, ddigits), acc,  round(result$scale, ddigits))
   colnames(optim) <- c("iterations", "initial", "precision", "scale")
   rownames(optim) <- ""
   attr(ratio, "optim") <- optim
   psi.function <- ifelse(asymmetric, "asym. Huber", "Huber")
   rob <- cbind(psi.function, round(mean(result$u), ddigits), k)
   colnames(rob) <- c("psi", "ave.weight", "k")
   rownames(rob) <- ""
   attr(ratio, "rob") <- rob
   attr(ratio, "design") <- original.design
   attr(ratio, "call") <-  match.call()
   attr(ratio, "ddigits") <- ddigits
   attr(ratio, "asymmetric") <- asymmetric
   class(ratio) <- c("svystat.rob", "plotable")
   if(!quietly){
      summary(ratio)
      return(ratio)
   }else{
      return(ratio)
   }
}
