msvymean <-
function(y, design, k, type="rht", na.rm=FALSE, control=rht.control(...), ...){
   # retrieve all parameters
   acc <- control$acc
   steps <- control$steps
   quietly <- control$quietly
   asymmetric <- control$asymmetric
   exact <- control$exact
   ddigits <- control$ddigits
   # remember the call
   call <- match.call()	
   # prepare the method
   method <- switch(type,
      rht = "robust Horvitz-Thompson estimator",
      rwm = "robust weighted mean")
   if (is.null(method)) stop("Type must be either 'rht' or 'rwm' !\n")   
   # add the one-step estimator as an explicit estimator
   add.method <- ifelse(steps == 1, " (One-step M-estimation)", " (M-estimation)")
   method <- paste(method, add.method, sep="")
   # work on the sampling design
   original.design <- design
   # check for survey.design and na
   check <- .rsvycheck(y, design, na.rm=na.rm)
   design <- check$design
   # allow only for one study variable
   if (dim(check$mat)[2] > 2) stop("only one 'y'-variable allowed")
   x <- switch(type,
      rht = mean(check$mat[,2]) / check$mat[,2],
      rwm = rep(1, NROW(check$mat)))
   # mat: matrix[, c(y, weights, x)]
   mat <- cbind(check$mat, as.matrix(x))
   nas <- check$na
   method <- paste(method, check$method, sep="")
   # compute the quantiles
   x.quant <- .svyquantiles(x=mat[,3], w=mat[,2], probs=c(0.5, 0.8))
   y.quant <- .svyquantiles(x=mat[,1], w=mat[,2], probs=c(0.5, 0.8))
   mad.y <- .svyquantiles(abs(mat[,1] - y.quant[1]), mat[,2], probs=0.5)*1.4826
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
   # call irls
   result <- .irls(k, mat, acc, steps, design, beta.med, mer=FALSE, true.mean=NULL, asymmetric=asymmetric)
   average <- result$m
   names(average) <- colnames(mat)[1]
   # IF exact=TRUE: the variance estimates considers the MAD as initial
   # scale estimate, ELSE: variance estimates comes from irls
   if(exact){
      linearized.values <- .inflmest(x=mat[,1], weight=mat[,2], est=result$m, k=k, S=mad.y, asymmetric=asymmetric, exact=exact)
      d <- update(design, zz=linearized.values)
      attr(average, "var") <- vcov(svymean(~zz, d)) 			
   }
   else{
      attr(average, "var") <- result$se^2
   }
   attr(average, "nobs") <- NROW(mat)
   attr(average, "na") <- nas
   attr(average, "method") <- method
   attr(average, "statistic") <- "mean"
   attr(average, "outliers") <- sum(result$u<1)
   attr(average, "robweights") <- result$u
   attr(average, "resid") <- as.vector(result$ures)
   attr(average, "design") <- original.design
   attr(average, "call") <-  call
   attr(average, "type") <- type
   psi.function <- ifelse(asymmetric, "asym. Huber", "Huber")
   rob <- cbind(psi.function, round(mean(result$u), ddigits), k)
   colnames(rob) <- c("psi function   ", "ave.weight", "k")
   rownames(rob) <- ""
   attr(average, "rob") <- rob
   optim <- cbind(result$niter, round(initial, ddigits), acc, round(result$scale, ddigits))
   colnames(optim) <- c("iterations", "initial", "precision", "scale")
   rownames(optim) <- ""
   attr(average, "optim") <- optim
   attr(average, "ddigits") <- ddigits
   attr(average, "asymmetric") <- asymmetric
   class(average) <- c("svystat.rob", "plotable")
   if(!quietly){
      summary(average)	
      return(average)
   }else{
      return(average)
   }
}
