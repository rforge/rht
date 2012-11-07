mer <-
function(object, init = 0.1, box.lo = 0.0001, tol = 1e-4){
   if(attr(object, "call")[[1]] == "msvyratio.survey.design"){ 
      stop("MER for ratio estimation not yet implemented")
   }
   if(attr(object, "call")[[1]] == "tsvymean.survey.design"){ 
      stop("MER for L-estimation not yet implemented")
   }
   y.variable <- attr(object, "call")$y 
   design <- attr(object, "design")
   acc <- attr(object, "optim")[,3]
   steps <- ifelse(attr(object, "optim")[, 1] == 1, 1, 50)
   type <- attr(object, "type")
   method <- attr(object, "method")
   ddigits <- attr(object, "ddigits")
   asymmetric <- attr(object, "asymmetric")
   #
   check <- .rsvycheck(eval(y.variable), design, na.rm=TRUE)
   x <- switch(type,
      rht = mean(check$mat[,2]) / check$mat[,2],
      rwm  = rep(1, NROW(check$mat)))
   mat <- cbind(check$mat, as.matrix(x))
   nas <- check$na
   x.quant <- .svyquantiles(x = mat[, 3], w = mat[, 2], probs = c(0.5, 0.8))
   y.quant <- .svyquantiles(x = mat[, 1], w = mat[, 2], probs = c(0.5, 0.8))
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
   true.mean <- coef(svymean(eval(y.variable), design, na.rm=TRUE, quitely=TRUE))
   opt <- optim(par = init, .irls, mat = mat, acc = acc, steps = steps, design = design, beta.med = beta.med, mer = TRUE, true.mean = true.mean, asymmetric = asymmetric, method = "L-BFGS-B", lower = box.lo, control = list(REPORT=2, trace=1))
   if(opt$convergence == 0 | opt$convergence == 51){
      mer.estimate <- .irls(k = opt$par, mat = mat, acc = acc, steps = steps, design = design, beta.med = beta.med, mer = FALSE, true.mean = NULL, asymmetric = asymmetric)
      average <- mer.estimate$m
      names(average) <- colnames(mat)[1]
      attr(average, "var") <- mer.estimate$se^2
      attr(average, "nobs") <- NROW(mat)
      attr(average, "na") <- nas
      attr(average, "method") <- paste("Minimum estimated risk estimation, based on \n", method, sep="")
      attr(average, "statistic") <- "mean"
      attr(average, "outliers") <- sum(mer.estimate$u < 1)
      attr(average, "robweights") <- mer.estimate$u
      attr(average, "resid") <- mer.estimate$ures
      attr(average, "design") <- design
      attr(average, "call") <-  match.call()
      attr(average, "type") <- type  
      rob <- cbind(attr(object, "rob")[1], round(mean(mer.estimate$u),5), round(opt$par,5))
      colnames(rob) <- c("psi", "ave.weight", "optimal k")
      rownames(rob) <- ""
      attr(average, "rob") <- rob
      optim <- cbind(opt$counts[1], round(opt$value, 1), acc, init, box.lo)
      colnames(optim) <- c("MER-iter", "mse", "precision", "initial", "low box.constr.")
      rownames(optim) <- "L-BFGS-B"
      attr(average, "optim") <- optim
      attr(average, "ddigits") <- ddigits
      attr(average, "asymmetric") <- asymmetric
      class(average) <- c("svystat.rob", "mer")
      summary(average)
      return(average)
   }
   else{
      cat("no output generated because not converged! \n")
      if(opt$convergence == 1) cat("maximum iterations reached, reduce tolerance < 1e-4 \n")
      if(opt$convergence > 5) cat("L-BFGS-B warning: ", opt$message, "\n")
      cat("see help for modifications of the inital estimate (init), \n the lower box constraint (box.lo), and the numerical convergence \n tolerance (tol)\n ")
   }
}
