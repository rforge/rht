msvytotal <-
function(y, design, k, type="rht", na.rm=FALSE, control=rht.control(...), ...){
   # retrieve all parameters
   acc <- control$acc
   steps <- control$steps
   asymmetric <- control$asymmetric
   ddigits <- control$ddigits
   quietly<- control$quietly
   # methods; note that the type here is different from the type argument in msvymean
   method <- switch(type,
      rht = "robust Horvitz-Thompson total",
      rwt = "robust weighted total")
   if (is.null(method)) stop("Type must be either 'rht' or 'rwt' !\n")   
   # add the one-step estimator as an explicit estimator
   add.method <- ifelse(steps == 1, " (One-step M-estimation)", " (M-estimation)")
   method <- paste(method, add.method, sep="")
   # remember the call
   call <- match.call()
   # call msvymean
   type <- ifelse(type == "rwt", "rwm", "rht")
   average <- msvymean(y, design, type=type, k, na.rm, acc=acc, steps=steps, quietly=TRUE, asymmetric=asymmetric, exact=FALSE, ddigits=ddigits)
   # retrieve the robustness weights from estimating the M-estimator
   u <- robweights(average)
   # compute y_i * u_i and update the design
   altered.design <- .updte(y, design, u)
   # compute the robust total
   tmp <-  svytotal(y, altered.design)
   total <- coef(tmp)
   # add the attributes from average
   attributes(total) <- attributes(average)   
   # and modifiy the relevant attributes
   attr(total, "var") <- vcov(tmp)
   attr(total, "call") <- call
   attr(total, "statistic") <- "total"
   attr(total, "method") <- method
   #
   if(!quietly){
      summary(total)	     
      invisible(total)
   }else{
      invisible(total)
   }
}
