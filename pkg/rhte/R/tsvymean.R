tsvymean <-
function(y, design, trim=c(0, 0), type = c("trim", "win"), na.rm = FALSE, control = rht.control(...), ...){
   # retrieve some parameters
   quietly <- control$quietly
   ddigits <- control$ddigits
   #
   type <- match.arg(type)
   method <- switch(type,
      trim = "trimmed weighted mean",
      win = "winsorized weighted mean")
   check <- .rsvycheck(y, design, na.rm=na.rm)
   if(dim(check$mat)[2] > 2) stop("Only one single 'y'-variable \n")
   method <- paste(method, check$method, sep="")
   x <- check$mat[, 1 : (dim(check$mat)[2] - 1)]
   w <- as.vector(check$mat[, dim(check$mat)[2]])
   n <- length(w)
   w.cumsum <- cumsum(w)
   # trim consists of the amount to be removed above and below 
   chk <- trim %% 1
   if (sum(chk) == 0){
      # trim is integer
      if (sum(trim) >= n) stop("argument 'trim' is not correctly specified!\n")
      # compute the relative amount of trimming
      if (trim[1] > 0){
	 trim[1] <- w.cumsum[trim[1]] / w.cumsum[n]
      }
      if (trim[2] > 0){
	 trim[2] <- 1 - (w.cumsum[n] - w.cumsum[(n - trim[2])]) / w.cumsum[n]
      }else{
	 trim[2] <- 1
      }
   }else{
      # trim is not integer
      if (trim[1] >= (1 - trim[2])) stop("argument 'trim' is not correctly specified!\n")
      if (trim[1] < 0) stop("argument 'trim' is not correctly specified!\n")
      if ((1 - trim[2]) > 1) stop("argument 'trim' is not correctly specified!\n")
      trim[2] <- 1 - trim[2]
   }
   #
   if(type == "trim"){
      tmp <- .trimmedmean(x, w, trim[1], trim[2])
      average <- tmp$tm
      u <- tmp$u
      #
      v <- svyrecvar((.infltrim(x, w, trim[1], trim[2])) * w / sum(w), design$cluster, design$strata, design$fpc, postStrata = design$postStrata)
   }
   if(type == "win"){
      average <- .winsorizedmean(x, w, trim[1], trim[2])
      v <- svyrecvar((.infltrim(x, w, trim[1], trim[2])) * w / sum(w), design$cluster, design$strata, design$fpc, postStrata = design$postStrata)
   v.win <- svyrecvar((.inflwinsor(x, w, trim[1], trim[2])) * w / sum(w), design$cluster, design$strata, design$fpc, postStrata = design$postStrata)
   }
   names(average) <- colnames(check$mat)[1 : (dim(check$mat)[2] - 1)]
   attr(average, "var") <- v
   attr(average, "nobs") <- NROW(x)
   attr(average, "na") <- check$na
   attr(average, "method") <- method
   attr(average, "statistic") <- "mean"
   attr(average, "outliers") <- floor(NROW(x) * trim[1]) + floor((1 - trim[2]) * NROW(x)) 
   attr(average, "design") <- design
   attr(average, "call") <-  match.call()
   if(type == "win"){
      rob <- cbind(trim[1], trim[2], round(sqrt(v.win), ddigits))
      colnames(rob) <- c("lo", "hi", "std.err.winsor")
   }
   else{
      rob <- cbind(trim[1], trim[2])
      colnames(rob) <- c("lo", "hi")
   }
   rob.names <- switch(type,
      trim = "trimming",
      win = "winsorized")
   rownames(rob) <- rob.names	
   attr(average, "rob") <- rob
   attr(average, "ddigits") <- ddigits
   attr(average, "robweights") <- u
   attr(average, "resid") <- u
   class(average) <- "svystat.rob"
   if(!quietly){
      summary(average)
      invisible(average)
   }else{
      invisible(average)
   }
}
