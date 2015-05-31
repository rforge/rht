.huberpsiw <-
function(x, k, w=FALSE, asymmetric=FALSE){
   # test whether k is valid
   if (k <= 0) stop("k must not be <= 0!\n")
   # weights
   if (w){
      if (asymmetric){
	 res <- pmin(1, k/x)
      }else{
	 res <- pmin(1, k / abs(x))
      }
   }else{
   # psi
      if (asymmetric){
         res <- pmin(k, x)
      }else{
	 res <- pmin(k, pmax(-k, x))
      }
   }
   return(res)
}
.inflmad <-
function(x, weight){
   est <- .svyquantiles(x, weight, probs=0.5)
   S <- .svyquantiles(abs(x - est), weight, probs=0.5) * 1.4826
   # density estimate at the median
   f <- .svydensity(x, weight, quant=0.5, ngrid=401)
   # quantile of T+S
   at.plus <-  sum( weight[1:length(which(x <= est+S))] ) / sum(weight)
   # quantile of T-S		
   at.minus <-  sum( weight[1:length(which(x <= est-S))] ) / sum(weight)   
   #density estimate at the quantile(T+S)
   f.plus <-  .svydensity(x, weight, quant=at.plus, ngrid=401)
   #density estimate at the quantile(TSS)
   f.minus <- .svydensity(x, weight, quant=at.minus, ngrid=401)
   numerator <- sign(abs(x - est) - S) - ((f.plus - f.minus) / f) * sign(x - est)
   denominator <- 2 * (f.plus + f.minus)
   res <- as.vector(numerator / denominator)
}
.inflmest <-
function(x, weight, est, k, S, asymmetric, exact){
   u <- (x - est) / S
   # psi
   psi.value <- .huberpsiw(u, k, w=FALSE, asymmetric=asymmetric)
   # psi prime
   if (asymmetric){
      psiprime.value <- 1 * (u <= k) 
   }else{
      psiprime.value <- 1 * (abs(u) <= k)
   }
   # exact or not?
   if(exact){
      numerator <- psi.value * S - .inflmad(x, weight) * (sum(weight * u * psiprime.value)) / (sum(weight))
   }
   else{
      numerator <- psi.value
   }
   denominator <- sum(weight * psiprime.value) / sum(weight)
   res <- numerator / denominator
   return(res)
}
.infltrim <-
function(x, weight, lo=0, hi=1){
	# x <- sort(x)
	if(is.null(weight)){
		weight <- rep(1, length(x)) 
	}	
	xilo <- .svyquantiles(x, weight, probs=lo)	
	xihi <- .svyquantiles(x, weight, probs=hi)	
	below <- floor(lo * length(x))
	above <- ceiling(hi * length(x))	
	
	mat <- c(rep((1 - lo)*xilo - (1 - hi)*xihi, below),
			rep(-lo*xilo - (1 - hi)*xihi, (above - below)),
			rep(hi*xihi - lo*xilo, (length(x) - above)))	
	if(below != 0){
		x[1:below] <- 0
	}
	if(above != length(x)){
		x[(above + 1):length(x)] <- 0
	}	
	phi <- (x + mat)*(1/(hi - lo)) - .trimmedmean(x, weight, lo, hi)$tm 
	return(phi)
}
.inflwinsor <-
function(x, weight, lo=0, hi=1){
	# x <- sort(x)				
	if(is.null(weight)){
		weight <- rep(1, length(x)) 
	}	
	xilo <- .svyquantiles(x, weight, probs=lo)	
	xihi <- .svyquantiles(x, weight, probs=hi)
	below <- floor(lo * length(x))
	above <- ceiling(hi * length(x))	
	dens.lo <- .svydensity(x, weight, quant=lo)
	dens.hi <- .svydensity(x, weight, quant=hi)
	winlo <- lo - c(rep(1, below),rep(0, (length(x) - below)))
	winhi <- hi - c(rep(1, above),rep(0, (length(x) - above)))
	mat <- c(rep((1 - lo)*xilo - (1 - hi)*xihi, below),
			rep(-lo*xilo - (1 - hi)*xihi, (above - below)),
			rep(hi*xihi - lo*xilo, (length(x) - above)))
	if(below != 0){
		x[1:below] <- 0
	}
	if(above != length(x)){
		x[(above + 1):length(x)] <- 0
	}	
	phi1 <- (x + mat) - .winsorizedmean(x, weight, lo, hi)*(hi - lo)
	if(dens.lo==0){
		phi <- phi1 + (1-hi)*(winhi/dens.hi)	
	}
	else{
		phi <- phi1 + lo*(winlo/dens.lo) + (1-hi)*(winhi/dens.hi)	
	}
	return(phi)
}
.irls <-
function(k, mat, acc, steps, design, beta.med, mer=FALSE, true.mean=NULL, asymmetric){
   niter  <- 0
   beta.0 <- beta.med + 2 * acc
   beta.r <- beta.med
   x.sqrt <- sqrt(mat[, 3])
   while ( (abs(beta.r - beta.0) > acc ) & ( niter < steps) ){
      beta.0 <- beta.r
      stand.res <- (mat[, 1] - beta.0 * mat[, 3]) / x.sqrt
      scale <- .svyquantiles(abs(stand.res), mat[, 2]) * 1.4829
      if (scale == 0) {
	 scale <- .svyquantiles(abs(stand.res), mat[, 2], probs=0.8) * 0.7803
	 cat("MAD is 0, trying quantile 0.8 of absolute deviations\n")
      }		
      u <- .huberpsiw(stand.res / scale, k, w=TRUE, asymmetric=asymmetric)
      beta.r <- sum(mat[, 2] * u * mat[, 1]) / sum(mat[, 2] * u * mat[, 3])
      niter <- niter + 1
   }
   ures <- u * (mat[, 1] - beta.r * mat[, 3])
   design1 <- update(design, ures)
   se <- as.vector(sum(mat[, 2]) * sqrt(vcov(svymean(~ures, design1))) / sum(mat[, 2] * u * mat[, 3]))
   if(mer){
      mse <- as.vector(se^2 + (beta.r - true.mean)^2)  
      return(mse)	
   }
   else{
      result <- list(m=beta.r, se=se, ures=ures, u=u, niter=niter, scale=scale)
      return(result)
   }	
}
.rsvycheck <-
function(x, design, na.rm=TRUE){
   if(!inherits(design, "survey.design2")){
      stop("Design must be an object of the 'survey.design (2)' class\n")
   }	
   method <- ""
   if(design$call[[1]] == "subset"){
      method <- paste(method, " (for domains)", sep="")
   }
   original.design <- design
   if(!inherits(x, "formula")){
      stop("y must be a formula object\n")
   }
   else{
      # ordering the data acording to the study variable
      # and generate a new svydesign with the ordered values
      ordereddata <- design$variables[order(design$variables[paste(x[2])]), ]
      design$call$data <- substitute(ordereddata)
      design <- eval(design$call)	
      # extract the variables from the svydesign -> xmat
      mf <- model.frame(x, design$variables, na.action = na.pass)
      xx <- lapply(attr(terms(x), "variables")[-1], function(tt) model.matrix(eval(bquote(~0 + .(tt))), mf))
      cols <- sapply(xx, NCOL)
      xmat <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))
      scols <- c(0, cumsum(cols))
      for (i in 1:length(xx)) {
	 xmat[, scols[i] + 1:cols[i]] <- xx[[i]]
      }
      colnames(xmat) <- do.call("c", lapply(xx, colnames))
      xmat <- cbind(xmat, as.matrix(weights(design)))
      colnames(xmat)[NCOL(xmat)] <- "w"
   }
   # check if xmat has NAs
   if(na.rm){
      nas <- rowSums(is.na(xmat))
      design <- design[nas == 0]
      xmat <- xmat[nas == 0, , drop=F]
      na <- length(which(nas!=0))
   }
   else{
      na <- NA
   }
   if(any(xmat[,NCOL(xmat)] == 0)) cat("There are weights=0\n")
   res <- list(mat=xmat, method=method, na=na, design=design)
   class(res) <- "rsvycheck"
   return(res)
}
.status <-
function(at, total){
	if(total >= 30){
		p.indicator <- diff(floor(30 * seq(1, total+1) / total))
		if(at==1){
			cat("  Progress \n")
			cat("  +------+------+------+------+ \n")
			cat("  0%    25%    50%    75%    100% \n")
			cat("  ")
		}
		else{
			if(p.indicator[at]==1){
				cat("*")	
			}
		}
	}
	else{
		repeating <- diff(floor(seq(1, 30, length.out=total+1)))
		if(at==1){
			cat("  Progress \n")
			cat("  +------+------+------+------+ \n")
			cat("  0%    25%    50%    75%    100% \n")
			cat("  ")
			for(j in 1:repeating[at]){
				cat("*")
			}
		}
		else{			
			for(j in 1:repeating[at]){
				cat("*")
			}
		}
	}
	if(at == total){
		cat("\n ...done \n")
	}
}
.svydensity <-
function(x, weight, quant=0.5, ngrid=401){
	bwd <- dpik(x, scalest="minim", level=2, kernel="normal",   
			canonical=FALSE, gridsize=ngrid, range.x=range(x), 
			truncate=TRUE)	
	at <- seq(min(x), max(x), length=ngrid)
	nx <- rowsum(c(rep(0,ngrid),weight), c(1:ngrid, findInterval(x,at)))
	dens <-locpoly(rep(1,ngrid),nx*ngrid/(diff(range(x))*sum(weight)),
			binned=TRUE, bandwidth=bwd, range.x=range(x))
	res <- dens$y[floor(length(dens$y)*quant)]
	if(length(res)==0) res <- 0
	return(res)
}
.svyquantiles <-
function(x, w, ties=NULL, probs = 0.5, method="constant", f=0){
	ord <- order(x)
	x <- x[ord]
	w <- w[ord]
	if(is.null(ties)){
		w.ord <- cumsum(w) / sum(w)
	}
	else{
		w <- rowsum(w, ties, reorder = F)	
		x <- sort(aggregate(x, by=list(ties), unique)$x)
		w.ord <- cumsum(w) / sum(w)
	}
	cdf <- approxfun(w.ord, x, method = method, f = f,			
			yleft = min(x), yright = max(x), rule=2)
	cdf(probs)	
}
.trimmedmean <-
function(x, weight, lo=0, hi=1){
   xilo <- .svyquantiles(x, weight, probs=lo)		
   xihi <- .svyquantiles(x, weight, probs=hi)
   w <- weight[x >= xilo & x <= xihi]
   sw <- sum(w)
   tm <- sum(w * x[x>=xilo & x<= xihi]) / sw 
   u <- rep(1, length(weight))
   u[x < xilo] <- 0
   u[x > xihi] <- 0
   u <- u * (sum(weight) / sw)
   res <- list(tm=tm, u=u)
   return(res)
}
.updte <-
function(y, design, multiply){
   variable.string <- paste(y[[2]])
   myvar <- design$variables[, variable.string]
   if (length(multiply) != length(myvar)) stop("dimension of 'multiply' is not correct!\n")
   design$variables[, variable.string] <- myvar * multiply
   return(design)
}
.winsorizedmean <-
function(x, weight, lo=0, hi=1){
	xilo <- .svyquantiles(x, weight, probs=lo)		
	xihi <- .svyquantiles(x, weight, probs=hi)	
	x[x < xilo] <- xilo
	x[x > xihi] <- xihi 
	res <- sum(weight*x) / sum(weight)
	return(res)
}
