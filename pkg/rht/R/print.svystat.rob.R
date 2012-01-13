#=============================================================================== 
# SUBJECT:  Influence Function of the MAD
# AUTHORS:  Tobias Schoch, February 8 2009; mod: October 26, 2011
# LICENSE:  GPL > 2
# COMMENT:  Huber (1981), p.137-38
#-------------------------------------------------------------------------------
#  INPUT
#  x	    robustly estimated survey mean
#  ellipse  [not used here; but ellipse is part of the generic method]
#-------------------------------------------------------------------------------
print.svystat.rob <- function(x, ...){
   ddigits <- attr(x, "ddigits")
   m <- cbind(x, sqrt(attr(x, "var")))
   colnames(m) <- c(attr(x, "statistic"), "SE")
   print(round(m, ddigits))
}

