print.svystat.rob <-
function(x, ...){
   ddigits <- attr(x, "ddigits")
   m <- cbind(x, sqrt(attr(x, "var")))
   colnames(m) <- c(attr(x, "statistic"), "SE")
   print(round(m, ddigits))
}
