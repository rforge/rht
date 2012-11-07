summary.svystat.rob <-
function(object, ...){
   ddigits <- attr(object, "ddigits")
   if(length(object)!=1) stop("no summary available\n")
   cat("\n")
   cat("SUMMARY FOR:", attr(object, "method"), "\n")
   cat("---\n")
   est <- cbind(round(object, ddigits), round(sqrt(attr(object, "var")), ddigits), attr(object, "outliers"), attr(object, "nobs"), attr(object, "na"))
   colnames(est) <- c(attr(object, "statistic"), "SE", "outliers", "nobs", "NA's")
   print(est) 
   cat("---\n")
   cat("ROBUSTNESS PROPERTIES\n")
   print(attr(object, "rob"), quote = FALSE)
   cat("---\n")
   if(!is.null(attr(object, "optim"))){
   cat("ALGORITHM PERFORMANCE \n")
   print(attr(object, "optim"))
   cat("---\n")
   }
   cat("SAMPLING DESIGN\n")
   print(attr(object, "design"))
}
