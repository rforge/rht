rht.control <-
function(acc=1e-6, steps=50, asymmetric=FALSE, quietly=FALSE, exact=FALSE, ddigits=3, ...){
   eps <- .Machine$double.eps^(1/4) 
   res <- list(acc=acc, steps=steps, asymmetric=asymmetric, quietly=quietly, exact=exact, eps=eps, ddigits=ddigits)
   return(res)
}
