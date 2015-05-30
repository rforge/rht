robweights <-
function(object){
   if(!inherits(object, "svystat.rob")){
      stop("robweights is not a valid method for this object!\n")
   }
   as.vector(attr(object, "robweights"))
}
