#=============================================================================== 
# SUBJECT:  resid method for class svystat.rob 
# AUTHORS:  October 27, 2011
# LICENSE:  GPL > 2
#-------------------------------------------------------------------------------
#  INPUT
#  object   robustly estimated survey mean
#  ellipse  [not used here; but ellipse is part of the generic method]
#-------------------------------------------------------------------------------
robweights <- function(object){
   if(!inherits(object, "svystat.rob")){
      stop("robweights is not a valid method for this object!\n")
   }
   as.vector(attr(object, "robweights"))
}

