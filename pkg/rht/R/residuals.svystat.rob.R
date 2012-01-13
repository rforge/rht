#=============================================================================== 
# SUBJECT:  resid method for class svystat.rob 
# AUTHORS:  October 27, 2011
# LICENSE:  GPL > 2
#-------------------------------------------------------------------------------
#  INPUT
#  object   robustly estimated survey mean
#  ellipse  [not used here; but ellipse is part of the generic method]
#-------------------------------------------------------------------------------
residuals.svystat.rob <- function(object, ...){
   attr(object, "resid")
}

