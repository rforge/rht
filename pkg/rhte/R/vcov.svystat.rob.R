vcov.svystat.rob <-
function(object, ...){
v <- as.matrix(attr(object, "var"))
rownames(v) <- names(object)
colnames(v) <- "Variance"
return(v)
}
