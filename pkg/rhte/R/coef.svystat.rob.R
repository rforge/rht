coef.svystat.rob <-
function(object, ...){
attr(object, "var") <- NULL
attr(object, "nobs") <- NULL
attr(object, "na") <- NULL
attr(object, "method") <- NULL
attr(object, "statistic") <- NULL
attr(object, "outliers") <- NULL
attr(object, "optim") <-  NULL
attr(object, "rob") <-  NULL
attr(object, "design") <- NULL
attr(object, "call") <- NULL
attr(object, "ddigits") <- NULL
attr(object, "type") <- NULL
unclass(object)
}
