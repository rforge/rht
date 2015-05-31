pkgname <- "rhte"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('rhte')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("mer")
### * mer

flush(stderr()); flush(stdout())

### Name: mer
### Title: Minimum Estimated Risk M-Estimation
### Aliases: mer

### ** Examples

## load the data
data(api)
## define "survey.design" for stratified sampling
dstrat <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, 
fpc=~fpc)
## compute the a robust Horvitz-Thompson mean
m1 <- msvymean(~api00, dstrat, type="rht", k=1.3)
## compute the minimum estimated risk (MER) estimator based on m1
m1.mer <- mer(m1)
summary(m1.mer)



cleanEx()
nameEx("msvymean")
### * msvymean

flush(stderr()); flush(stdout())

### Name: msvymean
### Title: Robust M-estimation of the mean for complex samples
### Aliases: msvymean

### ** Examples

## load "api" data set from "survey" package (a description of the data
## set can be found there)
data(api)
## define "survey.design" for stratified sampling
dstrat <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, 
fpc=~fpc)
## compute a robust Horvitz-Thompson estimate for the mean of the 
## variable "api00" (Academic Performance Index in 2000)
rht1 <- msvymean(~api00, dstrat, type="rht", k=1.2)
# get a summary of the estimation
summary(rht1)
## robust Horvitz-Thompson estimates for a domain of the variable. Here
## we are interessted in the robust mean for api00 for 
## (sch.wide == "Yes"). That is the average of the academic performance 
## in 2000 only for the schools that met the school-wide growth target.
msvymean(~api00, subset(dstrat, sch.wide == "Yes"), type="rht", k=1.2)
## to extract the estimate from the object 
coef(rht1) 
## to extract the variance from the object
vcov(rht1)



cleanEx()
nameEx("msvyratio")
### * msvyratio

flush(stderr()); flush(stdout())

### Name: msvyratio
### Title: Robust ratio M-estimation for complex samples
### Aliases: msvyratio

### ** Examples

## load "api" data set from "survey" package (a description of the data
## set can be found there)
data(api)
## define "survey.design" for stratified sampling
dstrat <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, 
fpc=~fpc)
## compute a robust Horvitz-Thompson estimate for the mean of the 
## variable api00 (Academic Performance Index in 2000)
ratio1 <- msvyratio(~api00, ~api99, dstrat, k=1.2, na.rm=TRUE)
## get a summary of the estimation
summary(ratio1)



cleanEx()
nameEx("msvytotal")
### * msvytotal

flush(stderr()); flush(stdout())

### Name: msvytotal
### Title: Robust M-estimation of the total for complex samples
### Aliases: msvytotal

### ** Examples

## load "api" data set from "survey" package (a description of the data
## set can be found there)
data(api)
## define "survey.design" for stratified sampling
dstrat <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, 
fpc=~fpc)
## compute a robust Horvitz-Thompson estimate for the total of the 
## variable "api00" (Academic Performance Index in 2000)
rht1 <- msvytotal(~api00, dstrat, k=1.2)
# get a summary of the estimation
summary(rht1)
## robust Horvitz-Thompson estimates for a domain of the variable. Here
## we are interessted in the robust total for api00 for 
## (sch.wide == "Yes"). That is the average of the academic performance 
## in 2000 only for the schools that met the school-wide growth target.
msvytotal(~api00, subset(dstrat, sch.wide == "Yes"), k=1.2)
## to extract the estimate from the object 
coef(rht1) 
## to extract the variance from the object
vcov(rht1)



cleanEx()
nameEx("rhteutils")
### * rhteutils

flush(stderr()); flush(stdout())

### Name: rhte-utils
### Title: rhte utility functions
### Aliases: summary.svystat.rob print.svystat.rob coef.svystat.rob
###   vcov.svystat.rob residuals.svystat.rob robweights

### ** Examples

## load "api" data set from "survey" package (a description of the data
## set can be found there)
data(api)
## define "survey.design" for stratified sampling
dstrat <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, 
fpc=~fpc)
## compute a robust Horvitz-Thompson estimate for the mean of the 
## variable api00 (Academic Performance Index in 2000)
rht1 <- msvymean(~api00, dstrat, type="rht", k=4)
# get a summary of the estimation
summary(rht1)
## robust Horvitz-Thompson estimates for a domain of the variable. Here
## we are interessted in the robust mean for api00 in case of 
## (sch.wide == "Yes"). That is the average of the academic performance
## in 2000 only for the schools that met the school-wide growth target.
msvymean(~api00, subset(dstrat, sch.wide == "Yes"), k=4, type="rht")
## to extract the estimate from the object 
coef(rht1) 
## to extract the variance from the object
vcov(rht1)



cleanEx()
nameEx("tsvymean")
### * tsvymean

flush(stderr()); flush(stdout())

### Name: tsvymean
### Title: Trimmed and winsorized weighted mean for complex samples
### Aliases: tsvymean

### ** Examples

## load "api" data set from "survey" package (a description of the data
## set can be found there)
data(api)
## define "survey.design" for stratified sampling
dstrat <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, 
fpc=~fpc)
## compute a robust Horvitz-Thompson estimate for the mean of the 
## variable "api00" (Academic Performance Index in 2000)
tm1 <- tsvymean(~api00, dstrat, trim=c(0.01, 0.09), type="trim")
# get a summary of the estimation
summary(tm1)
## robust estimates for a domain of the variable. Here we are 
## interessted in the trimmed mean for api00 in case of 
## (sch.wide == "Yes"). That is the average of the academic performance
## in 2000 only for the schools that met the school-wide growth target.
tsvymean(~api00, subset(dstrat, sch.wide == "Yes"), trim=c(0.01, 0.09), 
type="trim")
## to extract the estimate from the object use
coef(tm1) 
## to extract the variance from the object use
vcov(tm1)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
