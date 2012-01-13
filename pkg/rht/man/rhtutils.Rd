\name{rht-utils}
\alias{summary.svystat.rob}
\alias{print.svystat.rob}
\alias{coef.svystat.rob}
\alias{vcov.svystat.rob}
\alias{residuals.svystat.rob}
\alias{robweights}

\title{rht utility functions}
\description{
The \code{rht} package contains some useful utility functions to extract relevant information from objects of the class \code{"svystat.rob"}.
}

\usage{
robweights(object)
\method{print}{svystat.rob}(x, ...)
\method{summary}{svystat.rob}(object, ...)
\method{coef}{svystat.rob}(object, ...)
\method{vcov}{svystat.rob}(object, ...)
\method{residuals}{svystat.rob}(object, ...)
}

\arguments{
   \item{x}{(only used by the \code{print} method)}
   \item{object}{(only used by the \code{summary}, \code{coef}, and \code{vcov} method)}
   \item{\dots}{additional specifications}
}


\details{

   \itemize{
      \item \code{summary} is a method to summarize the object
      \item \code{coef} is a method to extract the estimates (coefficients) from a \code{svystat.rob} object
      \item \code{vcov} is a method to extract the variance from a \code{svystat.rob} object
      \item \code{print} is a print method
      \item \code{residuals} extracts the residuals from a \code{svystat.rob} object
      \item \code{robweights} extracts the robustness weights from an M-estimate
   }
}

\author{Beat Hulliger and Tobias Schoch}

\examples{
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
}
