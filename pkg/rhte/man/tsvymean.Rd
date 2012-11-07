\name{tsvymean}
\alias{tsvymean}

\title{Trimmed and winsorized weighted mean for complex samples}
\description{
	\code{tsvymean} computes either the trimmed or winsorized weighted mean for complex samples. 
}
\usage{
tsvymean(y, design, trim = c(0, 0), type = c("trim", "win"), 
         na.rm=FALSE, control = rht.control(...), ...)
}

\arguments{
   \item{y}{a formula object (only one variable)}
   \item{design}{a \code{\link[=svydesign]{survey.design}} object}
   \item{trim}{vector of size two consisting of the lower and upper amount of trimming (\code{[lo, hi]}). The fraction \code{lo} of observations is trimmed from the lower end and the fraction \code{hi} is trimmed from the upper end. The symmetrically 5\%-trimmed mean is obtained with \code{trim=c(0.05, 0.05)}. Alternatively, one may specify the amount of trimming in terms of integer values; e.g., \code{trim=c(0, 6)} trims only the largest 6 observations.}
   \item{type}{either \code{"trim"} for trimming (default), or \code{"win"} for winsorization}
   \item{na.rm}{should cases with missing values be dropped? (default \code{FALSE})}
   \item{control}{control object; see \code{\link{rht.control}}}
   \item{\dots}{additional specifications}
}

 \details{
  By default \code{trim} equals \code{c(0,0)} and the regular weighted mean is computed. The variance estimators are based on first-order linearizations using the design-based-estimation facilities of the \pkg{survey} package. For reasons of numerical stability, the variance of the winsorized weighted mean is computed using the variance estimator of the trimmed mean.
  	
\code{tsvymean} allows also the estimation for domains. Use the command \code{\link[=subset.survey.design]{subset}} and a \emph{design subset expression} instead of the original \code{\link[=svydesign]{survey.design}} object in \code{tsvymean} (see examples for more details).

}

\value{
Object of class \code{"svystat.rob"}
 
The following (S3) methods are defined for objects of the class \code{"svystat.rob"}:
   \itemize{
      \item \code{\link[=print.svystat.rob]{print}} method,
      \item \code{\link[=summary.svystat.rob]{summary}} method,
      \item \code{\link[=coef.svystat.rob]{coef}} method,
      \item \code{\link[=vcov.svystat.rob]{vcov}} method.
   }

}

\references{Hulliger, B. (1999): Simple and robust estimators for sampling, \emph{Proceedings of the Survey Research Methods Section}, American Statistical Association, 1999, pp. 54-63.

Hulliger, B. and T. Schoch (2011): Elementary robust estimators. In: \emph{Robust methodology for Laeken indicators: AMELI Deliverable D4.2}, ed. by B. Hulliger, A. Alfons, P. Filzmoser, A. Meraner, T. Schoch and M. Templ. AMELI Project.

}

\author{Tobias Schoch}

\seealso{\code{\link{svymean}}
}

\examples{
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
}
