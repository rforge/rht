\name{msvytotal}
\alias{msvytotal}

\title{Robust M-estimation of the total for complex samples}
\description{
	\code{msvytotal} computes robust Horvitz-Thompson estimates of the total or robust weighted total estimates for complex samples by means of M-estimation. 
}

\usage{
msvytotal(y, design, k, type = "rht", na.rm = FALSE, 
         control = rht.control(...), ...)
}

\arguments{
   \item{y}{a formula object (only one variable)}
   \item{design}{a \code{\link[=svydesign]{survey.design}} object}
   \item{k}{robustness tuning constant}
   \item{type}{either \code{"rht"} for robust Horvitz-Thompson estimator (default), or \code{"rwt"} for robust weighted total estimator}
   \item{na.rm}{should cases with missing values be dropped? (default \code{FALSE})}
   \item{control}{control object; see \code{\link{rht.control}}}
   \item{\dots}{(additional specifications which are delivered to \code{rht.control})}
}

\details{
	\code{msvytotal} performs (inverse probability-) weighted M-estimation (Huber psi-function or asymmetric Huber psi-function; \code{asymmetric=TRUE}). The \code{msvytotal} methods supports the following two methods (depending on the underlying survey design)

   \itemize{
      \item robust Horvitz-Thompson total estimator (\code{type="rht"}), 
      \item robust weighted total estimator (\code{type="rwt"}).
   }

   If \code{y} is positively correlated with the inclusion probabilities a "rht" type estimator should be used, and "rwm" otherwise. The initial value is a weighted median or a ratio of weighted medians. You may set \code{steps} equal to one in order to get a one-step estimation. Variance estimates are computed as first-order linearization using the design-based estimation facilities in the \pkg{survey} package.
	
   \code{msvytotal} allows also the estimation for domains. Use the command \code{\link[=subset.survey.design]{subset}} and  a \emph{design subset expression} instead of the original \code{\link[=svydesign]{survey.design}} object in \code{msvytotal} (see examples for more details).

Users may set \code{exact=TRUE} to compute an "exact" linearization-variance estimate, which takes into account that the MAD has been used as preliminary scale estimate. However, the estimates may become very unstable. 
}

\value{
Object of class \code{"svystat.rob"}

The following (S3) methods are defined for objects of the class \code{"svystat.rob"}:
   \itemize{
      \item \code{\link[=print.svystat.rob]{print}} method,
      \item \code{\link[=summary.svystat.rob]{summary}} method,
      \item \code{\link[=coef.svystat.rob]{coef}} method,
      \item \code{\link[=vcov.svystat.rob]{vcov}} method,
      \item \code{\link[=residuals.svystat.rob]{residuals}} method,
      \item \code{\link[=robweights]{robweights}} method.
   }
}

\references{
Hulliger, B. (1995): Outlier robust Horvitz-Thompson estimators, \emph{Survey Methodology} 21 (1), pp. 79-87.

Hulliger, B. (1999): Simple and robust estimators for sampling, \emph{Proceedings of the Survey 
Research Methods Section}, American Statistical Association, 1999, pp. 54-63.

Hulliger, B. and T. Schoch (2011): Elementary robust estimators. In: \emph{Robust methodology for Laeken indicators: AMELI Deliverable D4.2}, ed. by B. Hulliger, A. Alfons, P. Filzmoser, A. Meraner, T. Schoch and M. Templ. AMELI Project.


}

\author{Beat Hulliger and Tobias Schoch}

\seealso{\code{\link{svytotal}}
}

\examples{
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
}
