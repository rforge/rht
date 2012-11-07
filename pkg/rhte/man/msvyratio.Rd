\name{msvyratio}
\alias{msvyratio}

\title{Robust ratio M-estimation for complex samples}

\description{
	\code{msvyratio} computes a robust ratio estimate for complex samples by means of M-estimation ("rht" type).}              

\usage{
msvyratio(numerator, denominator, design, k, na.rm = FALSE, 
          control = rht.control(...), ...)
}

\arguments{
   \item{numerator}{a formula object (only one variable)}
   \item{denominator}{a formula object (only one variable)}
   \item{design}{a \code{\link[=svydesign]{survey.design}} object}
   \item{k}{robustness tuning constant}
   \item{na.rm}{should cases with missing values be dropped? (default \code{FALSE})}
   \item{control}{control object; see \code{\link{rht.control}}}
   \item{\dots}{(additional specifications which are delivered to \code{rht.control})}
}

\details{
	\code{msvyratio} computes a robust ratio estimate for complex samples by means of M-estimation (type "rht"; see \code{\link{msvymean}} for more details). Variance estimates are computed as first-order linearization using the design-based estimation facilities in the \pkg{survey} package.
	
   The initial value is a weighted median or a ratio of weighted medians. You may set \code{steps} equal to one in order to obtain an one-step estimator.
	
   \code{msvyratio} allows also the estimation for domains. Use the command \code{\link[=subset.survey.design]{subset}} and a \emph{design subset expression} instead of the original \code{\link[=svydesign]{survey.design}} object in \code{msvyratio} (see examples for more details).
 
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

\seealso{\code{\link{svyratio}}
}

\examples{
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
}
