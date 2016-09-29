\name{mer}
\alias{mer}

\title{Minimum Estimated Risk M-Estimation}

\description{
\code{mer} searches for the robustness-tuning parameter \code{k} (for M-estimation) that minimizes the (inverse-probability weighted) mean square error (MSE). Thus, MER-estimation is a strategy of adaptively choosing the optimal robustness tuning.
}

\usage{
mer(object, init = 0.1, box.lo = 1e-04, tol = 1e-04)
}

\arguments{
  \item{object}{an object of the class \code{"svystat.rob"} (i.e. an estimate of \code{\link{msvymean}} with a first guess of the robustness tuning parameter \code{k})} 
  \item{init}{an initial value of the parameter \code{k} to start the search for an optimum (default 0.1)}
  \item{box.lo}{lower bound (box-constraint) on the variables for the \code{L-BFGS-B} method (default 1e-4; must be \code{>0})}
  \item{tol}{numerical tolerance criterion (delivered to the IRLWS algorithm)}
}

\details{
\code{mer} searches for the robustness tuning parameter \code{k} (for a M-estimator) that minimizes the MSE. The function \code{mer} calls  \code{optim} (in the \pkg{stats} package) to search for an optimal tuning constant \code{k} that minimizes the estimated risk function. Minimization is computed by means of the \code{L-BFGS-B} method (Byrd et al., 1995; Nocedal and Wright, 2006), i.e. a limited-memory modification of the BFGS quasi-Newton method. By default, the following box-constraints are used: lower=1e-4, upper=inf. Note that in typical applications, neither the box-constraints nor the initial value for the parameters to be optimized over, need to be adapted. The algorithm usually converges in a couple of iterations, since it capitalizes (by means of a finite-difference approximation of the gradient) on the almost quadratic shape (at least for symmetric distributions) of the MSE.

\emph{Important notice:} In case of asymmetric distributions, mer-estimation tends to choose optimal tuning constants \code{k} that are far too large. Sometimes the global minimum of the MSE is at zero. In such a case, smaller \code{k} (i.e. downweighting a larger amount of observations) will always reduce the MSE and the optimal M-estimator may be, e.g., the median. 

\emph{Failure of convergence:} If the algorithm failed to converged, set the initial value (\code{init}) of \code{k} near the 'true' k. In addition, you may modify the numeric convergence criterion, \code{tol}.

}

\value{Object of the class(es) \code{"svystat.rob"} and \code{"mer"}. 

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
Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. (1995) A limited memory algorithm for bound constrained 
optimization. \emph{SIAM J. Scientific Computing}, 16, 1190--1208. 

Hulliger, B. (1995): Outlier robust Horvitz-Thompson estimators, \emph{Survey Methodology} 21 (1), pp. 79-87.

Hulliger, B. (1999): Simple and robust estimators for sampling, \emph{Proceedings of the Survey 
Research Methods Section}, American Statistical Association, 1999, pp. 54-63.

Nocedal, J. and Wright, S. J. (2006) \emph{Numerical Optimization}, 2nd. ed. Springer. 
}

\author{Beat Hulliger and Tobias Schoch}

\examples{
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
}

