\name{rht.control}
\alias{rht.control}
\title{
Control function for rht
}
\description{
This function specifies the parameters needed to compute robust estimators 
}
\usage{
rht.control(acc = 1e-06, steps = 50, asymmetric = FALSE, 
            quietly = FALSE, exact = FALSE, ddigits = 3, ...)
}

\arguments{
  \item{acc}{
scalar; determines the precision in the termination rule of the IRWLS algorithm 
}
  \item{steps}{
scalar (integer); determines the maximal number of IRWLS iterations
}
  \item{asymmetric}{
logical; \code{FALSE} gives symmetric (standard) Huber psi-function (default); \code{TRUE} gives the asymmetric Huber psi-function
}
  \item{quietly}{
logical; default is \code{FALSE}; if \code{quietly} is \code{TRUE}, the text/message output (stdout) of all algorithms is suppressed (this is useful in simulation studies)
}
  \item{exact}{
locial; if \code{TRUE} the exact asymptotic variance of the M-estimator is computed, where exact means that the variance of the scale estimate is considered; if \code{FALSE}, the contribution of the scale variance is neglected (default). 
}
  \item{ddigits}{
scalar (integer); determines the number of relevant decimal digits to be reported (default=3)
}
  \item{\dots}{additional specifications
}
}
\details{
The user may specify parameters differently if needed. Notably, if the algorithm fails to converge, one may consider modifying the number of iterations (\code{steps}) and/or the precision (\code{acc}).
}
\value{
A list consisting of the parameter values.
}
\author{
Tobias Schoch
}

