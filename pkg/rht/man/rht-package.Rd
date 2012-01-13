\name{rht-package}
\alias{rht-package}
\alias{rht}
\docType{package}
\title{Robust Mean and Ratio Estimation for Complex Samples}
\description{Robust mean and ratio estimation in complex samples based on M-estimation, trimming, and winsorization.}

\details{
\tabular{ll}{
Package: \tab rht\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2011-11-04\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}


   \subsection{Implemented methods}{
      \itemize{
	 \item Horvitz--Thompson M-estimator (total and mean)
	 \item Ratio M-estimator
	 \item Trimmed and winsorized estimation of the mean
      }
   }

   \subsection{Dependency}{
   The package is based on the \pkg{survey} package (Lumley, 2004, 2011).

   }
}

\author{Beat Hulliger and Tobias Schoch

Maintainer: Tobias Schoch <tobias.schoch@fhnw.ch>
}

\references{
Thomas Lumley (2011). \emph{survey: Analysis of complex survey samples}. R package version 3.26

Thomas Lumley (2004). Analysis of complex survey samples. \emph{Journal of Statistical Software} 9(1), 1--19.
}
\keyword{ survey }
\seealso{\code{survey}}
