\name{pmu284}
\alias{pmu284}
\docType{data}
\title{
The MU284 Population
}
\description{
The MU284 population comes from Sarndal, Swensson, and Wretman (1992, Appendix B) and can be downloaded from Stablib. The data consists of readings from 284 Swedish municipalites. Typically, a Swedish municipality consists of a town and the surrounding area. 
}
\usage{data(pmu284)}
\format{
  A data frame with 284 observations on the following 11 variables.
  \describe{
    \item{\code{LABEL}}{municipality identifier, running from 1 to 284}
    \item{\code{P85}}{1985 population in thousands}
    \item{\code{P75}}{1975 population on thousands}
    \item{\code{RMT85}}{revenues from the 1985 municipal taxation (in millions of kronor)}
    \item{\code{CS82}}{number of Conservative seats in municipal council}
    \item{\code{SS82}}{number of Social-Democratic seats in municipal council}
    \item{\code{S82}}{total number of seats in municipal council}
    \item{\code{ME84}}{number of municipal employees in 1984}
    \item{\code{REV84}}{real estate values according to 1984 assessement (in millions of kronor)}
    \item{\code{REG}}{geographic region indicator}
    \item{\code{CL}}{cluster indicator (a cluster consists of a set of neighboring municipalities)}
  }
}
\details{
See Sarndal, Swensson, and Wretman (1991, Appendix B) for more details
}
\source{
Carl Erik Sarndal, Bengt Swensson, and Jan Wretman (1991). \emph{Model Assisted Survey Sampling}. New York: Springer.
}
\references{
Full MU284 population at \url{http://lib.stat.cmu.edu/datasets/mu284}
}
\keyword{datasets}
