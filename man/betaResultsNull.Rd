\name{betaResultsNull}
\alias{betaResultsNull}
\docType{data}
\title{
The output of \code{betaRegression} for resampled data
}

\description{
  Please see the package vignette for description.
  }

\usage{data(betaResultsNull)}
\format{
  A data frame with 4276 observations on the following 10 variables.
  \describe{
    \item{\code{chr}}{a factor with levels \code{chr1} \code{chr2}}
    \item{\code{pos}}{a numeric vector}
    \item{\code{p.val}}{a numeric vector}
    \item{\code{meth.group1}}{a numeric vector}
    \item{\code{meth.group2}}{a numeric vector}
    \item{\code{meth.diff}}{a numeric vector}
    \item{\code{estimate}}{a numeric vector}
    \item{\code{std.error}}{a numeric vector}
    \item{\code{pseudo.R.sqrt}}{a numeric vector}
    \item{\code{cluster.id}}{a character vector}
  }
}
\examples{
data(betaResultsNull)
head(betaResultsNull)
}
\keyword{datasets}
