\name{estLocCor}
\alias{estLocCor}
\alias{estLocCor,list-method}

\title{
Estimates the correlations of the z-scores
}
\description{
For each location the correlation of this location's z-score to \eqn{\bar{Z}}
of its CpG cluster is estimated.
}
\usage{
estLocCor(vario.sm)
}
\arguments{
  \item{vario.sm}{Output of \code{smoothVariogram}.}
}

\value{A list:
  \item{variogram }{A variogram matrix, usually created by
    \code{smoothVariogram} beforehand.}
  \item{pValsList }{A list of the test results per CpG cluster.}
  \item{sigma.cluster }{The standard deviations of z-scores
    within each cluster.}
  \item{Z.cluster }{The arithmetic means of the z-scores for each cluster.}
  \item{length.cluster }{The widths (number of pase pairs) of each cluster.}
}
\references{
  Yoav Benjamini and Ruth Heller (2007): False Discovery Rates for Spatial
Signals. American Statistical Association, 102 (480): 1272-81.
}
\author{
Katja Hebestreit
}

\seealso{
  \code{\link{makeVariogram}}, \code{\link{smoothVariogram}}
  }

  \examples{
data(betaResultsNull)

vario <- makeVariogram(betaResultsNull)
vario.sm <- smoothVariogram(vario, sill = 1)

locCor <- estLocCor(vario.sm)
}