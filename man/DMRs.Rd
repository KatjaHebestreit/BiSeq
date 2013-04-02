\name{DMRs}
\alias{DMRs}
\docType{data}
\title{
The output of \code{findDMRs}
}

\description{
  Please see the package vignette for description.
  }

\usage{data(DMRs)}
\format{
A \code{GRanges} of the chromosomes, start and end positions of the
detected DMRs together with information (in the \code{elementMetadata}
slot) on DMRs: \code{median.p}, code{median.meth.group1},
code{median.meth.group2}, \code{median.meth.diff}.
  }
\examples{
data(DMRs)
head(DMRs)
}
\keyword{datasets}
