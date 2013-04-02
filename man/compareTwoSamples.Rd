\name{compareTwoSamples}
\alias{compareTwoSamples}
\alias{compareTwoSamples,BSrel,numeric,numeric,numeric,numeric-method}
\alias{compareTwoSamples,BSrel,character,character,numeric,numeric-method}
\alias{compareTwoSamples,BSraw,numeric,numeric,numeric,numeric-method}
\alias{compareTwoSamples,BSraw,character,character,numeric,numeric-method}

\title{
Detects DMRs by comparing two samples
}
\description{
Determines the differences of (smoothed) methylation levels between two
samples and aggregates the sites surpassing a minimum difference to DMRs.
}
\usage{
compareTwoSamples(object, sample1, sample2, minDiff, max.dist)
}

\arguments{
  \item{object}{
A \code{BSrel}.
}
  \item{sample1}{
A numeric or character specifying the first sample to be used.
}
  \item{sample2}{
A numeric or character specifying the second sample to be used.
}
  \item{minDiff}{
A numeric greater than 0 and smaller or equal to 1.
}
  \item{max.dist}{
Numeric. The maximum distance between two CpG sites (or grid points)
with absolute methylation differences greater or equal
than \code{minDiff} in a DMR. If grid points are used: should be the same as grid.dist in
\code{predictMeth}.
}
}

\details{
This function determines the differences between the methylation levels of sample1 and sample2 for each site. Successive sites with methylation differences smaller or equal to \code{minDiff} are summarized.
}
\value{
A GRanges object.
}

\author{
Katja Hebestreit
}

\seealso{
\code{\link{predictMeth}}
}

\examples{
data(rrbs)
rrbs <- rrbs[, c(1,6)]
CpG.clusters <- clusterSites(object = rrbs, perc.samples = 1,
                             min.sites = 20, max.dist = 100)
predictedMeth <- predictMeth(object = CpG.clusters)
DMRs <- compareTwoSamples(predictedMeth, sample1 = 1, sample2 = 2,
                          minDiff = 0.3, max.dist = 100)
}
