\name{filterByCov}
\alias{filterByCov}
\alias{filterByCov,BSraw,missing,logical-method}
\alias{filterByCov,BSraw,missing,missing-method}
\alias{filterByCov,BSraw,numeric,logical-method}
\alias{filterByCov,BSraw,numeric,missing-method}



\title{
Filters regions (or single CpGs) of a \code{BSraw} object with a minimum
coverage
}
\description{
This method reduces a \code{BSraw} object to its regions (or single CpGs) with a minimum number of reads.
}
\usage{
filterByCov(object, minCov, global)
}
\arguments{
  \item{object}{A \code{BSraw}.}
  \item{minCov}{Minimum number of reads overlapping the CpG sites.}
  \item{global}{A logical indicating whether the regions should achieve
    the minimum coverage in each sample. If \code{global = TRUE} the
    filtered object will consist of the regions achieving the minimum
    coverage in all samples. If \code{global = FALSE} (default) this function filters the regions for each sample separately, irrespectively of the coverages in other samples. \code{totalReads} and \code{methReads} are set to zero, if the minimum coverage is not obtained. Regions covered too sparse in all samples are dropped.}
}

\value{A \code{BSraw} object containing the CpGs or regions achieving
  the minimum coverage in all (if \code{global=TRUE}) or at least one
  (if \code{global=FALSE}) samples.}

\author{
Katja Hebestreit
}

\seealso{\code{\link{filterBySharedRegions}}}

\examples{
data(rrbs)
rrbs.reduced <- filterByCov(object=rrbs, minCov=10, global=TRUE)
}
