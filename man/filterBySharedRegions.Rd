\name{filterBySharedRegions}
\alias{filterBySharedRegions}
\alias{filterBySharedRegions,BSraw,ANY,missing,numeric,missing-method}
\alias{filterBySharedRegions,BSraw,ANY,missing,numeric,numeric-method}
\alias{filterBySharedRegions,BSraw,ANY,numeric,missing,missing-method}
\alias{filterBySharedRegions,BSraw,ANY,numeric,missing,numeric-method}
\alias{filterBySharedRegions,BSraw,missing,missing,missing,missing-method}
\alias{filterBySharedRegions,BSraw,missing,missing,missing,numeric-method}
\alias{filterBySharedRegions,BSraw,missing,missing,numeric,missing-method}
\alias{filterBySharedRegions,BSraw,missing,missing,numeric,numeric-method}
\alias{filterBySharedRegions,BSraw,missing,numeric,missing,missing-method}
\alias{filterBySharedRegions,BSraw,missing,numeric,missing,numeric-method}
\alias{filterBySharedRegions,BSrel,ANY,missing,numeric,ANY-method}
\alias{filterBySharedRegions,BSrel,ANY,numeric,missing,ANY-method}
\alias{filterBySharedRegions,BSrel,missing,missing,missing,ANY-method}
\alias{filterBySharedRegions,BSrel,missing,missing,numeric,ANY-method}
\alias{filterBySharedRegions,BSrel,missing,numeric,missing,ANY-method}


\title{Reduces a \code{BSraw} or \code{BSrel} object to regions (or single CpGs) shared by
  a fraction of samples}

\description{This method determines the regions which are covered in a
  given fraction of samples and reduces the \code{BSraw} or \code{BSrel} object to these regions.}

\usage{filterBySharedRegions(object, groups, perc.samples, no.samples, minCov)}

\arguments{
  \item{object}{A \code{BSraw} or \code{BSrel} object.}
  \item{groups}{OPTIONAL. A factor specifying two or more groups within
  the given object. See Details.}
  \item{perc.samples}{A numeric vector with elements between 0 and 1 of
  length 1 or of the same length as levels in group. Default is 1.}
  \item{no.samples}{Alternative to \code{perc.samples}. An integer vector of length 1 or of the same length as levels
    in group.}
  \item{minCov}{A numeric: If \code{object} is a \code{BSraw} object the
  minimum coverage may be set. Default is 1.}
}

\details{If argument \code{group} is given \code{perc.samples} or
  \code{no.samples} are applied for all group levels.}

\value{An object of the same class as \code{object} storing methylation information 
solely for regions or single CpGs covered in at least \code{round(perc.samples*ncol(object))} 
samples, if \code{perc.samples} is given. Alternatively, the number of
samples can be given directly by using \code{no.samples}.}

\author{Katja Hebestreit}

\seealso{\code{\link{filterByCov}}}

\examples{
data(rrbs)
rrbs.reduced <- filterBySharedRegions(object = rrbs, perc.samples = 1)

# Or filter CpG sites where at least 50\% of cancer and 50\% of normal samples have coverage:
rrbs.reduced <- filterBySharedRegions(object=rrbs, groups = colData(rrbs)$group,
                                      perc.samples = 0.5)
}
