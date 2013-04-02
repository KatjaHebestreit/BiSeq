\name{limitCov}
\alias{limitCov}
\alias{limitCov,BSraw,numeric-method}

\title{Limits the coverage of a \code{BSraw} object}
\description{Number of methylated and unmethylated reads 
of a CpG site with coverage above \code{maxCov} are reduced such that
the methylation level remains unchanged.}
\usage{
limitCov(object, maxCov)
}

\arguments{
  \item{object}{A \code{BSraw}.}
  \item{maxCov}{The maximum number of reads a CpG should have. All
  coverages above this threshold are limited. (Default is 50)}
}
\details{This function might be useful prior to the use of 
\code{predictMeth} to limit the weights of CpGs with extremly
  high coverages. See \code{binomLikelihoodSmooth} for details.}
\value{A \code{BSraw} object.}

\author{
Katja Hebestreit
}

\seealso{
  \code{\link{predictMeth}}, \code{\link{binomLikelihoodSmooth}}
  }


\examples{
data(rrbs)

rrbs.clust.unlim <- clusterSites(object = rrbs,
                                 groups = colData(rrbs)$group,
                                 perc.samples = 4/5,
                                 min.sites = 20, max.dist = 100)
covBoxplots(rrbs.clust.unlim)

# 90\% quantile of coverage is 39x
quantile(totalReads(rrbs.clust.unlim)[totalReads(rrbs.clust.unlim)>0],
0.9)

rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = 39)
covBoxplots(rrbs.clust.lim)
}
