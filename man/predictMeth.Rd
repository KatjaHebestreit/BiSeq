\name{predictMeth}
\alias{predictMeth}
\alias{predictMeth,BSraw,missing,missing,missing-method}
\alias{predictMeth,BSraw,missing,missing,numeric-method}
\alias{predictMeth,BSraw,numeric,missing,missing-method}
\alias{predictMeth,BSraw,numeric,missing,numeric-method}
\alias{predictMeth,BSraw,numeric,numeric,numeric-method}



\title{Predicts methylation levels along CpG sites or for a grid of sites in CpG clusters.}
\description{Uses local regression to predict methylation levels per sample.}
\usage{
predictMeth(object, h, grid.dist, mc.cores)
}

\arguments{
  \item{object}{A \code{BSraw} with a \code{cluster.id} column in the
    \code{elementMetadata}, usually the output of \code{clusterSites}.}
  \item{h}{Bandwidth in base pairs. Large values produce a smoother
  curve. Default is 80.}
  \item{grid.dist}{OPTIONAL. If \code{numeric}, than methylation values are
  predicted at intervals of \code{grid.dist} base pairs. By default,
  methylation is smoothed at each CpG site.}
  \item{mc.cores}{
    Passed to \code{mclapply}. Default is 1.
  }
}
\details{Uses \code{binomLikelihoodSmooth} with \code{pos} = CpG position, \code{m} = number
  methylated reads and \code{n} = number of reads. \code{pred.pos}
  corresponds to all CpG positions, or to the grid sites respectively, within the CpG clusters.}
\value{A \code{BSrel} object containing the predicted methylation levels in
  the \code{methLevel} slot.}

\author{
Katja Hebestreit
}

\seealso{
\code{\link{clusterSites}}, \code{\link{binomLikelihoodSmooth}}, \code{mclapply}
}


\examples{
data(rrbs)

rrbs.clust.unlim <- clusterSites(object = rrbs,
                                 groups = colData(rrbs)$group,
                                 perc.samples = 4/5,
                                 min.sites = 20, max.dist = 100)

ind.cov <- totalReads(rrbs.clust.unlim) > 0
quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov], 0.9)
rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)

# with a small subset to save calculation time:
rrbs.part <- rrbs.clust.lim[1:100,]

predictedMeth <- predictMeth(object=rrbs.part)
}
