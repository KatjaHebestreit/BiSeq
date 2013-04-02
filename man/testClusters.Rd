\name{testClusters}
\alias{testClusters}
\alias{testClusters,list,missing-method}
\alias{testClusters,list,numeric-method}


\title{
Tests CpG clusters
}
\description{
CpG clusters are tested with a cluster-wise FDR level.
}
\usage{
testClusters(locCor, FDR.cluster)
}
\arguments{
  \item{locCor}{Output of \code{estLocCor}.}
  \item{FDR.cluster}{A \code{numeric}. The WFDR (weighted FDR) level at which the CpG clusters
    should be tested. Default is 0.05.}
}
\details{
CpG clusters containing at least one differentially methylated location are detected.
}
\value{A list is returned:
  \item{FDR.cluster }{Chosen WFDR (weighted FDR) for clusters.}
  \item{CpGs.clust.reject }{A \code{list} of the CpG sites together
    with test results within clusters that were rejected.}
  \item{CpGs.clust.not.reject }{A \code{list} of the CpG sites together
    with test results within clusters that were not rejected.}
  \item{clusters.reject }{A \code{GRanges} of the clusters that were rejected.}
  \item{clusters.not.reject }{A \code{GRanges} of the clusters that were
    not rejected.}
  \item{sigma.clusters.reject }{The standard deviations for z-scores
    within each rejected cluster.}
  \item{variogram }{The variogram matrix.}
  \item{m }{Number of clusters tested.}
  \item{k }{Number of clusters rejected.}
  \item{u.1 }{Cutoff point of the largest P value rejected.}
}
\references{
  Yoav Benjamini and Ruth Heller (2007): False Discovery Rates for Spatial
Signals. American Statistical Association, 102 (480): 1272-81.
}
\author{
Katja Hebestreit
}


\seealso{
\code{\link{estLocCor}}, \code{\link{trimClusters}}
}
\examples{
## Variogram under Null hypothesis (for resampled data):
data(vario)

plot(vario$variogram)
vario.sm <- smoothVariogram(vario, sill=0.9)

# auxiliary object to get the pValsList for the test
# results of interest:
data(betaResults)
vario.aux <- makeVariogram(betaResults, make.variogram=FALSE)

# Replace the pValsList slot:
vario.sm$pValsList <- vario.aux$pValsList

## vario.sm contains the smoothed variogram under the Null hypothesis as
## well as the p Values that the group has an effect on DNA methylation.

locCor <- estLocCor(vario.sm)

clusters.rej <- testClusters(locCor, FDR.cluster = 0.1)

}
