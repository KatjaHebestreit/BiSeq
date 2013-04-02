\name{trimClusters}
\alias{trimClusters}
\alias{trimClusters,list,missing-method}
\alias{trimClusters,list,numeric-method}


\title{
Trims CpG clusters
}
\description{
CpG clusters rejected in a previous step are trimmed.
}
\usage{
trimClusters(clusters.rej, FDR.loc)
}

\arguments{
  \item{clusters.rej}{Output of \code{testClusters}.}
  \item{FDR.loc}{Location-wise FDR level. Default is 0.02.}
}
\details{
Not differentially methylated CpG sites are removed within the CpG
  clusters rejected by \code{testClusters}.
}
\value{A \code{data.frame} containing the differentially methylated CpG sites.
}
\references{
    Yoav Benjamini and Ruth Heller (2007): False Discovery Rates for Spatial
Signals. American Statistical Association, 102 (480): 1272-81.
}
\author{
Katja Hebestreit
}

\seealso{
\code{\link{testClusters}}
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

clusters.trimmed <- trimClusters(clusters.rej, FDR.loc = 0.05)

}
