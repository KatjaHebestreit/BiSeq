\name{findDMRs}
\alias{findDMRs}
\alias{findDMRs,data.frame,missing,numeric,logical-method}
\alias{findDMRs,data.frame,missing,numeric,missing-method}
\alias{findDMRs,data.frame,numeric,numeric,logical-method}
\alias{findDMRs,data.frame,numeric,numeric,missing-method}

\title{Aggregates CpG sites to DMRs}

\description{This function aggregates CpG sites to DMRs on the basis of
  their P values.}

\usage{
findDMRs(test.out, alpha, max.dist, diff.dir)
}

\arguments{
  \item{test.out}{An object returned by \code{betaRegression}.}
  \item{alpha}{OPTIONAL. A DMR contains CpG sites with P values smaller or equal than \code{alpha}.}
  \item{max.dist}{Numeric. The maximum distance between two P values smaller than 
\code{alpha} in a DMR. Should be the same as \code{grid.dist} in
\code{predictMeth}.}
 \item{diff.dir}{Logical. Should DMRs be seperated if the direction of methylation
  differences changes? If \code{TRUE} (default), than resulting DMRs
  will consist of sites which are all hypomethylated, or hypermethylated respectively.}
}


\value{A \code{GRanges} object storing the start and end positions of the DMRs with information in \code{elementMetadata}:
  \item{median.p }{median of P values}
  \item{median.meth.group1 }{median of modeled methylation level of group1.}
  \item{median.meth.group1 }{median of modeled methylation level of group2.}
  \item{median.meth.diff }{median of difference of modeled methylation levels of group1 and group2.}

}



\author{Katja Hebestreit}

\seealso{
  \code{\link{predictMeth}}, \code{\link{betaRegression}}
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

DMRs <- findDMRs(clusters.trimmed, max.dist=100, diff.dir=TRUE)

}
