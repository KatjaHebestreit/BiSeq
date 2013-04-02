\name{summarizeRegions}
\alias{summarizeRegions}
\alias{summarizeRegions,BSraw,GRanges,logical-method}
\alias{summarizeRegions,BSraw,GRanges,missing-method}
\alias{summarizeRegions,BSrel,GRanges,logical-method}
\alias{summarizeRegions,BSrel,GRanges,missing-method}
  
\title{Aggregates methylation information of single CpG sites}

\description{This method summarizes the methylation
states of single CpG sites to a single methylation state for
a given genomic region.}

\usage{summarizeRegions(object, regions, outputAll)}

\arguments{
  \item{object}{An \code{BSraw} or \code{BSrel} object.}
  \item{regions}{A \code{GRanges} object storing the genomic regions.}
  \item{outputAll}{A logical. If \code{outputAll = TRUE}, all regions
    will be returned. If \code{FALSE} (default), regions are dropped if
    their coverage is zero.}
}

\details{When the given object is of class \code{\link{BSraw-class}}, all
  (methylated) reads of all CpG site lying within a region are summed up and
  assign as total number of (methylated) reads to that region. It is recommended
  to use \code{\link{limitCov}} before applying \code{summarizeRegions} to an
  \code{\link{BSraw-class}} object in order to avoid an excessive influence of a
  single CpG site on the methylation value of a region. When the
  given object is of class \code{\link{BSrel-class}}, the mean relative methylation
  of all CpGs within a region is assign to that region.

  The \code{rowData} slot of the returned object is the given object
  \code{regions} with all columns preserved.
}

\value{An \code{BSraw} or an \code{BSrel} object storing methylation information about
the given regions.}

\author{Hans-Ulrich Klein}

\seealso{
\code{\link{BSraw-class}}, \code{\link{BSrel-class}}, \code{limitCov}
}

\examples{
data(rrbs)
rrbs.clustered <- clusterSites(rrbs)
regions <- clusterSitesToGR(rrbs.clustered)

rrbs <- limitCov(rrbs, maxCov=50)
rrbsRegion <- summarizeRegions(rrbs, regions)
totalReads(rrbsRegion)
}