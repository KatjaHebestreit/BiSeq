\name{clusterSitesToGR}
\alias{clusterSitesToGR}
\alias{clusterSitesToGR,BSraw-method}
\alias{clusterSitesToGR,BSrel-method}

\title{
A function to obtain a \code{GRanges} object of CpG clusters from  \code{BSraw} and \code{BSrel} objects
}
\description{
This function allows to get the start and end positions of CpG clusters
from a \code{BSraw} or \code{BSrel} object, when there is a \code{cluster.id}
column in the \code{rowData} slot.
}
\usage{
clusterSitesToGR(object)
}
\arguments{
  \item{object}{A \code{BSraw} or \code{BSrel} object with a
    \code{cluster.id} column in the \code{rowData} slot. Usually the
    output of \code{clusterSites}.}
}

\value{
An object of class \code{GRanges} is returned.
}

\author{
Katja Hebestreit
}


\seealso{
  \code{\link{clusterSites}}
}

\examples{
data(rrbs)
rrbs.clustered <- clusterSites(rrbs)

clusterSitesToGR(rrbs.clustered)
}