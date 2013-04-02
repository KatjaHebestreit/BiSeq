\name{plotSmoothMeth}
\alias{plotSmoothMeth}
\alias{plotSmoothMeth,BSrel,GRanges,ANY,logical-method}
\alias{plotSmoothMeth,BSrel,GRanges,ANY,missing-method}


\title{
Plots smoothed methylation values for a bunch of samples and a given region
}
\description{
This function plots the smoothed methylation data as lines for a given region and all given
samples. It is also possible to average the data for groups of
samples.
}
\usage{
plotSmoothMeth(object.rel, region, groups, group.average, ...)
}
\arguments{
  \item{object.rel}{
A \code{BSrel}.
}
  \item{region}{
A \code{GRanges} of length one.
}
  \item{groups}{
OPTIONAL. A \code{factor} defining two or more sample groups within the given
object.
}
  \item{group.average}{
OPTIONAL. A \code{logical}. If \code{TRUE}, then the data is averaged for the groups
given in \code{groups}. Default is \code{FALSE}.. 
}
  \item{\dots}{
Other graphical parameters passed to the \code{plot} function.
}
}

\author{
Katja Hebestreit
}



\seealso{
  \code{\link{plotMeth}}, \code{plot}
}
\examples{

data(predictedMeth)
data(DMRs)

plotSmoothMeth(object.rel = predictedMeth,
               region = DMRs[3] + 200,
               groups = colData(predictedMeth)$group,
               col=c("magenta", "blue"))
legend("topright",
       lty=1,
       legend=levels(colData(predictedMeth)$group),
       col=c("magenta", "blue"))

}
