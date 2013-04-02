\name{plotMeth}
\alias{plotMeth}
\alias{plotMeth,BSraw,BSrel,GRanges-method}

\title{
Plots raw and smoothed methylation data for a given region
}
\description{
This function plots the raw and the smoothed methylation
data for one sample and a given region. The smoothed data is shown as a line (one line
per CpG cluster) and the raw data is shown as points with color
intensities proportional to the coverage.
}
\usage{
plotMeth(object.raw, object.rel, region, col.lines, lwd.lines, col.points, ...)
}
\arguments{
  \item{object.raw}{
A \code{BSraw} with only one sample.
}
  \item{object.rel}{
A \code{BSrel} with only one sample.
}
  \item{region}{
A \code{GRanges} of length one.
}
  \item{col.lines}{
OPTIONAL. The color for the line representing the smoothed
methylation values.
}
  \item{lwd.lines}{
OPTIONAL. The line width for the line representing the smoothed
methylation values.
}
  \item{col.points}{
OPTIONAL. The color for the points representing the raw methylation levels.
}
  \item{\dots}{
Other graphical parameters passed to the \code{plot} function.
}
}


\author{
Katja Hebestreit
}


\seealso{
\code{\link{plotSmoothMeth}}, \code{plot}
}
\examples{

data(rrbs)
data(predictedMeth)

region <- GRanges(seqnames="chr1", 
                  ranges=IRanges(start = 875200, 
                                 end = 875500))
plotMeth(object.raw = rrbs[,6],
         object.rel = predictedMeth[,6],
         region = region)

}
