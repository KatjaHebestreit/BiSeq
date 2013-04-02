\name{plotMethMap}
\alias{plotMethMap}
\alias{plotMethMap,BSraw,GRanges,factor,logical-method}
\alias{plotMethMap,BSraw,GRanges,factor,missing-method}
\alias{plotMethMap,BSraw,GRanges,missing,logical-method}
\alias{plotMethMap,BSraw,GRanges,missing,missing-method}
\alias{plotMethMap,BSrel,GRanges,factor,logical-method}
\alias{plotMethMap,BSrel,GRanges,factor,missing-method}
\alias{plotMethMap,BSrel,GRanges,missing,logical-method}
\alias{plotMethMap,BSrel,GRanges,missing,missing-method}


\title{Plots methylation values of multiple samples in a given region}

\description{A heatmap like plot is generated showing the relative
methylation of single CpG sites. Samples are clustered hierarchically.}

\usage{plotMethMap(object, region, groups, intervals, ...)}

\arguments{
  \item{object}{A \code{BSraw} or \code{BSrel} object storing the methylation values.}
  \item{region}{A \code{GRanges} object giving the region of interest.}
  \item{groups}{OPTIONAL. A \code{factor} that will be encoded by a color
    bar.}
  \item{intervals}{OPTIONAL. A \code{logical} indicating whether
    neighboured CpG sites should be placed side by side (if \code{FALSE})
    or whether the intervals between CpG sites should be preserved
    (if \code{TRUE}).}
  \item{\dots}{Further arguments passed to the \code{heatmap} function.}
}

\details{
The relative methylation values are passed to the heatmap function. Default
colors are green (not methylated), black and red (methylated). To ensure
that a relative methylation of 0 corresponds to green, 0.5 to black and
1 to red, the default value for the \code{zlim} argument of the
\code{\link{heatmap}} function is set to \code{c(0,1)}. And the default
for the \code{scale} parameter is set to \code{"none"}.\cr
If argument \code{intervals} is set to \code{TRUE}, region should not
be too large (< 1kb) and respect the resolution of your screen.
}

\author{Hans-Ulrich Klein}

\seealso{
\code{heatmap}, \code{\link{BSraw-class}}, \code{\link{BSrel-class}}, \code{\link{filterBySharedRegions}},
\code{\link{filterByCov}}
}

\examples{
data(rrbs)
data(predictedMeth)
data(DMRs)

plotMethMap(rrbs, region = DMRs[4], groups = colData(rrbs)[, "group"])

plotMethMap(predictedMeth, region = DMRs[4],
            groups = colData(rrbs)[,"group"], intervals = FALSE)
}
