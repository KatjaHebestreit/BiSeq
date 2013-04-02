\name{plotBindingSites}
\alias{plotBindingSites}
\alias{plotBindingSites,BSraw,GRanges-method}
\alias{plotBindingSites,BSrel,GRanges-method}

\title{Plots the mean methylation of given regions}

\description{
\code{plotBindingSites} takes several genomic regions (e.g. protein binding
sites), centers them such that the position 0 refers to the center of each region
and finally calculates the mean methylation of all regions for each given sample.
If several samples are given, the median of the samples' methylation values and
optionally other quantiles are plotted.
}

\usage{
plotBindingSites(object, regions, width, groups, quantiles, bandwidth, ...)
}

\arguments{
  \item{object}{An object of class \code{BSraw} or \code{BSrel}.}
  \item{regions}{Regions given a \code{GRanges} object. The regions may have different widths.}
  \item{width}{The width of the genomic region that is plotted. Default value is the width of the largest given region.}
  \item{groups}{An optional factor defining two or more groups within the given object. The mean methylation is than plotted for each group separately.}
  \item{quantiles}{Other quantiles to be plotted besides the median. Default are the 25\% and the 75\% quantiles.}
  \item{bandwidth}{The bandwidth of the kernel smoother used for smoothing methylation values. Default value is 1/8 \code{width}.}
  \item{\dots}{Other graphical parameters passed to the \code{plot} function.}
}

\details{
First, all regions were expanded or shrinked to the given \code{width} by adding
or removing base pairs symmetrically at both ends of the regions (not by scaling).
A new coordinate system is centered at the middle of the equally sized regions.
Next, the relative methylation values for each sample are averaged accross all
regions. That means, if there are several CpGs from different regions lying the
same position, the mean methylation value is calculated for that position. Then,
the median of these methylation values across all samples is calculated.
Optionally, other quantiles are calculated, too. The median of the methylation is
then plotted for each position after smoothing using a gaussian kernel with the
given bandwidth.

If the given regions correspond to binding sites of a certain protein, the plot
can be used to discover whether the protein induces changes in the DNA methylation
in the proximity of its binding sites.
}

\author{Hans-Ulrich Klein}

\seealso{\code{\link{BSraw-class}}, \code{\link{BSrel-class}}}

\examples{
data(rrbs)
data(promoters)
plotBindingSites(object=rrbs,
                 regions=promoters,
                 width=4000, 
                 groups=colData(rrbs)$group)
}
