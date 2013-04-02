\name{betaRegression}
\alias{betaRegression}
\alias{betaRegression,formula,character,BSrel,numeric-method}
\alias{betaRegression,formula,character,BSrel,missing-method}

\title{
A function to estimate and test a group factor within a beta regression
}
\description{
This function models the methylation level within a beta regression. The
first independent variable in \code{formula} is tested to be
unequal to zero.
}
\usage{
betaRegression(formula, link, object, mc.cores, ...)
}

\arguments{
  \item{formula}{
    Symbolic description of the model. For the first independent
  variable the P value (Wald test) and the effect on methylation is returned. For
  details see below.
  }
  \item{link}{
    A character specifying the link function in the mean
    model (mu). Currently, \code{"logit"}, \code{"probit"}, \code{"cloglog"},
    \code{"log"}, \code{"loglog"} are supported.
}
  \item{object}{
    A \code{BSrel} object.
}
  \item{mc.cores}{
    Passed to \code{mclapply}.
}
  \item{\dots}{
    Other parameters passed to the \code{betareg} function.
}
}
\details{
  See \code{betareg} function for details.

  \code{mclapply}
}
\value{
  A \code{data.frame} containing the position, chromosome, P value, estimated
methylation level in group 1 and group 2 and methylation difference of
group 1 and group 2.
}
\references{
  Hebestreit, K., Dugas, M., Mueller-Tidow, C., and Klein HU
(2013). Detection of Significantly Differentially Methylated Regions in
Reduced Representation Bisulfite Sequencing Data. In preparation.
  
  See also reference list in the documentation of \code{\link{betareg}}.
}

\author{
  Katja Hebestreit
}


\seealso{
\code{betareg}
}
\examples{
# load RRBS data, subset to save time, find CpG clusters and smooth methylation data:
data(rrbs)
rrbs.small <- rrbs[1:1000,]
rrbs.clust.unlim <- clusterSites(object = rrbs.small,
                                 groups = colData(rrbs)$group,
                                 perc.samples = 4/5,
                                 min.sites = 20,
                                 max.dist = 100)

ind.cov <- totalReads(rrbs.clust.unlim) > 0
quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov], 0.9)
rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)

# with a small subset to save calculation time:
rrbs.part <- rrbs.clust.lim[1:100,]

predictedMeth <- predictMeth(object=rrbs.part)

betaResults <- betaRegression(formula = ~group, link = "probit",
                              object = predictedMeth, type="BR")
}
