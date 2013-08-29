\name{globalTest}
\alias{globalTest}
\alias{globalTest,ANY,BSrel-method}

\title{Test whether at least one CpG is differentially methylated in a
  given genomic region}

\description{This method is a wrapper for conveniently invoking the
  globaltest method \code{gt} on a \code{\link{BSrel-class}} object. The
  globaltest can be applied to test against a high dimensional
  alternative in various regression models. E.g., it can be used to test
  whether at least one CpG is differentially methylated between two
  groups.}

\usage{globalTest(response, alternative, ...)}

\arguments{
  \item{response}{The response vector of the regression model. May be
    supplied as a vector or as a \code{formula} object. In the latter
    case, the right hand side of \code{response} defines the null
    hypothesis. The default null hypopthesis is \code{~1}, i.e. only an
    intercept.}
  \item{alternative}{An object of \code{\link{BSrel-class}} defining the
    alternative. The CpGs are used as explanatory variable in the
    alternative regression model. The null hypothesis is that the
    coefficients of all CpGs are zero. If \code{response} is given as
    \code{formula}, \code{colData(alternative)} is used to obtain the
    respective data.}
  \item{\dots}{Other arguments passed to the \code{\link{gt}}
    method. The argument \code{subsets} can be given as
    \code{\link{GRanges-class}} object. Then, the globaltest is applied
    for each region using only the CpGs lying within the respective
    region. This is useful for, e.g., testing all promoter regions
    within function call.}
}

\details{
  For details see the documentation of the \code{\link{gt}} method in
  package globaltest.
}

\value{
  The function returns an object of class \code{gt.object}. Several
  operations and diagnostic plots for this class are provided by the
  globaltest package.
}

\references{
Goeman, J. J., van de Geer, S. A., and van Houwelingen, J. C. (2006). Testing against
a high-dimensional alternative. Journal of the Royal Statistical Society Series B-
Statistical Methodology, 68(3):477-493.
}

\author{Hans-Ulrich Klein}

\seealso{\code{link{gt}}, \code{link{BSrel}}}

\examples{
data(rrbs)
rrbs <- rawToRel(rrbs)
regions <- GRanges(IRanges(start=c(850000, 1920000, 500), end=c(879000, 1980000, 600)),
                   seqnames=c("chr1", "chr2", "chr3"))

globalTest(group~1, rrbs)
globalTest(group~1, rrbs, subsets=regions)
}