\name{makeVariogram}
\alias{makeVariogram}
\alias{makeVariogram,data.frame,logical-method}
\alias{makeVariogram,data.frame,missing-method}


\title{
Variogram estimator.
}
\description{
A function which estimates the variogram of the z-scores in the given data frame.
}
\usage{
makeVariogram(test.out, make.variogram)
}
\arguments{
  \item{test.out}{
A \code{data.frame}. Usually the output of \code{betaRegression}. Must contain columns chr, pos, p.val and cluster.id.
}
  \item{make.variogram}{A \code{logical}. Default is \code{TRUE}.}
}
\details{
For each CpG site the z-score is determined by \code{qnorm}(1 -
  P value). The variogram of the z-scores of locations \eqn{k}
and \eqn{l} within one cluster is estimated robustly by
\deqn{
2 \hat{\gamma}(h) = [median{(Z_{k}-Z_{l})^2: (s_{k}, s_{l}) \in
  N(h)}] / .455 
  }.
}
\value{
  A list:
  \item{variogram }{A \code{matrix} with columns \code{h} and \code{v}.}
  \item{pValsList }{A \code{list} of data frames. Each data frame corresponds
    to a CpG cluster and contains same information as \code{test.out}
    plus the columns z.score and pos.new (position corresponding to the
    respective CpG cluster).}
}
\references{
Yoav Benjamini and Ruth Heller (2007): False Discovery Rates for Spatial
Signals. American Statistical Association, 102 (480): 1272-81.
}
\author{
Katja Hebestreit
}

\seealso{
  \code{\link{betaRegression}}
  }

\examples{
data(betaResults)

vario <- makeVariogram(betaResults)

plot(vario$variogram)
}
