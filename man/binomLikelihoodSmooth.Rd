\name{binomLikelihoodSmooth}
\alias{binomLikelihoodSmooth}
\alias{binomLikelihoodSmooth,ANY-method}

\title{
Calculates local likelihood estimations for binomial random variables
}

\description{
For a given set of binomial random variables with 1-dimensional
coordinates, this function calculates the local likelihood estimation of
the success probability p at a given point. For this purpose, a weighted
likelihood estimation with weights obtained by a triangular kernel with
given bandwidth is used. This can be used to predict values at points
where no variable has been observed and/or to smooth observations using
neighboured observations.
}

\usage{
binomLikelihoodSmooth(pred.pos, pos, m, n, h)
}

\arguments{
\item{pred.pos}{A vector of positions where p should be estimated.}
\item{pos}{A vector of positions where binomial variables have been observed.}
\item{m}{A vector of length \code{pos} with the number of successfull experiments.}
\item{n}{A vector of length \code{pos} with the number of experiments.}
\item{h}{The bandwidth of the kernel.}
}

\details{
For a given position x, the weighted likelihood for parameter \eqn{p}
\deqn{L(p; m, n, w) = \prod_{i=1}^k B(m_i|n_i,p)^{w_i}}
is maximized. B denotes the binomial probability function. The weights
\eqn{w_i} are calculated using a triangular kernel with bandwidth
\eqn{h}:
\deqn{w_i = K(x_i) = (1 - (|x - x_i|)/h) \mathbf{1}_{((|x - x_i|)/h) \le 1}}
}

\value{
A vector of length \code{pred.pos} giving the local likelihood estimation of
the success probability p at the given positions.
}

\author{
Hans-Ulrich Klein
}

\seealso{
\code{\link{predictMeth}}
}

\examples{
n = rpois(100, lambda=10)
E = c(rep(0.4, 30), rep(0.8, 40), rep(0.1, 30))
m = rbinom(100, n, E)
pos = 1:100
p_10 = binomLikelihoodSmooth(pos, pos, m, n, h=10)
p_20 = binomLikelihoodSmooth(pos, pos, m, n, h=20)

\dontrun{plot(x=pos, y=m/n)
points(x=pos, y=p_10, col="green")
lines(x=pos, y=p_10, col="green")
points(x=pos, y=p_20, col="red")
lines(x=pos, y=p_20, col="red")}
}

\keyword{~smooth}
