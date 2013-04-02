\name{smoothVariogram}
\alias{smoothVariogram}
\alias{smoothVariogram,list,numeric,missing-method}
\alias{smoothVariogram,list,numeric,numeric-method}
\alias{smoothVariogram,matrix,numeric,missing-method}
\alias{smoothVariogram,matrix,numeric,numeric-method}

\title{
Smoothes variogram
}
\description{
Nonparametric smoothing with kernel regression estimators and adaptable
bandwidth for variogram smoothing.
}
\usage{
smoothVariogram(variogram, sill, bandwidth)
}

\arguments{
  \item{variogram}{
A \code{list} or a \code{matrix}. Usually the output of \code{makeVariogram}.
}
  \item{sill}{
A \code{numeric}. The sill (upper bound) of the variogram. See Details.
}
  \item{bandwidth}{
A numeric vector of same length as the variogram (number of
rows). Default: \code{seq(10,1000, length.out=nrow(variogram))}. See Details. 
}
}
\details{
It is necessary to smooth the variogram. Especially for greater
\code{h} the variogram tends to oscillate strongly. This is the reason
why the default bandwidth increases with increasing
\code{h}. Nevertheless, the smoothed variogram may further increase or
decrease after a horizontal part (sill). This is mostly due to the small number
of observations for high distances. To wipe out this bias it is useful to
set the smoothed variogram to a fixed value above a certain \code{h}, usually the mean value of the
horizontal part. If a smoothed value \eqn{v.sm} is greater than
\code{sill} for distance \eqn{h_{range}}, this \eqn{v.sm} and all other
smoothed values with \eqn{h > h_{range}} are set to \code{sill}. Internally, the function \code{lokerns} from package \code{lokerns}
is used for smoothing.
}
\value{The variogram matrix (or a list with the variogram matrix) with an additional column of the smoothed
  \code{v} values.
}
\author{
Katja Hebestreit
}

\seealso{
  \code{\link{makeVariogram}}, \code{lokerns}
  }

\examples{
data(vario)

# Find out the sill (this is more obvious for larger data sets):
plot(vario$variogram)

vario.sm <- smoothVariogram(vario, sill = 0.9)

plot(vario$variogram)
lines(vario.sm$variogram[,c("h", "v.sm")],
col = "red")

}