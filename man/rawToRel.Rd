\name{rawToRel}
\alias{rawToRel}
\alias{rawToRel,BSraw-method}

\title{Converts a \code{BSraw} object to a \code{BSrel} object}
\description{Determines the \code{methLevel} matrix via: \code{methReads(object) / totalReads(object)}.}
\usage{
rawToRel(object)
}

\arguments{
  \item{object}{A \code{BSraw}.}
}
\value{A \code{BSrel}.}

\author{
Katja Hebestreit
}

\seealso{
\code{\link{BSraw-class}}
\code{\link{BSrel-class}}
}

\examples{
data(rrbs)
rrbs.rel <- rawToRel(rrbs)
}
