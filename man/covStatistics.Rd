\name{covStatistics}
\alias{covStatistics}
\alias{covStatistics,BSraw-method}
\alias{covStatistics,BSrel-method}

\title{
Prints a short summary of coverage statistics per sample
}
\description{
This function produces information per samples about 1.) the covered CpG-sites 2.) the median of their coverages.
}
\usage{
covStatistics(object)
}

\arguments{
  \item{object}{
A \code{BiSeq} object.
}
}

\author{
Katja Hebestreit
}


\examples{
data(rrbs)
covStatistics(rrbs)
}
