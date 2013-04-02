\name{covBoxplots}
\alias{covBoxplots}
\alias{covBoxplots,BSraw-method}

\title{
Creates boxplots of coverages per sample
}
\description{
A boxplot per sample is plotted for the coverages of CpG-sites. 
It is constrained to CpG-sites which are covered in the respective sample (coverage \code{!= 0} and not \code{NA}).
}
\usage{R
covBoxplots(object, ...)
}

\arguments{
  \item{object}{
A \code{BSraw}.
}
  \item{\dots}{
Other graphical parameters passed to the \code{boxplot} function.
}
}

\author{
Katja Hebestreit
}

\seealso{\code{boxplot}}

\examples{
data(rrbs)
covBoxplots(rrbs)
}
