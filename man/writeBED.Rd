\name{writeBED}
\alias{writeBED,BSraw,character,character-method}
\alias{writeBED,BSraw,missing,character-method}
\alias{writeBED,BSraw,missing,missing-method}
\alias{writeBED,BSraw,character,missing-method}
\alias{writeBED,BSrel,character,character-method}
\alias{writeBED,BSrel,missing,character-method}
\alias{writeBED,BSrel,missing,missing-method}
\alias{writeBED,BSrel,character,missing-method}

\alias{writeBED}

\title{Writes \code{BSraw} and \code{BSrel} data to a bed file suitable for the IGV}

\description{The created bed files contains an entry for each
CpG site. Strand information, relative methylation and
absolute number of reads covering the CpG sites are stored.
The relative methylation is indicated by colors:
green via black to red for unmethylated to methylated.}

\usage{writeBED(object, name, file)}

\arguments{
  \item{object}{A \code{BSraw} or \code{BSrel} object.}
  \item{name}{Track names (sample names) written to the bed file's header.}
  \item{file}{Character vector with names of the bed file.}
}
\details{The written bed file contains the following
extra information:
\enumerate{
\item score: the relative methylation of the CpG site
\item name: the coverage of the CpG site
\item itemRgb: a color value visualizing the methylation score
}
A separate bed file is created for each sample in the
given object. The lengths of the arguments
\code{name} and \code{file} should equal the number of
samples.
}

\value{Nothing. Bed files are written.}

\author{Hans-Ulrich Klein}

\seealso{
\code{\link{readBismark}}
}

\examples{
data(rrbs)
s1 <- rrbs[,1]
out <- tempfile(, fileext = ".bed")
writeBED(s1, name = colnames(s1), file = out)

}
