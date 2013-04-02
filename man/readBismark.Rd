\name{readBismark}
\alias{readBismark}
\alias{readBismark,character,character,DataFrame-method}
\alias{readBismark,character,character,character-method}
\alias{readBismark,character,character,data.frame-method}


\title{Reads cytosine methylation stati determined by Bismark}

\description{Bismark is a tool to map bisulfite converted sequence reads
and to determine cytosine methylation states. This method reads Bismark's
output files and returns a \code{BSraw} object.}

\usage{readBismark(otFiles, obFiles, colData)}

\arguments{
  \item{otFiles}{A \code{character} pointing to Bismark's
    output file for the plus strand (original top).}
  \item{obFiles}{A \code{character} pointing to Bismark's
    output file for the minus strand (original bottom).}
  \item{colData}{Samples' names plus additional sample information
    as \code{character}, \code{data.frame} or \code{DataFrame}.}
}

\details{The arguments \code{otFiles} and \code{obFiles} point to
Bismark's output files for the positive (original top) and the negative
(original bottom) strand respectively. Both \code{character} vectors must be of equal
length and each vector element corresponds to one sample.

The \code{colData} argument should specify the sample names as
\code{character}. Alternatively, a \code{data.frame} or
\code{DataFrame} can be given. Then, the row names are used as
sample names and the data frame is passed to the final
\code{BSraw} object.}

\value{A \code{BSraw} object storing coverage and methylation information.}

\references{http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/}

\author{Hans-Ulrich Klein}

\seealso{\code{\link{BSraw-class}}}

\examples{
ot <- system.file("extdata", "bismark_OT.txt", package = "BiSeq")
ob <- system.file("extdata", "bismark_OB.txt", package = "BiSeq")
rrbs <- readBismark(otFiles = ot,
                    obFiles = ob,
                    colData= DataFrame(row.names="sample_1"))
}
