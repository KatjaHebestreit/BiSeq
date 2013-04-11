\name{BSrel-class}
\Rdversion{1.1}
\docType{class}
\alias{BSrel-class}
\alias{BSrel}
\alias{BSrel,matrix,GRanges-method}
\alias{combine,BSrel,BSrel-method}
\alias{methLevel}
\alias{methLevel<-}
\alias{methLevel,BSrel-method}
\alias{methLevel<-,BSrel,matrix-method}

\title{
  Class to contain Bisulfite Sequencing (BiSeq) Data
}
\description{
  The \code{BSrel} class is derived from
  \code{SummarizedExperiment} and contains a \code{SimpleList} of one
  matrix named \code{methLevel} as \code{assays}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{BSrel(exptData = SimpleList(),
    rowData,
    colData = DataFrame(row.names=colnames(methLevel)),
    methLevel,
    ...)}.

  However, one will most likely create a \code{BSraw} object when use
  \code{\link{readBismark}} to load data.
}
\section{Slots}{
  \describe{
    \item{\code{exptData}:}{An optional \code{SimpleList} of arbitrary
      content describing the overall experiment.}
    \item{\code{rowData}:}{Object of class \code{"GRanges"}
      containing the genome positions of CpG-sites covered by bisulfite
      sequencing.}
    \item{\code{colData}:}{Object of class \code{"DataFrame"}
      containing information on variable values of the samples.}
    \item{\code{assays}:}{Object of class \code{SimpleList} of a
      matrix, named \code{methLevel} containing the methylation levels
      (between 0 and 1) per CpG site. The rows represent the CpG sites
      in \code{rowData} and the columns represent the samples in \code{colData}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SummarizedExperiment}"}, directly.
}
\section{Methods}{
  \describe{

    \item{methLevel}{\code{signature(x = "BSrel")}: Gets the \code{methLevel} slot.}
    \item{methLevel<-}{\code{signature(x = "BSrel", value = "matrix")}:
      Sets the \code{methLevel} slot.}

     \item{combine}{\code{signature(x = "BSrel", y = "BSrel")}: Combines two
      \code{BSrel} objects.}

  }
}

\author{
Katja Hebestreit
}

\seealso{
\code{SummarizedExperiment}, \code{\link{BSraw-class}}, \code{\link{readBismark}}
}
\examples{
showClass("BSrel")

## How to create a BSrel object by hand:
exptData <- SimpleList(Sequencer = "Sequencer", Year = "2013")
rowData <- GRanges(seqnames = "chr1",
                  ranges = IRanges(start = c(1,2,3), end = c(1,2,3)))
colData <- DataFrame(group = c("cancer", "control"),
                     row.names = c("sample_1", "sample_2"))
methLevel <- matrix(c(rep(0.5, 3), rep(1, 3)), ncol = 2)
BSrel(exptData = exptData,
      rowData = rowData,
      colData = colData,
      methLevel = methLevel)

# Or get a BSrel object out of a BSraw object:
data(rrbs)
rrbs.rel <- rawToRel(rrbs)
}
\keyword{classes}
