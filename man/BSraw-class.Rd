\name{BSraw-class}
\Rdversion{1.1}
\docType{class}
\alias{BSraw-class}
\alias{BSraw}
\alias{BSraw,matrix,matrix,GRanges-method}
\alias{combine,BSraw,BSraw-method}
\alias{methReads}
\alias{methReads<-}
\alias{totalReads}
\alias{totalReads<-}
\alias{totalReads,BSraw-method}
\alias{totalReads<-,BSraw,matrix-method}
\alias{methReads,BSraw-method}
\alias{methReads<-,BSraw,matrix-method}

\title{Class to contain raw Bisulfite Sequencing (BiSeq) Data}
\description{
  The \code{BSraw} class is derived from
  \code{SummarizedExperiment} and contains a \code{SimpleList} of
  matrices named \code{methReads} and \code{totalReads} as \code{assays}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{BSraw(exptData = SimpleList(),
    rowData,
    colData = DataFrame(row.names=colnames(methReads)),
    methReads,
    totalReads,
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
    \item{\code{assays}:}{Object of class \code{SimpleList} of two
      matrices, named \code{totalReads} and \code{methReads}. The matrix
      \code{totalReads} contains the number of reads spanning a
      CpG-site. The rows represent the CpG sites in \code{rowData} and
      the columns represent the samples in \code{colData}. The matrix
      \code{methReads} contains the number of methylated reads spanning
      a CpG-site.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SummarizedExperiment}"}, directly.
}
\section{Methods}{
  \describe{
    
    \item{totalReads}{\code{signature(x = "BSraw")}: Gets the \code{totalReads} slot.}
    \item{totalReads<-}{\code{signature(x = "BSraw", value = "matrix")}:
      Sets the \code{totalReads} slot.}

    \item{methReads}{\code{signature(x = "BSraw")}: Gets the \code{methReads} slot.}
    \item{methReads<-}{\code{signature(x = "BSraw", value = "matrix")}:
      Sets the \code{methReads} slot.}

    \item{combine}{\code{signature(x = "BSraw", y = "BSraw")}: Combines two
      \code{BSraw} objects.}
  }
  }

\author{
Katja Hebestreit
}

\seealso{
  \code{SummarizedExperiment}, \code{\link{BSrel-class}}, \code{\link{readBismark}}
}
\examples{
showClass("BSraw")

## How to create a BSraw object by hand:
exptData <- SimpleList(Sequencer = "Sequencer", Year = "2013")
rowData <- GRanges(seqnames = "chr1",
                  ranges = IRanges(start = c(1,2,3), end = c(1,2,3)))
colData <- DataFrame(group = c("cancer", "control"),
                     row.names = c("sample_1", "sample_2"))
totalReads <- matrix(c(rep(10L, 3), rep(5L, 3)), ncol = 2)
methReads <- matrix(c(rep(5L, 3), rep(5L, 3)), ncol = 2)
BSraw(exptData = exptData,
      rowData = rowData,
      colData = colData,
      totalReads = totalReads,
      methReads = methReads)

## A more realistic example can be loaded:
data(rrbs)
rrbs

head(totalReads(rrbs))
head(methReads(rrbs))

}
\keyword{classes}
