\name{readBismark}
\alias{readBismark}
\alias{readBismark,character,DataFrame-method}
\alias{readBismark,character,character-method}
\alias{readBismark,character,data.frame-method}


\title{Reads cytosine methylation stati determined by Bismark}

\description{Bismark is a bisulfite read mapper and methylation caller. This method reads Bismark's
output files and returns a \code{BSraw} object.}

\usage{readBismark(files, colData)}

\arguments{
  \item{files}{A \code{character} pointing to \code{cov} files created by
  Bismark's \code{methylation_extractor} and \code{bismark2bedGraph}; see Details.
  This can be a compressed file (see \link{file}).}
  \item{colData}{Samples' names plus additional sample information
    as \code{character}, \code{data.frame} or \code{DataFrame}.}
}

\details{Input files are created with Bismark as follows (from the
  command line):

  \code{bismark_methylation_extractor -s --comprehensive test_sample.sam}

  \code{bismark2bedGraph -o CpG_context_test_sample.bedGraph
    CpG_context_test_sample.txt}

  This will output two files, a \code{.bedGraph} and a \code{.cov} file. 
  We will import the \code{CpG_context_test_sample.cov} using \code{readBismark}.  

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
file <- system.file("extdata", "CpG_context_test_sample.cov", package = "BiSeq")
rrbs <- readBismark(file,
                    colData= DataFrame(row.names="sample_1"))
}
