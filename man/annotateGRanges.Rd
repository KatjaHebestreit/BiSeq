\name{annotateGRanges}
\alias{annotateGRanges}
\alias{annotateGRanges,GRanges,GRanges,character,character-method}
\alias{annotateGRanges,GRanges,GRanges,character,integer-method}
\alias{annotateGRanges,GRanges,GRanges,character,missing-method}

\title{Annotates a \code{GRanges} object by means of a second \code{GRanges} object}
\description{
Each genomic location of \code{object} is checked for overlapping with
genomic ranges of \code{regions}. In case of an overlapping, this
genomic location is marked as \code{TRUE}, or with the identifier of
respective the \code{regions} object (if any).
}
\usage{
annotateGRanges(object, regions, name, regionInfo)
}
\arguments{
  \item{object}{
A \code{GRanges}.
}
  \item{regions}{
A \code{GRanges}.
}
  \item{name}{
A string specifying the name of the \code{elementMetadata}
column with the overlapping information to be added to
\code{object}. Usually the type of the \code{regions} object, e.g. 'Promoter'.
}
  \item{regionInfo}{
OPTIONAL: A string or integer specifying the \code{elementMetadata}
column of \code{regions} object containing the identifiers,
e.g. entrez gene IDs of the promoters. If not specified, the genomic
locations of \code{object} are labeled as
\code{TRUE} (overlap) or \code{FALSE} (no overlap).
}
}

\details{
If multiple ranges of \code{regions} overlap with a genomic
region in \code{object}, the identifier names of the overlapping regions
are seperated by ','.
}
\value{
A \code{GRanges} object similar to \code{object} containing an
additional \code{elementMetadata} column with the overlapping information.
}


\author{
Katja Hebestreit
}

\seealso{\code{GRanges-class}}


\examples{
# load detected DMRs:
data(DMRs)

# annotate the DMRs with a GRanges object:
data(promoters)
DMRs.anno <- annotateGRanges(object = DMRs,
                             regions = promoters,
                             name = 'Promoter',
                             regionInfo = 'acc_no')
DMRs.anno
}
