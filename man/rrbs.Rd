\name{rrbs}
\alias{rrbs}
\docType{data}
\title{
RRBS data of APL patient samples and controls.
}
\description{
RRBS data of the CpG sites CpG sites from genomic regions on p arms of chromosome 1 and 2
covered in at least one sample. Data was obtained from 5 APL patient
samples and 5 control samples (APL in remission). RRBS data was
preprocessed with the Bismark software version 0.5.
}
\usage{rrbs}
\format{
A \code{\link{BSraw-class}} object.
}
\source{
  Schoofs T, Rohde C, Hebestreit K, Klein HU, Goellner S, Schulze I,
  Lerdrup M, Dietrich N, Agrawal-Singh S, Witten A, Stoll M, Lengfelder
  E, Hofmann WK, Schlenke P, Buechner T, Hansen K, Berdel WE, Rosenbauer
  F, Dugas M, Mueller-Tidow C (2012). DNA methylation changes are a late
  event in acute promyelocytic leukemia and coincide with loss of
  transcription factor binding. Blood.

}
\references{
Krueger F, Andrews SR. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics. 2011;27:1571-1572.
}
\examples{
data(rrbs)
show(rrbs)
}
\keyword{datasets}
