% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simGenotypeInternal.R
\name{setGeno}
\alias{setGeno}
\title{setGeno}
\usage{
setGeno(infoSNV, genotype)
}
\arguments{
\item{infoSNV}{a \code{list} of snv = \code{data.frame} with
column : chr, pos, ref, alt, gtype
and
listSubset = index of the snv keep from the snv input}

\item{genotype}{\code{data.frame}
one column format 0|0, 0|1, 1|0 or 1|1
one row per snv same row then snv}

\item{fileList}{a \code{list} of \code{GRanges}, the segments from multiple
files.}
}
\value{
a \code{list} of snv = \code{data.frame} with
column : chr, pos, ref, alt, gtype
and
listSubset = index of the snv keep from the snv input
}
\description{
Generate a subset of snv with genotype
}
\examples{

# TODO

}
\author{
Pascal Belleau, Astrid Deschenes and
Alexander Krasnitz
}
\keyword{internal}
