% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GBins-class.R
\docType{class}
\name{GBins}
\alias{GBins}
\alias{.GBins}
\title{An S4 class to represent genomic bins}
\usage{
GBins(
  seqnames = NULL,
  ranges = NULL,
  strand = NULL,
  ...,
  seqinfo = NULL,
  seqlengths = NULL,
  reduce = FALSE
)
}
\arguments{
\item{seqnames}{\code{NULL}, or an \code{\link[S4Vectors]{Rle}} object,
charactor vector or factor containing the sequence names.}

\item{ranges}{\code{NULL}, or an \code{\link[IRanges]{IRanges}} object
containing the ranges.}

\item{strand}{\code{NULL}, or an \code{\link[S4Vectors]{Rle}} object,
character vector or factor containing the strand information. Only '*' is
allowed as the strand.}

\item{...}{Metatdata columns to set on the GBins object. All the metadata
columns must be vector-like object of the same length as the object to
construct. They cannot be named \code{start}, \code{end}, \code{width}, or
\code{element}. If \code{reduce} is \code{TRUE}, then all metadata columns
will ignored.}

\item{seqinfo}{\code{NULL}, or a \code{\link[GenomeInfoDb]{Seqinfo}} object,
character vector of unique sequence names, or a named numeric vector of
sequence lengths. When not \code{NULL}, \code{seqinfo} must be compatible
with the sequence names in \code{seqnames}, that is, it must have one
entry for each unique sequence name in \code{seqnames}. Note that it can
have additional entries not in \code{seqnames}.}

\item{seqlengths}{\code{NULL}, or a named integer vector containing the
length or NA of each sequence in \code{levels(seqnames)}.}

\item{reduce}{Should any overlapping ranges be reduced? If \code{FALSE},
passing overlapping ranges is an error.}
}
\description{
A \code{GBins} class is an extension of the \code{GRanges} class that
enforces sorted and disjoint ranges and no strand.
}
\section{Sorting}{

If the given ranges are not sorted, then the constructor will silently sort
them when creating the object. This means that the sorted ranges may not be
in the same order as the input ranges.
}

\section{Modification}{

\code{GBins} objects should not be modified or subset because it is easy
to create an invalid object.
}

\seealso{
\code{\link[GenomicRanges]{GRanges}}
}
