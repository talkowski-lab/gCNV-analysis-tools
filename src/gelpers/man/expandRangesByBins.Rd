% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GBins-class.R
\name{expandRangesByBins}
\alias{expandRangesByBins}
\alias{expandRangesByBins,GRanges,GBins-method}
\title{Expand genomic ranges by fraction of genomic bins overlapped}
\usage{
expandRangesByBins(x, y, pad = 0.2)

\S4method{expandRangesByBins}{GRanges,GBins}(x, y, pad = 0.2)
}
\arguments{
\item{x}{A \code{GRanges} object.}

\item{y}{A \code{GBins} object.}

\item{pad}{Fraction of overlapped bins to add upstream and downstream.}
}
\value{
A \code{GRanges} object with expanded ranges.
}
\description{
Expand the ranges in a \code{GRanges} object by a fraction of the number
ranges in a \code{GBins} object that they overlap.
}
\details{
Each range in \code{x} will be checked for overlap of each range in \code{y}
to find the set of overlapping bins. Then each set will be expanded by
adding upstream and downstream bins, given that there are enough such bins
in \code{y} (as many bins as possible will be added if there are not
enough). The number of bins added is proportional to the size of the set
(the proportion being determined by \code{pad}). Then each range in \code{x}
will be expanded to the genomic region covered by its expanded set of
overlapping bins. Ranges in \code{x} that do not overlap a bin in \code{y}
will be kept unmodified.
}
