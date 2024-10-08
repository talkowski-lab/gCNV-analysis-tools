% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_trio_dcr.R
\name{get_trio_dcr}
\alias{get_trio_dcr}
\title{Retrieve a dCR region for a family trio}
\usage{
get_trio_dcr(
  x,
  samples,
  batches,
  dcr_map,
  include_bg = FALSE,
  keep_all_ranges = TRUE,
  ...
)
}
\arguments{
\item{x}{A \code{gregion} object for the target region.}

\item{samples}{\code{character(3)} Sample IDs of the trio. The order of the
samples must be offspring, father, and mother.}

\item{batches}{\code{character(3)} Batch IDs of the trio. The vector must be
parallel to \code{samples} i.e. \code{batches[[i]]} must be the batch for
\code{samples[[i]]}.}

\item{dcr_map}{\code{hashmap} Map between batch IDs and paths to the dCR
matrices of the batches.}

\item{include_bg}{\code{logical(1)} If \code{TRUE}, all of the samples in
the same batch as the offspring sample will be included in the return
value.}

\item{keep_all_ranges}{\code{logical(1)} If \code{TRUE}, the genomic
intervals in the returned dCR matrix will include any interval present in
the matrix for any sample in the trio. Otherwise, only genomic intervals
shared by all members of the trio will be included.}

\item{...}{Any additional arguments to pass to
\code{\link{get_samples_dcr}}.}
}
\value{
\code{data.frame} representing the dCR matrix. The columns will be
  "chr", "start", and "end", followed by the IDs in \code{samples} and then
  background samples, if requested.
}
\description{
Retrieve the dCR matrix for a father, mother and offspring trio, subset to
the ranges overlapping a genomic region.
}
\seealso{
\code{\link{get_samples_dcr}} and \code{\link{get_batch_dcr}}
}
