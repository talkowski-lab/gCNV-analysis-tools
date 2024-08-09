#' Read gCNV bins file
#'
#' @description
#' Read the file of genomic ranges used by gCNV to bin the genome.
#'
#' @details
#' The file should be a three-column tab-delimited file containing the
#' chromosome, start, and end of each interval, in that order. Coordinates are
#' assumed to be 1-start, closed. Any columns beyond the first three will be
#' ignored. The intervals must be disjoint.
#'
#' @export
#' @param path Path to the file.
#' @param reduce Should overlapping ranges be merged?
#' @returns A \code{GBins} object.
read_gcnv_bins <- function(path, reduce = FALSE) {
    assert(is_string(path))
    assert(is_flag(reduce))
    assert(!is.na(reduce))

    x <- utils::read.table(
        path,
        header = FALSE,
        sep = "\t",
    )

    gb <- GBins(
        seqnames = x[[1]],
        ranges = IRanges::IRanges(x[[2]], x[[3]]),
        reduce = reduce
    )

    gb
}
