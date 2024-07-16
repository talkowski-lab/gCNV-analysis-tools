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
#' @returns A \code{GBins} object.
read_gcnv_bins <- function(path) {
    assert(is_string(path))
    x <- data.table::fread(
        path,
        header = FALSE,
        sep = "\t",
        select = 1:3,
        col.names = c("chr", "start", "end"),
        showProgress = FALSE,
        key = c("chr", "start")
    )

    gb <- GBins(
        seqnames = x$chr,
        ranges = IRanges::IRanges(x$start, x$end)
    )

    gb
}
