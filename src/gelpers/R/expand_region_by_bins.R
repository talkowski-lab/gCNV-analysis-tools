#' Expand a genomic region based on overlap of genomic bins
#'
#' @description
#' Find the genomic bins that a genomic region overlaps and expand the region
#' by a proportion of the overlapped bins.
#'
#' @details
#' If \code{gregion} overlaps \var{N} bins in \code{bins},
#' then \eqn{\lceil N \times pad \rceil}{ceil(N * pad)} bins will be added
#' upstream and downstream of the overlapping set and the total region covered
#' by the bins will become the new expanded region.
#'
#' @export
#' @param region A \code{gregion} object.
#' @param bins A \code{GBins} object.
#' @param pad The fraction of overlapped bins to use to expand \code{region}.
#' @returns A \code{gregion} object representing the expanded region. If
#'   \code{region} does not overlap any bins, it will be returned unmodified.
#'
#' @seealso expandRangesByBins
expand_region_by_bins <- function(region, bins, pad = 0.2) {
    assert(inherits_from(region, "gregion"))
    assert(inherits_from(bins, "GBins"))

    bin_chr <- as(GenomicRanges::seqnames(bins), "character")
    bin_start <- GenomicRanges::start(bins)
    bin_end <- GenomicRanges::end(bins)

    ol_idx <- which(
        region$chr == bin_chr &
            overlaps(bin_start, bin_end, region$start, region$end)
    )

    flankn <- ceiling(length(ol_idx) * pad)
    if (flankn == 0) {
        return(region)
    }
    min_idx <- max(1, ol_idx[[1]] - flankn)
    max_idx <- min(length(bins), ol_idx[[length(ol_idx)]] + flankn)

    # Make sure flanking bins don't spill into other sequences
    bin_ol_idx <- seq.int(min_idx, max_idx)
    chr_ol <- bin_ol_idx[bin_chr[bin_ol_idx] == region$chr]
    
    gregion(
        region$chr,
        bin_start[[chr_ol[[1]]]],
        bin_end[[chr_ol[[length(chr_ol)]]]]
    )
}
