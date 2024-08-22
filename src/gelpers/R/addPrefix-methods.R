#' Add prefixes to sequence names of a \code{GRanges} object
#'
#' @name addPrefix
#' @export
#' @param x \code{GRanges} object
#' @param prefix \code{character()} The prefixes to prepend. Must be either a
#'   single string or a character vector with length equal to the length of
#'   \code{x}.
#' @returns \code{GRanges} object with the modified sequence names. Only the
#'   ranges, strands, and metadata columns will be retained from \code{x}.
#'
#' @seealso \code{\link{rmPrefix}}
setGeneric("addPrefix",
           function(x, prefix) standardGeneric("addPrefix"),
           signature = "x")

#' @rdname addPrefix
setMethod("addPrefix", c("GRanges"), function(x, prefix) {
    assert(is.character(prefix))
    assert(length(prefix) == 1 || length(prefix) == length(x),
           msg = "length of `prefix` must be 1 or equal to length of `x`")
    assert(!any(is.na(prefix)), msg = "no prefix can be `NA`")

    s <- as.character(GenomicRanges::seqnames(x))
    out <- GenomicRanges::GRanges(paste0(prefix, s),
                                  ranges = GenomicRanges::ranges(x),
                                  strand = GenomicRanges::strand(x))
    GenomicRanges::mcols(out) <- GenomicRanges::mcols(x)

    out
})
