#' Remove prefixes from sequence names of a \code{GRanges} object
#'
#' @name rmPrefix
#' @export
#' @param x \code{GRanges} object
#' @param prefix \code{character()} The prefixes to prepend. Must be either a
#'   single string or a character vector with length equal to the length of
#'   \code{x}. \code{NA} prefixes are not allowed.
#' @returns \code{GRanges} object with the modified sequence names. Only the
#'   ranges, strands, and metadata columns will be retained from \code{x}. If
#'   removing a prefix would result in an empty string, an error will be
#'   signalled.
#'
#' @seealso \code{\link{addPrefix}}
setGeneric("rmPrefix",
           function(x, prefix) standardGeneric("rmPrefix"),
           signature = "x")

#' @rdname rmPrefix
setMethod("rmPrefix", c("GRanges"), function(x, prefix) {
    assert(is.character(prefix))
    assert(length(prefix) == 1 || length(prefix) == length(x),
           msg = "length of `prefix` must be 1 or equal to length of `x`")
    assert(!any(is.na(prefix)), msg = "no prefix can be `NA`")

    s <- as.character(GenomicRanges::seqnames(x))
    if (length(prefix) == 1) {
        prefix <- rep.int(prefix, length(s))
    }
    m <- startsWith(s, prefix)
    if (any(m)) {
        iden_prefix_idx <- which(m & nchar(s) == nchar(prefix))
        if (length(iden_prefix_idx) > 0) {
            stop(sprintf("prefix and sequence name at index %d are identical",
                         iden_prefix_idx[[1]]),
                 call. = FALSE)
        }
        s[m] <- substr(s[m], nchar(prefix[m]) + 1, nchar(s[m]))
    }

    out <- GenomicRanges::GRanges(s,
                                  ranges = GenomicRanges::ranges(x),
                                  strand = GenomicRanges::strand(x))
    GenomicRanges::mcols(out) <- GenomicRanges::mcols(x)

    out
})
