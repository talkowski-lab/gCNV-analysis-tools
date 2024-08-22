#' Compute coverage between two \code{GRanges} objects
#'
#' @description
#' For every range in \code{x}, compute the proportion of the range that is
#' overlapped by ranges in \code{y}. This function is intended to be used to
#' get the coverage between sets of CNVs.
#'
#' @name cnvCoverage
#' @export
#' @param x,y \code{GRanges} objects. Because the intent is to compute
#'   coverage of CNVs, both objects should be stranded according to SV type.
#'   "DUP"s should be "+" and "DEL"s should be "-".
#' @returns A \code{double} vector equal to the length of \code{x} with the
#'   coverage of each range in \code{x}.
setGeneric("cnvCoverage",
           function(x, y) standardGeneric("cnvCoverage"),
           signature = c("x", "y"))

#' @rdname cnvCoverage
setMethod("cnvCoverage", c("GRanges", "GRanges"), function(x, y) {
    xstrands <- S4Vectors::runValue(GenomicRanges::strand(x))
    if ("*" %in% xstrands) {
        warning("`x` has ranges with no strand", call. = FALSE)
    }
    ystrands <- S4Vectors::runValue(GenomicRanges::strand(y))
    if ("*" %in% ystrands) {
        warning("`y` has ranges with no strand", call. = FALSE)
    }

    out <- double(length(x))
    ol <- suppressWarnings(GenomicRanges::findOverlaps(x, y))
    if (length(ol) == 0) {
        return(out)
    }

    qh <- S4Vectors::queryHits(ol)
    sh <- S4Vectors::subjectHits(ol)
    int <- suppressWarnings(GenomicRanges::pintersect(x[qh], y[sh]))
    cov <- by(IRanges::width(int) / IRanges::width(x[qh]), qh, sum)
    out[as.integer(names(cov))] <- cov

    out
})
