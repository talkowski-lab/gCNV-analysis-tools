#' An S4 class to represent genomic bins
#'
#' @description
#' A \code{GBins} class is an extension of the \code{GRanges} class that
#' enforces sorted and disjoint ranges and no strand.
#'
#' @section Sorting:
#' If the given ranges are not sorted, then the constructor will silently sort
#' them when creating the object. This means that the sorted ranges may not be
#' in the same order as the input ranges.
#'
#' @section Modification:
#' \code{GBins} objects should not be modified or subset because it is easy
#' to create an invalid object.
#'
#' @name GBins
#' @importClassesFrom GenomicRanges GRanges
#'
#' @seealso \code{\link[GenomicRanges]{GRanges}}
.GBins <- setClass("GBins", contains = "GRanges")

setValidity("GBins", function(object) {
    errs <- .valid_GBins(object)
    if (is.null(errs)) TRUE else errs
})

#' @rdname GBins
#' @param seqnames \code{NULL}, or an \code{\link[S4Vectors]{Rle}} object,
#'   charactor vector or factor containing the sequence names.
#' @param ranges \code{NULL}, or an \code{\link[IRanges]{IRanges}} object
#'   containing the ranges.
#' @param strand \code{NULL}, or an \code{\link[S4Vectors]{Rle}} object,
#'   character vector or factor containing the strand information. Only '*' is
#'   allowed as the strand.
#' @param ... Metatdata columns to set on the GBins object. All the metadata
#'   columns must be vector-like object of the same length as the object to
#'   construct. They cannot be named \code{start}, \code{end}, \code{width}, or
#'   \code{element}. If \code{reduce} is \code{TRUE}, then all metadata columns
#'   will ignored.
#' @param seqinfo \code{NULL}, or a \code{\link[GenomeInfoDb]{Seqinfo}} object,
#'   character vector of unique sequence names, or a named numeric vector of
#'   sequence lengths. When not \code{NULL}, \code{seqinfo} must be compatible
#'   with the sequence names in \code{seqnames}, that is, it must have one
#'   entry for each unique sequence name in \code{seqnames}. Note that it can
#'   have additional entries not in \code{seqnames}.
#' @param reduce Should any overlapping ranges be reduced? If \code{FALSE},
#'   passing overlapping ranges is an error.
#' @param seqlengths \code{NULL}, or a named integer vector containing the
#'   length or NA of each sequence in \code{levels(seqnames)}.
#' @export
GBins <- function(seqnames = NULL,
                  ranges = NULL,
                  strand = NULL,
                  ...,
                  seqinfo = NULL,
                  seqlengths = NULL,
                  reduce = FALSE) {
    assert(is_flag(reduce))
    assert(!is.na(reduce))

    gr <- GenomicRanges::GRanges(
        seqnames, ranges, strand, ..., seqinfo, seqlengths
    )

    if (reduce) {
        gr <- GenomicRanges::reduce(gr, ignore.strand = TRUE)
    }

    gr <- GenomicRanges::sort(gr)

    .GBins(gr)
}

#' Map genomic ranges to genomic bins
#'
#' @description
#' Convert a set of genomic ranges to coordinates defined by spans of genomic
#' bins.
#'
#' @details
#' \code{toBinSpace} converts each genomic bin to its relative position (first
#' bin is 1, second bin is 2, and so on) and then maps each genomic range to
#' its span of bins. Thus, if a genomic range overlapped bins 7 through 13, its
#' coordinates would become 7 to 13. Strand is ignored when determining
#' overlap, but sequence name is respected.
#'
#' @name toBinSpace
#' @export
#' @param x A \code{GRanges} object.
#' @param y A \code{GBins} object.
#' @returns A \code{GRanges} object with mapped coordinates.
setGeneric("toBinSpace",
    function(x, y) standardGeneric("toBinSpace"),
    signature = c("x", "y")
)

#' @rdname toBinSpace
setMethod("toBinSpace", c("GRanges", "GBins"), function(x, y) {
    if (length(x) == 0) {
        return(x)
    }

    if (length(y) == 0) {
        return(x[0])
    }

    ol <- GenomicRanges::findOverlaps(x, y, ignore.strand = TRUE)
    hits <- tapply(S4Vectors::subjectHits(ol), S4Vectors::queryHits(ol), range) |>
        S4Vectors::zipdown()
    idx <- as.integer(names(hits))
    if (length(idx) < length(x)) {
        warning("some ranges in `x` don't have any overlaps in `y`")
    }
    out <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(x[idx]),
        ranges = IRanges::IRanges(
            S4Vectors::first(hits), S4Vectors::second(hits)
        ),
        strand = GenomicRanges::strand(x[idx])
    )
    GenomicRanges::mcols(out) <- GenomicRanges::mcols(x[idx])

    out
})

#' Expand genomic ranges by fraction of genomic bins overlapped
#'
#' @description
#' Expand the ranges in a \code{GRanges} object by a fraction of the number
#' ranges in a \code{GBins} object that they overlap.
#'
#' @details
#' Each range in \code{x} will be checked for overlap of each range in \code{y}
#' to find the set of overlapping bins. Then each set will be expanded by
#' adding upstream and downstream bins, given that there are enough such bins
#' in \code{y} (as many bins as possible will be added if there are not
#' enough). The number of bins added is proportional to the size of the set
#' (the proportion being determined by \code{pad}). Then each range in \code{x}
#' will be expanded to the genomic region covered by its expanded set of
#' overlapping bins. Ranges in \code{x} that do not overlap a bin in \code{y}
#' will be kept unmodified.
#'
#' @name expandRangesByBins
#' @export
#' @param x A \code{GRanges} object.
#' @param y A \code{GBins} object.
#' @param pad Fraction of overlapped bins to add upstream and downstream.
#' @returns A \code{GRanges} object with expanded ranges.
setGeneric("expandRangesByBins",
    function(x, y, pad = 0.2) standardGeneric("expandRangesByBins"),
    signature = c("x", "y")
)

#' @rdname expandRangesByBins
setMethod("expandRangesByBins", c("GRanges", "GBins"), function(x, y, pad = 0.2) {
    assert(is.numeric(pad))
    assert(length(pad) == 1)
    assert(pad >= 0 && pad <= 1)
    ol <- suppressWarnings(GenomicRanges::findOverlaps(x, y, ignore.strand = TRUE))
    if (length(ol) == 0) {
        return(x) 
    }

    hits <- tapply(S4Vectors::subjectHits(ol), S4Vectors::queryHits(ol), range) |>
        S4Vectors::zipdown()
    idx <- as.integer(names(hits))
    ol_bins <- data.frame(
        query = idx,
        first = S4Vectors::first(hits),
        last = S4Vectors::second(hits)
    )
    ol_bins$nbins <- ol_bins$last - ol_bins$first + 1
    ol_bins$pad <- ceiling(ol_bins$nbins * pad)
    ol_bins$first_expanded <- pmax(1, ol_bins$first - ol_bins$pad)
    ol_bins$last_expanded <- pmin(length(y), ol_bins$last + ol_bins$pad)

    queries <- rep(ol_bins$query, ol_bins$last_expanded - ol_bins$first_expanded + 1) 
    subjects <- unlist(mapply(seq.int, ol_bins$first_expanded, ol_bins$last_expanded))

    m <- as.vector(GenomicRanges::seqnames(x)[queries], mode = "character")  ==
        as.vector(GenomicRanges::seqnames(y)[subjects], mode = "character")
    queries <- queries[m]
    subjects <- subjects[m]

    hits <- tapply(subjects, queries, range) |> S4Vectors::zipdown() 
    idx <- as.integer(names(hits))
    new_starts <- GenomicRanges::start(x)
    new_ends <- GenomicRanges::end(x)
    new_starts[idx] <- GenomicRanges::start(y)[S4Vectors::first(hits)]
    new_ends[idx] <- GenomicRanges::end(y)[S4Vectors::second(hits)]
    GenomicRanges::start(x) <- new_starts
    GenomicRanges::end(x) <- new_ends
    
    x
})

.valid_GBins <- function(x) {
    if (length(x) == 0) {
        return(NULL)
    }

    errs <- c(
        .valid_GBins_nostrand(x),
        .valid_GBins_disjoint(x),
        .valid_GBins_sorted(x)
    )

    errs
}

.valid_GBins_nostrand <- function(x) {
    s <- S4Vectors::unique(GenomicRanges::strand(x))
    if (length(s) == 1 && s == "*") {
        return(NULL)
    }

    "GBins ranges are unstranded (strand can only be '*')"
}

.valid_GBins_disjoint <- function(x) {
    if (GenomicRanges::isDisjoint(x)) {
        return(NULL)
    }

    "GBins ranges are disjoint"
}

.valid_GBins_sorted <- function(x) {
    if (!GenomicRanges::is.unsorted(x)) {
        return(NULL)
    }

    "GBins ranges are sorted"
}
