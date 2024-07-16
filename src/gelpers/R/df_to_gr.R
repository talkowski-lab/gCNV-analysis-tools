#' Make a \code{GRanges} object from a \code{data.frame}
#'
#' Make a \code{GRanges} object from a \code{data.frame}
#'
#' @export
#' @param x A \code{data.frame} with at least the columns "chr", "start", and
#'   "end". A "strand" column will be used for strand if available. All other
#'   columns will become the metadata columns.
#' @returns A \code{GRanges} object.
df_to_gr <- function(x) {
    UseMethod("df_to_gr")
}

#' @export
#' @method df_to_gr data.frame
df_to_gr.data.frame <- function(x) {
    gr <- .df_to_gr(x)
    GenomicRanges::mcols(gr$gr) <- x[, gr$other_cols]

    gr$gr
}

#' @export
#' @method df_to_gr data.table
df_to_gr.data.table <- function(x) {
    gr <- .df_to_gr(x)
    GenomicRanges::mcols(gr$gr) <- x[, gr$other_cols, with = FALSE]

    gr$gr
}

.df_to_gr <- function(x) {
    chr <- NULL
    start <- NULL
    end <- NULL
    assert(df_has_columns(x, chr, start, end))
    s <- if ("strand" %in% colnames(x)) x$strand else "*"
    gr <- GenomicRanges::GRanges(
        seqnames = x$chr,
        ranges = IRanges::IRanges(x$start, x$end),
        strand = s
    )
    cols <- colnames(x)[!colnames(x) %in% c("chr", "start", "end", "strand")]
    
    list(gr = gr, other_cols = cols)
}
