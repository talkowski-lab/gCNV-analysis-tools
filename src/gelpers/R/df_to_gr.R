#' Create a \code{GRanges} object from a \code{data.frame}
#'
#' Create a \code{GRanges} object from a \code{data.frame} containing genomic
#' ranges.
#'
#' @name df_to_gr
#' @export
#' @param x A \code{data.frame} with at least the columns \code{chr},
#'   \code{start}, and \code{end}. A \code{strand} column will be used for
#'   strand if available and \code{cnv} is \code{FALSE}. All columns that are
#'   not \code{chr}, \code{start}, or \code{end} will become metadata columns.
#' @param cnv If \code{TRUE}, the genomic ranges will be treated as CNVs. In
#'   this case, there should be an \code{svtype} column in \code{x} indicating
#'   the CNV type (either 'DUP' or 'DEL'). DUPs will be given a strand of '+'
#'   and DELs will given a strand of '-'.
#' @param ... Other arguments passed to methods.
#' @returns A \code{GRanges} object.
df_to_gr <- function(x, ...) {
    UseMethod("df_to_gr")
}

#' @rdname df_to_gr
#' @export
#' @method df_to_gr data.frame
df_to_gr.data.frame <- function(x, cnv = FALSE, ...) {
    req_cols <- c("chr", "start", "end")
    assert(df_has_columns(x, .cols = req_cols))
    if (cnv) {
        if (!"svtype" %in% colnames(x)) {
            stop("`svtype` is not in `x`, but `cnv == TRUE`", call. = FALSE)
        }
        if (any(!x$svtype %in% c("DUP", "DEL"))) {
            stop("only 'DUP' and 'DEL' are allowed in `svtype`", call. = FALSE)
        }
        s <- ifelse(x$svtype == "DUP", "+", "-")
    } else {
        s <- if ("strand" %in% colnames(x)) x$strand else "*"
    }
    gr <- GenomicRanges::GRanges(
        seqnames = x$chr,
        ranges = IRanges::IRanges(x$start, x$end),
        strand = s
    )

    
    other_cols <- colnames(x)[!colnames(x) %in% req_cols]
    if (length(other_cols) > 0) {
        GenomicRanges::mcols(gr) <- x[, other_cols]
    }

    gr
}
