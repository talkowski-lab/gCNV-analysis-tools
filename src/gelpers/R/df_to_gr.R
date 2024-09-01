GRANGES_MCOLS_BL <- c(
    "seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular",
    "start", "end", "width", "element"
)

GRANGES_REQ_COLS <- c("chr", "start", "end")


#' Create a \code{GRanges} object from a \code{data.frame}
#'
#' Create a \code{GRanges} object from a \code{data.frame} containing genomic
#' ranges.
#'
#' @details
#' The columns 'seqnames', 'ranges', 'strand', 'seqlevels', 'seqlengths',
#' 'isCircular', 'start', 'end', 'width', and 'element' are not allowed names
#' for a \code{GRanges} metadata column so columns with these names will be
#' removed from \code{x} with a warning before assigning them to metadata
#' columns.
#'
#' @name df_to_gr
#' @export
#' @param x A \code{data.frame} with at least the columns \code{chr},
#'   \code{start}, and \code{end}. A \code{strand} column will be used for
#'   strand if available and \code{cnv} is \code{FALSE}. All columns that are
#'   not \code{chr}, \code{start}, \code{end}, or otherwise prohibited will
#'   become metadata columns.
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
    gr <- .df_to_gr(x, cnv)
    mcol_names <- .df_mcol_names(x)
    if (length(mcol_names) > 0) {
        GenomicRanges::mcols(gr) <- x[, mcol_names, drop = FALSE]
    }

    gr
}

#' @rdname df_to_gr
#' @export
#' @method df_to_gr data.table
df_to_gr.data.table <- function(x, cnv = FALSE, ...) {
    gr <- .df_to_gr(x, cnv)
    mcol_names <- .df_mcol_names(x)
    if (length(mcol_names) > 0) {
        GenomicRanges::mcols(gr) <- x[, .SD, .SDcols = mcol_names]
    }

    gr
}

.df_to_gr <- function(x, cnv) {
    assert(is_flag(cnv))
    assert(!is.na(cnv))
    assert(df_has_columns(x, .cols = GRANGES_REQ_COLS))
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

    gr
}

.df_mcol_names <- function(x) {
    other_cols <- colnames(x)[!colnames(x) %in% GRANGES_REQ_COLS]
    if (length(other_cols) > 0) {
        bad_cols <- other_cols[other_cols %in% GRANGES_MCOLS_BL]
        if (length(bad_cols) > 0) {
            warning(
                paste0(
                    "removing prohibited GRanges metadata columns: ",
                    paste0("'", bad_cols, "'", collapse = ", ")
                ),
                call. = FALSE
            )
            other_cols <- other_cols[!other_cols %in% GRANGES_MCOLS_BL]
        }
    }

    other_cols
}
