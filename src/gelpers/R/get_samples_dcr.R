#' Retrieve a dCR region for samples in a gCNV batch
#'
#' @description
#' Retrieve the dCR matrices for samples in a single gCNV batch, subset to the
#' ranges overlapping a given genomic region.
#'
#' @export
#' @inheritParams get_batch_dcr
#' @param samples A character vector of sample IDs. Must contain at least one
#'   ID and all samples must be from the same batch.
#' @param include_bg If \code{TRUE}, all of the other samples not among
#'   \code{samples} will be included in the returned \code{data.table}.
#'   Otherwise, only the requested samples are returned.
#' @returns A \code{data.table} with the portion of the dCR matrix overlapping
#'   the query region. The first three columns will be \code{'chr'},
#'   \code{'start'}, and \code{'end'}, corresponding to the dCR ranges. The
#'   next \var{N} columns (where \var{N} is the number samples requested) will
#'   be the dCR columns for the requested samples. The order of these columns is
#'   not guaranteed. The remaining columns will be the remaining samples in the
#'   batch, if requested. If any sample ID in \code{samples} is not found in
#'   the dCR matrix, this function will signal an error.
#'
#' @seealso get_batch_dcr
get_samples_dcr <- function(x,
                            paths,
                            samples,
                            include_bg = FALSE,
                            bound_dcr = TRUE) {
    UseMethod("get_samples_dcr")
}

#' @export
get_samples_dcr.gregion <- function(x,
                                    paths,
                                    samples,
                                    include_bg = FALSE,
                                    bound_dcr = TRUE) {
    assert(is_flag(include_bg))
    assert(is.character(samples))
    assert(length(samples) > 0)

    dcr <- get_batch_dcr(x, paths, bound_dcr)
    if (nrow(dcr) == 0) {
        stop("no ranges in the dCR overlapped the requested region")
    }

    samples <- unique(samples)
    has_all_samples <- df_has_columns(dcr, .cols = samples)
    if (!has_all_samples) {
        stop(paste0(
            "dCR is missing samples:\n",
            paste0(attr(has_all_samples, "missing_cols"), collapse = ", ")
        ))
    }

    cols <- c("chr", "start", "end", samples)
    if (include_bg) {
        cols <- c(cols, colnames(dcr)[!colnames(dcr) %in% cols])
    }

    dcr[, cols, with = FALSE]
}
