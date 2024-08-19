#' Retrieve a dCR region for a family trio
#'
#' @description
#' Retrieve the dCR matrix for a father, mother and offspring trio, subset to
#' the ranges overlapping a genomic region.
#'
#' @export
#' @param x A \code{gregion} object for the target region.
#' @param samples \code{character(3)} Sample IDs of the trio. The order of the
#'   samples must be offspring, father, and mother.
#' @param batches \code{character(3)} Batch IDs of the trio. The vector must be
#'   parallel to \code{samples} i.e. \code{batches[[i]]} must be the batch for
#'   \code{samples[[i]]}.
#' @param dcr_map \code{hashmap} Map between batch IDs and paths to the dCR
#'   matrices of the batches.
#' @param include_bg \code{logical(1)} If \code{TRUE}, all of the samples in
#'   the same batch as the offspring sample will be included in the return
#'   value.
#' @param keep_all_ranges \code{logical(1)} If \code{TRUE}, the genomic
#'   intervals in the returned dCR matrix will include any interval present in
#'   the matrix for any sample in the trio. Otherwise, only genomic intervals
#'   shared by all members of the trio will be included.
#' @param ... Any additional arguments to pass to
#'   \code{\link{get_samples_dcr}}.
#' @return \code{data.frame} representing the dCR matrix. The columns will be
#'   "chr", "start", and "end", followed by the IDs in \code{samples} and then
#'   background samples, if requested.
#'
#' @seealso \code{\link{get_samples_dcr}} and \code{\link{get_batch_dcr}}
get_trio_dcr <- function(x,
                         samples,
                         batches,
                         dcr_map,
                         include_bg = FALSE,
                         keep_all_ranges = TRUE,
                         ...) {
    UseMethod("get_trio_dcr")
}

#' @export
get_trio_dcr.gregion <- function(x,
                                 samples,
                                 batches,
                                 dcr_map,
                                 include_bg = FALSE,
                                 keep_all_ranges = TRUE,
                                 ...) {
    assert(is.character(samples))
    assert(is.character(batches))
    assert(length(samples) == 3L)
    assert(length(batches) == length(samples))
    assert(inherits_from(dcr_map, "hashtab"))
    assert(is_flag(include_bg))
    assert(!is.na(include_bg))
    assert(is_flag(keep_all_ranges))
    assert(!is.na(keep_all_ranges))

    coord_cols <- c("chr", "start", "end")
    dcr_groups <- split(samples, batches)
    dcrs <- vector(mode = "list", length = length(dcr_groups))
    for (i in seq_along(dcr_groups)) {
        batch_i <- names(dcr_groups)[[i]]
        samples_i <- dcr_groups[[i]]

        paths <- utils::gethash(dcr_map, batch_i)
        if (is.null(paths)) {
            stop(sprintf("batch '%s' is not a key in the dCR paths map", batch_i),
                 call. = FALSE)
        }
        dcr <- get_samples_dcr(x,
                               paths,
                               samples_i,
                               include_bg = include_bg && (samples[[1]] %in% samples_i),
                               ...)

        others <- setdiff(samples, samples_i)
        if (length(others) > 0) {
            dcr[others] <- NULL
        }
        dcrs[[i]] <- dcr
    }

    merged_dcr <- Reduce(\(x, y) merge(x, y, by = coord_cols, all = keep_all_ranges),
                         dcrs)

    bg <- setdiff(colnames(merged_dcr), c(coord_cols, samples))

    merged_dcr[c(coord_cols, samples, bg)]
}
