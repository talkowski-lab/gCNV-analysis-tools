#' Retrieve a dCR region for a gCNV batch
#'
#' @description
#' Retrieve all the dCR matrices for a single gCNV batch, subset to the ranges
#' overlapping a given genomic region.
#'
#' @details
#' It is critical that all of the dCR files passed to this function are actually
#' for the same batch because the ranges extracted from each matrix are
#' assumed to be exactly the same and in the same order (which is the case for
#' matrices from the same batch). This function will only check that the
#' number, not the order or coordinates, of the extracted regions from each
#' sub-matrix is the same.
#'
#' @export
#' @param x A \code{gregion} object for the target region.
#' @param paths A character vector of paths. The length must be at least 1.
#' @param bound_dcr If \code{TRUE}, all dCR values less than 0 will be set to
#'   0 and all dCR values greater than 5 will be set to 5. Otherwise the dCR
#'   matrix will be left as is.
#' @returns A \code{data.table} with the portion of the dCR matrix overlapping
#'   the query region. The first three columns will be \code{'chr'},
#'   \code{'start'}, and \code{'end'} corresponding to the dCR regions. The
#'   remaining columns will be the samples and their dCR values.
get_batch_dcr <- function(x, paths, bound_dcr = TRUE) {
    UseMethod("get_batch_dcr")
}

#' @export
get_batch_dcr.gregion <- function(x, paths, bound_dcr = TRUE) {
    assert(is.character(paths))
    assert(length(paths) > 0)
    assert(is_flag(bound_dcr))

    dcrs <- vector(mode = "list", length = length(paths))
    for (i in seq_along(paths)) {
        dcr <- .tabix(paths[[i]], x$chr, x$start, x$end)
        if (bound_dcr) {
            dcr <- .bound_dcr(dcr)
        }

        if (i == 1) {
            dcrs[[i]] <- dcr
        } else {
            dcr[, `:=`(chr = NULL, start = NULL, end = NULL)]
            dcrs[[i]] <- dcr
        }
    }
    row_counts <- vapply(dcrs, nrow, integer(1))
    if (max(row_counts) != min(row_counts)) {
        stop(
            "dCR matrices from a batch have different numbers of rows",
            call. = FALSE
        )
    }

    dcr <- do.call(cbind, dcrs)
    chr <- NULL
    start <- NULL
    data.table::setkey(dcr, chr, start)

    dcr
}

.tabix <- function(path, chr, start, end) {
    tabix_con <- Rsamtools::TabixFile(path)
    open(tabix_con)
    on.exit(close(tabix_con), add = TRUE, after = FALSE)

    param <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
    header <- Rsamtools::headerTabix(tabix_con)$header |>
        strsplit(split = "\t", fixed = TRUE) |>
        unlist()
    records <- Rsamtools::scanTabix(tabix_con, param = param)[[1]]

    out <- data.table::fread(
        text = records,
        header = FALSE,
        sep = "\t"
    )
    colnames(out) <- header
    data.table::setnames(out, "#chr", "chr")

    out
}

.bound_dcr <- function(x, lower = 0, upper = 5) {
    fn <- function(column) {
        column[column < lower] <- lower
        column[column > upper] <- upper
        column
    }

    x[, names(.SD) := lapply(.SD, fn), .SDcols = !c("chr", "start", "end")]

    x
}
