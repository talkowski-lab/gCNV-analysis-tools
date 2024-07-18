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
#' @returns A \code{data.frame} with the portion of the dCR matrix overlapping
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
        if (ncol(dcr) == 0) {
            stop("dCR is missing a header", call. = FALSE)
        }

        assert(df_has_columns(dcr, .cols = c("#chr", "start", "end")))
        colnames(dcr)[colnames(dcr) == "#chr"] <- "chr"

        dcr <- .fix_dcr_coltypes(dcr)

        if (bound_dcr) {
            dcr <- .bound_dcr(dcr)
        }

        if (i == 1) {
            dcrs[[i]] <- dcr
        } else {
            dcr[c("chr", "start", "end")] <- list(NULL, NULL, NULL)
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

    out_dcr <- do.call(cbind, dcrs)
    order_idx <- order(out_dcr$chr, out_dcr$start)

    out_dcr[order_idx, ]
}

.tabix <- function(path, chr, start, end) {
    tabix_con <- Rsamtools::TabixFile(path)
    open(tabix_con)
    on.exit(close(tabix_con), add = TRUE, after = FALSE)

    param <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
    header <- Rsamtools::headerTabix(tabix_con)$header
    if (is.null(header) || length(header) == 0) {
        header <- character()
    } else {
        header <- strsplit(header, split = "\t", fixed = TRUE) |>
            unlist()
    }
    records <- Rsamtools::scanTabix(tabix_con, param = param)[[1]] |>
        strsplit(split = "\t", fixed = TRUE)

    if (length(records) == 0) {
        if (length(header) == 0) {
            return(data.frame())
        }
        out <- lapply(header, \(x) character())
        names(out) <- header
        out <- as.data.frame(out)
        return(out)
    }

    out <- do.call("rbind", records) |>
        as.data.frame()
    if (length(header) == ncol(out)) {
        colnames(out) <- header
    }

    out
}

.bound_dcr <- function(x, lower = 0, upper = 5) {
    cols <- which(!colnames(x) %in% c("chr", "start", "end"))
    if (length(cols) == 0) {
        return(x)
    }
    x[cols] <- lapply(x[cols], .bound_dcr_col)

    x
}

.bound_dcr_col <- function(x, lower, upper) {
    x[x < lower] <- lower
    x[x > upper] <- upper
    x
}

.fix_dcr_coltypes <- function(x) {
    x$start <- as.integer(x$start)
    x$end <- as.integer(x$end)

    other_cols <- which(!colnames(x) %in% c("chr", "start", "end"))
    x[other_cols] <- lapply(x[other_cols], as.double)

    x
}
