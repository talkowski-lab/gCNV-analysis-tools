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
#' @param squeeze If \code{TRUE}, all dCR values less than 0 will be set to
#'   0 and all dCR values greater than 5 will be set to 5. Otherwise the dCR
#'   matrix will be left as is.
#' @returns A \code{data.frame} with the portion of the dCR matrix overlapping
#'   the query region. The first three columns will be \code{'chr'},
#'   \code{'start'}, and \code{'end'} corresponding to the dCR regions. The
#'   remaining columns will be the samples and their dCR values.
get_batch_dcr <- function(x, paths, squeeze = TRUE) {
    UseMethod("get_batch_dcr")
}

#' @export
get_batch_dcr.gregion <- function(x, paths, squeeze = TRUE) {
    assert(is.character(paths))
    assert(length(paths) > 0)
    assert(!any(is.na(paths)), msg = "`NA` paths are not allowed")
    assert(is_flag(squeeze))
    assert(!is.na(squeeze))

    dcrs <- vector(mode = "list", length = length(paths))
    rows <- 0L
    for (i in seq_along(paths)) {
        dcr <- .tabix(paths[[i]], x$chr, x$start, x$end)

        if (squeeze && nrow(dcr) > 0) {
            dcr <- .squeeze_dcr(dcr)
        }

        if (i == 1) {
            dcrs[[i]] <- dcr
            rows <- nrow(dcr)
        } else {
            if (nrow(dcr) != rows) {
                stop(
                    "dCR matrices for a batch have different row counts",
                    call.= FALSE
                )
            }
            dcr[c("chr", "start", "end")] <- list(NULL, NULL, NULL)
            dcrs[[i]] <- dcr
        }
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
    meta <- Rsamtools::headerTabix(tabix_con)
    if (length(meta[["header"]]) == 0) {
        stop(.dcr_parse_error(path, "missing header"))
    }

    if (!chr %in% meta[["seqnames"]]) {
        stop(.dcr_parse_error(path, sprintf("no sequence '%s'", chr)))
    }

    if (length(meta[["comment"]]) > 0 && startsWith(meta[["header"]], meta[["comment"]])) {
        header <- substr(meta[["header"]], 1 + nchar(meta[["comment"]]), nchar(meta[["header"]]))
    } else {
        header <- meta[["header"]]
    }

    header <- strsplit(header, split = "\t", fixed = TRUE)[[1]]
    # header must have columns for chr, start, and end
    if (length(header) < 3 || header[[1]] != "chr" || header[[2]] != "start" || header[[3]] != "end") {
        stop(
            .dcr_parse_error(
                path,
                "first three columns must be 'chr', 'start', and 'end'"
            )
        )
    }

    if (anyDuplicated(header)) {
        stop(.dcr_parse_error(path, "duplicate samples"))
    }

    records <- Rsamtools::scanTabix(tabix_con, param = param)[[1]]
    coltypes <- c(
        list(character(), integer(), integer()),
        lapply(seq_len(length(header) - 3), \(x) double())
    )
    if (length(records) == 0) {
        names(coltypes) <- header
        return(as.data.frame(coltypes, check.names = FALSE))
    }

    out <- scan(
        text = records,
        what = coltypes,
        nmax = length(records),
        sep = "\t",
        quiet = TRUE
    )
    names(out) <- header

    as.data.frame(out, check.names = FALSE)
}

.squeeze_dcr <- function(x, lower = 0, upper = 5) {
    if (ncol(x) <= 3) {
        return(x)
    }

    for (i in seq.int(4, ncol(x), 1)) {
        x[[i]][x[[i]] < lower] <- lower
        x[[i]][x[[i]] > upper] <- upper
    }

    x
}

.dcr_parse_error <- function(path, error, call = NULL) {
    msg <- paste0(
        sprintf("parsing dCR in '%s' failed", path),
        "\n- ",
        error
    )
    errorCondition(msg, class = "dcr_parse_error", call = call)
}
