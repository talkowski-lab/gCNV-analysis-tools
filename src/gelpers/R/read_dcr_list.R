#' Read list of dCR matrix paths
#'
#' @description
#' Read a list of dCR matrix paths and create a hash table from gCNV batch to
#' dCR matrix path.
#'
#' @details
#' There are two options for the input file. If the file is a one-column table
#' (when a tab is used as the separator), then it will be treated as a list of
#' file paths, one per line. The paths must be of the form
#' "path/to/\{BATCH\}_((CASE)|(COHORT)).dcr.bed.gz" where "\{BATCH\}" is the batch
#' ID. If a path is not of this form, this function will signal an error.
#'
#' If the file is a two-column table, then the first column must be the batch
#' ID and the second columm must be the path to a dCR matrix for that batch.
#' It is assumed that all rows with the same batch ID are actually from the
#' same batch. If there are any empty strings in either column, this function
#' will signal an error.
#'
#' Any other file format is an error.
#'
#' @export
#' @param path Path to the file of dCR matrix paths.
#' @returns \code{hashtab} in which the keys are batch IDs and the values are
#'   character vectors of paths to the dCR matrices for the batch.
read_dcr_list <- function(path) {
    assert(is_string(path))

    df <- utils::read.table(path, sep = "\t", header = FALSE)
    if (ncol(df) > 2) {
        stop(
            "dCR matrix paths file must be either one or two columns",
            call. = FALSE
        )
    }

    if (ncol(df) == 1) {
        paths <- df[[1]]
        batches <- sub(
            "_((CASE)|(COHORT))\\.dcr\\.bed\\.gz$",
            "",
            basename(paths)
        )
        parse_fails <- .valid_dcr_batch_parse(batches, paths)
        if (is.character(parse_fails)) {
            stop(parse_fails, call. = FALSE)
        }
    } else {
        batches <- df[[1]]
        paths <- df[[2]]
        if (any(!nzchar(batches) | !nzchar(paths))) {
            stop(
                "dCR path table cannot contain empty batches or paths",
                call. = FALSE
            )
        }
    }

    groups <- split(paths, batches)

    h <- utils::hashtab("identical", length(groups))
    set_hash(h, names(groups), groups)

    h
}

.valid_dcr_batch_parse <- function(batches, paths) {
    fail_idx <- which(batches == basename(paths) | !nzchar(batches))
    if (length(fail_idx) == 0) {
        return(TRUE)
    }

    fail_paths <- paths[fail_idx]
    base_msg <- paste0("could not parse batch ID from ", fail_paths[[1]])
    fail_counts <- length(fail_paths)
    if (fail_counts == 1) {
        return(base_msg)
    } else if (fail_counts == 2) {
        tail_msg <- "and 1 other"
    } else {
        tail_msg <- sprintf("and %d others", fail_counts - 1)
    }

    paste0(base_msg, " ", tail_msg)
}
