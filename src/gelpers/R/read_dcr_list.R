#' Read list of dCR matrix paths
#'
#' @description
#' Read a list of dCR matrix paths and create a hash table from gCNV batch to
#' dCR matrix path.
#'
#' @details
#' The paths must be of the form
#' "path/to/\{BATCH\}_((CASE)|(COHORT)).dcr.bed.gz" which is the default for
#' gCNV. The CASE and COHORT files will be combined for a batch such that
#' querying the hash table for a batch will return a character vector with the
#' paths to both the CASE and COHORT matrices, if they are available.
#'
#' @export
#' @param path Path to the file of dCR matrix paths. There should be one dCR
#'   matrix path per line.
#' @returns \code{hashtab} in which the keys are batch IDs and the values are
#'   paths to the dCR matrices for the batch.
read_dcr_list <- function(path) {
    assert(is_string(path))
    paths <- unique(readLines(path))

    if (length(paths) == 0) {
        return(utils::hashtab("identical"))
    }
    batches <- sub("_((CASE)|(COHORT)).dcr.bed.gz", "", basename(paths))
    groups <- split(paths, batches)

    h <- utils::hashtab("identical", length(groups))
    set_hash(h, names(groups), groups)

    h
}
