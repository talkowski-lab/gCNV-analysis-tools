#' Read a gCNV callset
#'
#' @description
#' Read a callset produced by the gCNV pipeline and ensure that it has the
#' required columns for workflows in this package.
#'
#' @details
#' Many of the functions in this package that operate on a gCNV callset
#' \code{data.table} assume that there are certain columns available because
#' the columns are required for the function to work. The most common ones are
#' "chr", "start", "end", "sample", "batch", "svtype", and "variant_name". gCNV
#' should include all of these columns as part of its output, but sometimes
#' column names change or callset processing drops some columns. Therefore,
#' \code{read_callset()} checks for these columns and throws an error if any
#' are missing.
#'
#' The callset should be a tab-delimited file with a header line. The file can
#' be compressed, but then the \pkg{R.utils} package is required to read it.
#'
#' @export
#' @param path \code{character(1)}. Path to the callset.
#' @returns \code{data.table} or signals an error if required columns are
#'   missing.
read_callset <- function(path) {
    assert(is_string(path))
    x <- data.table::fread(
        file = path,
        sep = "\t",
        header = TRUE,
        showProgress = FALSE
    )

    req_cols <- c(
        "chr", "start", "end", "sample", "batch", "svtype", "variant_name", "CN"
    )
    assert(df_has_columns(x, .cols = req_cols))
    data.table::setkey(x, sample)

    x
}
