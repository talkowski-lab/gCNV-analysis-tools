#' Read a gCNV callset
#'
#' @description
#' Read a callset produced by the gCNV pipeline and check for required columns.
#'
#' @details
#' There are a few columns that every gCNV callset should have in order to do
#' any useful analysis on it. As the gCNV pipeline changes and postprocessing
#' is done on callsets, the columns can change as well. To ensure that a
#' callset has these required columns, we check that they exist after reading
#' the file (which is a tab-delimited file). Currently the required columns are
#' "chr", "start", "end", "sample", "batch", "svtype", "variant_name", and
#' "CN". If any of these columns don't exist, an error will be signaled.
#' The file can be compressed with one of the compression methods supported by
#' base R.
#'
#' @export
#' @param path \code{character(1)}. Path to the callset.
#' @returns \code{data.table} keyed on the "sample" column.
read_callset <- function(path) {
    assert(is_string(path))
    callset <- data.table::fread(file = path,
                                 sep = "\t",
                                 header = TRUE)

    req_cols <- c("chr", "start", "end", "sample", "batch", "svtype",
                  "variant_name", "CN")
    assert(df_has_columns(callset, .cols = req_cols))

    data.table::setkey(callset, sample)

    callset
}
