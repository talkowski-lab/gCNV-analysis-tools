#' Read a pedigree file
#'
#' @description
#' A read a pedigree file and peform minor cleaning.
#'
#' @details
#' The file must be tab-delimited with six columns in the following order:
#' \enumerate{
#'   \item family_id
#'   \item sample_id
#'   \item paternal_id
#'   \item maternal_id
#'   \item sex
#'   \item phenotype
#' }
#' The file must have a header, but it will be ignored and columns will be
#' named according to the above. Only the phenotypes -9, 0, 1, and 2 are
#' accepted with -9 and 0 meaning unknown, 1 meaning unaffected and 2 meaning
#' affected. All other phenotypes will be silently converted to 0 and then -9
#' and 0 will be converted to \code{NA}. Sex can be 1 for male and 2 for female.
#' All other sex codes will be converted to \code{NA}. Any empty or duplicate
#' sample IDs will be silently removed. Any empty family, paternal, or maternal
#' IDs will be converted to \code{NA}. It is better for the user to properly
#' sanitize the file than to rely on this function to perform the cleaning.
#'
#' @export
#' @param path \code{character(1)} Path to the pedigree file.
#' @returns \code{data.frame}.
read_pedigree <- function(path) {
    assert(is_string(path))
    x <- utils::read.table(
        file = path,
        header = TRUE,
        sep = "\t",
        col.names = c(
            "family_id", "sample_id", "paternal_id", "maternal_id", "sex",
            "phenotype"
        ),
        colClasses = rep("character", 6)
    )
    x <- x[nzchar(x$sample_id), ]

    valid_phen <- c("-9", "0", "1", "2")
    x[!x$phenotype %in% valid_phen, "phenotype"] <- "0"
    x$phenotype <- as.integer(x$phenotype)
    x[x$phenotype == 0L | x$phenotype == -9L, "phenotype"] <- NA

    x[!x$sex %in% c("1", "2"), "sex"] <- NA
    x$sex <- as.integer(x$sex)

    # Some pedigrees have duplicate sample IDs with different attributes. We
    # want to remove duplicates, but keep families together so we sort on both
    # family ID and sample ID before taking a unique.
    x <- x[order(x$family_id, x$sample_id), ]
    x <- x[!duplicated(x["sample_id"]), ]

    x
}
