#' Read a pedigree file
#'
#' @description
#' A read a pedigree file and peform minor cleaning.
#'
#' @section File format:
#' The file must be tab-delimited with six columns in the following order:
#' \enumerate{
#'   \item family_id
#'   \item sample_id
#'   \item paternal_id
#'   \item maternal_id
#'   \item sex
#'   \item phenotype
#' }
#' The file must have a header, but it will be ignored and columns will be named
#' according to the above. Only the phenotypes -9, 0, 1, and 2 are accepted with
#' -9 and 0 meaning unknown, 1 meaning unaffected and 2 meaning affected. All
#' other phenotypes will be silently converted to 0 and then -9 and 0 will be
#' converted to \code{NA}. Sex can be 1 for male and 2 for female. All other sex
#' codes will be converted to \code{NA}.
#'
#' @section Sample IDs:
#' Any empty or duplicate sample IDs will be removed with a warning. Any empty
#' family, paternal, or maternal IDs will be converted to \code{NA}. It is
#' better for the user to properly sanitize the file than to rely on this
#' function to perform the cleaning.
#'
#' @export
#' @param path \code{character(1)} Path to the pedigree file.
#' @returns \code{data.frame}.
read_pedigree <- function(path) {
    assert(is_string(path))

    ped <- data.table::fread(file = path,
                             header = TRUE,
                             sep = "\t",
                             col.names = c("family_id", "sample_id",
                                           "paternal_id", "maternal_id",
                                           "sex", "phenotype"),
                             colClasses = rep("character", 6),
                             na.strings = "NA")

    missing_ids <- is.na(ped$sample_id) | !nzchar(ped$sample_id)
    if (any(missing_ids)) {
        warning(
            sprintf("removing %d missing sample(s) from pedigree", sum(missing_ids)),
            call. = FALSE
        )
    }
    ped <- ped[!missing_ids, ]

    family_id <- NULL
    sample_id <- NULL
    paternal_id <- NULL
    maternal_id <- NULL
    sex <- NULL
    phenotype <- NULL

    valid_phen <- c("-9", "0", "1", "2")
    ped[!phenotype %in% valid_phen, phenotype := "0"]
    ped[, phenotype := as.integer(phenotype)]
    ped[phenotype == 0L | phenotype == -9L, phenotype := NA_integer_]

    ped[!sex %in% c("1", "2"), sex := NA_character_]
    ped[, sex := as.integer(sex)]

    ped[!nzchar(family_id), family_id := NA_character_]
    ped[!nzchar(paternal_id), paternal_id := NA_character_]
    ped[!nzchar(maternal_id), maternal_id := NA_character_]

    # Some pedigrees have duplicate sample IDs with different attributes. We
    # want to remove duplicates, but keep families together so we sort on both
    # family ID and sample ID before taking a unique.
    data.table::setorder(ped, family_id, sample_id, na.last = TRUE)
    dups <- duplicated(ped, by = "sample_id")
    if (any(dups)) {
        warning(
            sprintf("removing %d duplicate sample(s) from pedigree", sum(dups)),
            call. = FALSE
        )
        ped <- ped[!dups, ]
    }

    ped
}
