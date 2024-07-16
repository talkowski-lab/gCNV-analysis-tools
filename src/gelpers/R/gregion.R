#' Create a \code{gregion} object
#'
#' @description
#' Create a \code{gregion} object.
#'
#' @details
#' A \code{gregion} object represents a single genomic region. It is useful for
#' working with genomic coordinates as a single object because it is named list
#' with \code{chr}, \code{start}, and \code{end} components.
#'
#' @export
#' @param chr Name of the sequence.
#' @param start Start of the region. Should be 1-based.
#' @param end End of the region. Should be closed.
#' @returns A new \code{gregion} object.
gregion <- function(chr, start, end) {
    validate_gregion(new_gregion(chr, start, end))
}

#' Create a string representation of an object.
#'
#' @description
#' S3 generic to create strings from objects.
#'
#' @details
#' Methods should accept a single object and return a single string (character
#' vector of length 1).
#'
#' @export
#' @param x An object.
#' @returns A string representing the object.
to_string <- function(x) {
    UseMethod("to_string")
}

#' @export
to_string.gregion <- function(x) {
    sprintf("%s:%d-%d", x$chr, x$start, x$end)
}

#' @export
print.gregion <- function(x, ...) {
    cat(to_string(x), "\n")
}

new_gregion <- function(chr, start, end) {
    assert(is_string(chr))
    assert(is_integerish(start))
    assert(is_integerish(end))

    structure(
        list(chr = chr, start = as.integer(start), end = as.integer(end)),
        class = "gregion"
    )
}

validate_gregion <- function(x) {
    if (length(x$chr) != 1) {
        stop("`x$chr` must be a string", call. = FALSE)
    }

    if (length(x$start) != 1) {
        stop("`x$start` must be a single integer", call. = FALSE)
    }

    if (length(x$end) != 1) {
        stop("`x$end` must be a single integer", call. = FALSE)
    }

    if (x$start < 1) {
        stop("`x$start` must be >= 1", call. = FALSE)
    }

    if (x$end < 1) {
        stop("`x$end` must be >= 1", call. = FALSE)
    }

    if (x$end < x$start) {
        stop("`x$start` must be >= `x$end`", call. = FALSE)
    }

    x
}
