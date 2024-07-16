get_hash <- function(h, keys, nomatch = NULL) {
    v <- lapply(keys, \(k) utils::gethash(h, k, nomatch))

    v
}

set_hash <- function(h, keys, values) {
    mapply(
        \(k, v) utils::sethash(h, k, v),
        keys, values,
        SIMPLIFY = FALSE, USE.NAMES = FALSE
    )

    invisible(values)
}
