cleanup_file <- function(x) {
    out <- tryCatch(
        {
            file.remove(x)
        },
        error = function(cnd) TRUE,
        warning = function(cnd) TRUE
    )

    invisible(out)
}

overlaps <- function(x1, y1, x2, y2) {
    x1 <= y2 & x2 <= y1
}
