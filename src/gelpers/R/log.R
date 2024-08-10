#' Log to stderr
#'
#' @description
#' Write logging messages of different levels to stderr.
#'
#' @name logging
#'
#' @details
#' These functions enable logging messages of different levels (i.e. INFO,
#' WARN, and ERROR) to stderr without using signals. The different levels
#' refers to what tags the messages will have and not in the sense of
#' application verbosity or directing messages of different levels to different
#' places. Every message will be written to stderr with a string prefix
#' consisting of the current date, time and log level tag. In addition to the
#' different log levels, there are two different forms of log messages: single
#' line and block. A single line message will be written on a single line. A
#' block comprises a opening tag, the message to log, and a closing tag. The
#' opening tag is a single line consisting of the usual prefix followed by the
#' string "LOGBLOCK". The closing tag is a single line consisting of the string
#' "LOGBLOCK".
#'
#' @param msg A character vector of the message to log. \code{msg} can also be
#'   an R object for which \code{as.character()} is defined and works, in which
#'   case the result of \code{as.character(msg)} will be used as the message.
#'   If \code{msg} is a string without any newlines, the message will be
#'   written as a single line. If \code{msg} contains newlines or has a length
#'   greater than one, then the message will be written as a block. See
#'   details.
#'
#' @export
log_info <- function(msg) {
    .log(msg, "INFO")
}

#' @rdname logging
#' @export
log_warn <- function(msg) {
    .log(msg, "WARN")
}

#' @rdname logging
#' @export
log_error <- function(msg) {
    .log(msg, "ERROR")
}

.log <- function(msg, level) {
    if (!is.character(msg)) {
        msg <- tryCatch(
            suppressWarnings(as.character(msg)),
            error = function(cnd) NULL
        )

        if (!is.character(msg)) {
            stop("`as.character()` must work on `msg` in order to log it")
        }
    }

    msg <- msg[!is.na(msg) & nzchar(msg)]
    if (length(msg) == 0) {
        return()
    }

    is_block <- .has_newline(msg) || length(msg) > 1
    if (is_block) {
        .start_block_log(level)
        if (length(msg) == 1) {
            .cat2(msg)
        } else {
            cat(msg, file = stderr(), sep = "\n")
        }
        .end_block_log()
    } else {
        .cat2(.msg_prefix(level), msg)
    }
}

.start_block_log <- function(level) {
    .cat2(.msg_prefix(level), "LOGBLOCK")
}

.end_block_log <- function() {
    .cat2("LOGBLOCK")
}

.current_time <- function() {
    format(Sys.time(), format = "%Y-%m-%d %H:%M:%S %Z")
}

.has_newline <- function(x) {
    any(grepl("\n", x, fixed = TRUE))
}

.msg_prefix <- function(level) {
    sprintf("%s [%s]: ", .current_time(), level)
}

.cat2 <- function(...) {
    cat(..., "\n", file = stderr(), sep = "")
}
