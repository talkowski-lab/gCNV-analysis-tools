DEFAULT_OPTS <- list(cpus = 1L)

DEFAULT_ARGS <- list(DCRS = NULL, OUTPUT = NULL)

usage <- function() {
'Usage: Rscript check_sex.R [OPTIONS] DCRS OUTPUT

DCRS       A file with either paths, one per line, to the dCR matrices or two
           tab-separated coluns with the batch ID in the first column and path
           to the dCR matrix for that batch in the second column.
OUTPUT     Where to write the output.

Options:
  -h,--help           Print this message and exit.
  -p,--cpus[=]<cpus>  Number of processors to use. [1]
'
}

#' Write to stderr and exit
msg_and_exit <- function(msg, status = 1) {
    cat(msg, sep = "", file = stderr())
    quit(save = "no", status)
}

#' Parse a command-line option requiring an argument
#'
#' Three option types are supported:
#' * short option (e.g. -o <arg>)
#' * long option (e.g. --opt <arg>)
#' * long option with '=' (e.g. --opt=<arg>)
#'
#' @param x A character vector of all the command-line arguments i.e. as
#'   returned by `commandArgs(trailingOnly = TRUE)`.
#' @param m The match object returned by `regexec` from comparing the
#'   option against the option-specific regular expression.
#' @param i The index in `x` of the current option.
#' @returns A list of length two in which the first element is the argument of
#'   the parsed option and the second is the index of the next element of `x`
#'   parse. If an option does not provide an argument, the script will exit
#'   with an error message.
opt_arg <- function(x, m, i) {
    matches <- regmatches(x[[i]], m)[[1]]
    can_inc <- i + 1 <= length(x)
    opt <- ""
    if (length(matches) == 1) { # ^-a$
        opt <- matches[[1]]
    } else if (length(matches) == 3) { # ^(--long)(=(.*)?)?$
        if (!nzchar(matches[[3]])) {
            opt <- matches[[2]]
        } else {
            return(list(matches[[3]], i + 1))
        }
    } else { # ^((-a)|((--long)(=(.*)?)?))$
        if (!nzchar(matches[[7]])) {
            opt <- matches[[5]]
        } else {
            return(list(matches[[7]], i + 1))
        }
    }

    if (!can_inc) {
        msg_and_exit(sprintf("option %s requires an argument", opt))
    }

    list(x[[i + 1]], i + 2)
}

parse_args <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    if (length(argv) == 0) {
        msg_and_exit(usage(), 0)
    }

    opts <- DEFAULT_OPTS
    args <- DEFAULT_ARGS

    j <- 1
    i <- 1
    while (i <= length(argv)) {
        if (grepl("^((-h)|(--help))$", argv[[i]])) {
            msg_and_exit(usage(), 0)
        }

        if (all((m <- regexec("^((-p)|((--cpus)(=(.*)?)?))$", argv[[i]]))[[1]] != -1)) {
            x <- opt_arg(argv, m, i)
            opts$cpus <- tryCatch(
                as.integer(x[[1]]),
                warning = \(cnd) msg_and_exit(sprintf("could not convert '%s' to integer\n", x[[1]])),
                error = \(cnd) msg_and_exit(sprintf("could not convert '%s' to integer\n", x[[1]]))
            )
            if (opts$cpus < 1) {
                msg_and_exit("number of processors must be greater than 0\n")
            }
            i <- x[[2]]
            next
        }

        if (grepl("^-", argv[[i]])) {
            msg_and_exit(sprintf("unknown option %s\n", argv[[i]]))
        } else {
            if (j <= length(args)) {
                args[[j]] <- argv[[i]]
                j <- j + 1
                i <- i + 1
            } else {
                msg_and_exit("too many arguments\n")
            }
        }
    }

    missing_args <- Filter(is.null, args)
    if (length(missing_args) > 0) {
        msg_and_exit(sprintf("required argument(s) %s missing\n",
                             paste0(names(missing_args), collapse = ", ")))
    }

    append(args, opts)
}

args <- parse_args()

suppressPackageStartupMessages(library(gelpers))
suppressPackageStartupMessages(library(parallel))

# Min and max expected average chrX dCR values for XX and XY samples
XX_MIN <- 1.9
XX_MAX <- 2.1
XY_MIN <- 0.9
XY_MAX <- 1.1

# Maximum sequence length supported by tabix. A hack to get an entire
# chromosome from a dCR matrix through Rsamtools.
TABIX_MAX_SEQLEN <- 536870912L

#' Estimate chrX and chrY ploidy for a single batch
ploidy_from_batch <- function(batch, dcr_map) {
    hg38_x <- gregion("chrX", 1, TABIX_MAX_SEQLEN)
    hg38_y <- gregion("chrY", 1, TABIX_MAX_SEQLEN)
    hg19_x <- gregion("X", 1, TABIX_MAX_SEQLEN)
    hg19_y <- gregion("Y", 1, TABIX_MAX_SEQLEN)

    dcr_paths <- gethash(dcr_map, batch)
    dcr_x <- tryCatch(
        get_batch_dcr(hg38_x, dcr_paths),
        dcr_parse_error = function(cnd) cnd,
        error = function(cnd) cnd
    )

    if (inherits(dcr_x, "dcr_parse_error")) {
        dcr_x <- get_batch_dcr(hg19_x, dcr_paths)
        dcr_y <- get_batch_dcr(hg19_y, dcr_paths)
    }

    if (!inherits(dcr_x, "error")) {
        dcr_y <- get_batch_dcr(hg38_y, dcr_paths)
    } else {
        stop(dcr_x)
    }

    if (nrow(dcr_x) == 0) {
        stop(sprintf("batch %s missing chrX dCR regions", batch))
    }

    if (nrow(dcr_y) == 0) {
        stop(sprintf("batch %s missing chrY dCR regions", batch))
    }

    setDT(dcr_x)
    setDT(dcr_y)
    mean_x <- dcr_x[, lapply(.SD, mean), .SDcols = !c("chr", "start", "end")]
    ploidy_x <- melt(mean_x,
                     measure.vars = colnames(mean_x),
                     variable.name = "sample",
                     value.name = "chrX_ploidy",
                     variable.factor = FALSE)
    mean_y <- dcr_y[, lapply(.SD, mean), .SDcols = !c("chr", "start", "end")]
    ploidy_y <- melt(mean_y,
                     measure.vars = colnames(mean_y),
                     variable.name = "sample",
                     value.name = "chrY_ploidy",
                     variable.factor = FALSE)

    ploidy_x[ploidy_y, on = "sample"]
}

#' Estimate chrX and chrY ploidy across cohort
#'
#' @param x `character` of batch IDs
#' @param dcr_map A hash map of batch ID to dCR matrix paths.
#' @param nproc Number of cpus to use.
#' @returns A `data.table` of all the chrX and chrY ploidy estimates.
estimate_ploidy <- function(x, dcr_map, nproc = 1L) {
    x <- unique(x)
    old_dtthreads <- getDTthreads()
    setDTthreads(1)
    ploidy <- mclapply(x, ploidy_from_batch, dcr_map = dcr_map, mc.cores = nproc)
    setDTthreads(old_dtthreads)

    rbindlist(ploidy)
}

# Read inputs -----------------------------------------------------------------
setDTthreads(args$cpus)

log_info("reading dCR paths")
dcrs <- read_dcr_list(args$DCRS)
i <- 0
batch_ids <- character(numhash(dcrs))
maphash(dcrs, \(k, v) { i <<- i + 1; batch_ids[i] <<- k })

log_info("estimating ploidy from dCR")
ploidy <- estimate_ploidy(batch_ids, dcrs, args$cpus)
ploidy[, c("chrX_ploidy", "chrY_ploidy") := list(sprintf("%0.2f", chrX_ploidy), sprintf("%0.2f", chrY_ploidy))]

log_info("writing output")
fwrite(ploidy, file = args$OUTPUT, sep = "\t", quote = FALSE)
