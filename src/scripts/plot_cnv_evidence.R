# CNV plotting pipeline

SV_COLORS <- c(DEL = "#D43925", DUP = "#2376B2")

usage <- function() {
'Usage: Rscript plot_cnv_evidence.R CNVS BINS DCRS OUTDIR

CNVS    CNVs to plot
BINS    Genomic intervals file used by the gCNV pipeline
DCRS    A file with either the paths, one per line, to the dCR matrices or
        two tab-separated columns with the batch ID in the first column and
        path to the dCR matrix for that batch in the second column
OUTDIR  Path to the output directory. The directory must not exist
'
}

parse_args <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    if (length(argv) != 4) {
        cat(usage(), sep = "", file = stderr())
        quit(save = "no", 2)
    }

    args <- as.list(argv)
    names(args) <- c("CNVS", "BINS", "DCRS", "OUTDIR")

    args
}

# Create the output directory for the plots. Signal error if it already exists
setup_outdir <- function(path) {
    if (dir.exists(path)) {
        stop("Output directory already exists", call. = FALSE)
    }

    success <- dir.create(path)
    if (!success) {
        stop("Failed to create output directory", call. = FALSE)
    }

    invisible(path)
}

# Given a logical vector parallel to the rows of a dCR matrix, indicating
# whether or not each row is part of the flanking region, return a list of
# indicies marking each row as part of the left flank, the SV or the right
# flank. It is assumed that the vector has contiguous blocks of logicals with
# at most two blocks of TRUE's and blocks of TRUE's only at the head and tail
# of the vector.
split_dcr_blocks <- function(x) {
    bls <- rle(x)
    if (length(bls$lengths) == 3) {
        left_flank <- c(1L, bls$lengths[[1]])
        sv <- c(left_flank[[2]] + 1L, left_flank[[2]] + bls$lengths[[2]])
        right_flank <- c(sv[[2]] + 1L, sv[[2]] + bls$lengths[[3]])
    } else if (length(bls$lengths) == 2) {
        if (bls$values[[1]]) {
            left_flank <- c(1L, bls$lengths[[1]])
            sv <- c(left_flank[[2]] + 1L, left_flank[[2]] + bls$lengths[[2]])
            right_flank <- integer(0)
        } else {
            left_flank <- integer(0)
            sv <- c(1L, bls$lengths[[1]])
            right_flank <- c(sv[[2]] + 1L, sv[[2]] + bls$lengths[[2]])
        }
    } else {
        left_flank <- integer(0)
        sv <- c(1L, length(x))
        right_flank <- integer(0)
    }

    list(left_flank = left_flank, sv = sv, right_flank = right_flank)
}

# Convert a genomic range to kilobases
range2kb <- function(start, end) {
    round((end - start + 1) / 1000, 1)
}

# Test if a vector has at least two consecutive non-missing values
has_consecutive_values <- function(x) {
    if (length(x) < 2) {
        return(FALSE)
    }

    idx <- which(!is.na(x))
    if (length(idx) == 0) {
        return(FALSE)
    }

    any(idx[-1] - idx[-length(idx)] == 1)
}

# Plot the case where the carrier has only one dCR interval overlapping its SV
plot_single_interval <- function(x, dcr, main) {
    dcr_point <- dcr[!is.na(carrier) & !in_flank, ]
    boxplot(as.double(dcr_point[, .SD, .SDcols = patterns("^bg_")]),
            horizontal = TRUE,
            outline = FALSE,
            xlab = "Denoised Coverage",
            ylab = "",
            ylim = c(0, 5),
            main = main)

    points(dcr_point$carrier,
           1 + runif(1, -0.1, 0.1),
           col = SV_COLORS[[x$svtype]],
           pch = 19,
           cex = 1.5)
}

plot_sample_dcr_line <- function(x, svtype) {
    color = SV_COLORS[[svtype]]
    lines(seq_along(x), x, lwd = 3, lty = "solid", col = color)

    idx <- which(!is.na(x))
    na_idx <- which(is.na(x))
    na_after <- idx[(idx + 1) %in% na_idx]
    na_before <- idx[(idx - 1) %in% na_idx]
    singleton_idx <- sort(intersect(na_after, na_before))
    points(singleton_idx, x[singleton_idx], cex = 1, pch = 19, col = color)

    # Handle intervals at the start and end of vector
    # Only guarantee is that there are at least two intervals because this
    # function is only called if the sample has at least two non-NA, non-flanking
    # intervals
    if (!is.na(x[[1]]) && is.na(x[[2]])) {
        points(1, x[[1]], cex = 1, pch = 19, col = color)
    }
    if (!is.na(x[[length(x)]]) && is.na(x[[length(x) - 1]])) {
        points(length(x), x[[length(x)]], cex = 1, pch = 19, col = color)
    }
}

plot_mult_interval <- function(x, dcr, main) {
    blocks <- split_dcr_blocks(dcr$in_flank)
    # If there are at least two consecutive intervals with values, then matplot
    # will be able to plot at least one line and because all the background
    # samples come from the same batch as the carrier, if the carrier has two
    # consecutive intervals, so will the background samples.
    sample_has_lines <- has_consecutive_values(dcr$carrier)
    matplot(data.matrix(dcr[, .SD, .SDcols = patterns("^bg_")]),
            type = if (sample_has_lines) "l" else "p",
            lty = 1,
            pch = 19,
            ylim = c(0, 5),
            col = "#77777777",
            xlab = "Interval",
            ylab = "Denoised Coverage",
            xaxt = "n",
            main = main)

    plot_sample_dcr_line(dcr$carrier, x$svtype)
    if (length(blocks$left_flank) > 0) {
        rect(1, 0, blocks$left_flank[[2]] + 0.5, 5, col = "#33333311", border = NA)
        text(1, 0.3,
            labels = paste0(blocks$left_flank[[2]] - blocks$left_flank[[1]] + 1,
                            " intervals: ",
                            range2kb(dcr[blocks$left_flank[[1]], start],
                                     dcr[blocks$left_flank[[2]], end]),
                            "kb"),
            pos = 4,
            offset = 0.3,
            cex = 0.9)
    }
    if (length(blocks$right_flank) > 0) {
      rect(blocks$right_flank[[1]] - 0.5, 0,
           blocks$right_flank[[2]], 5,
           col = "#33333311",
           border = NA)
      text(blocks$right_flank[[2]], 0.3,
           labels = paste0(blocks$right_flank[[2]] - blocks$right_flank[[1]] + 1,
                           " intervals: ",
                           range2kb(dcr[blocks$right_flank[[1]], start],
                                    dcr[blocks$right_flank[[2]], end]),
                           "kb"),
           pos = 2,
           offset = 0.3,
           cex = 0.9)
    }

    xticks <- axTicks(1)
    xticks <- xticks[trunc(xticks) == xticks]
    axis(1, at = xticks, labels = sprintf("%d", xticks))
}

plot_dcr <- function(x, dcr, outfile) {
    png(filename = paste0(outfile, ".png"), width = 647, height = 400)
    main <- paste0(x$chr, ":", x$start, "-", x$end, "\n",
                   x$sample, ";", x$variant_name, ";", x$svtype)
    old_par <- par(font.lab = 2, mar = c(4, 4, 4, 2) + 0.1)
    if (nrow(dcr[!is.na(carrier) & !in_flank, ]) == 1) {
        plot_single_interval(x, dcr, main)
    } else {
        plot_mult_interval(x, dcr, main)
    }
    par(old_par)
    dev.off()
}

get_dcr <- function(x,
                    dcr_paths,
                    bins,
                    bg = 50,
                    pad = 0.2) {
    sv_region <- gregion(x$chr, x$start, x$end)
    expanded_region <- expand_region_by_bins(sv_region, bins, pad)
    coord_cols <- c("chr", "start", "end")

    dcr <- get_samples_dcr(expanded_region,
                           gethash(dcr_paths, x$batch),
                           x$sample,
                           include_bg = TRUE) |>
        as.data.table()
    bg_samples <- colnames(dcr)[!colnames(dcr) %in% c(coord_cols, x$sample)]
    if (length(bg_samples) == 0) {
        stop("batch has no background samples", call. = FALSE)
    }
    bg_samples <- sample(bg_samples, min(length(bg_samples), bg))
    keep_cols <- c(coord_cols, x$sample, bg_samples)
    dcr <- dcr[, ..keep_cols]
    setnames(dcr, x$sample, "carrier")
    setnames(dcr, bg_samples, paste0("bg_", seq_along(bg_samples)))
    setkey(dcr, chr, start, end)

    dcr[, in_flank := end < sv_region$start | start > sv_region$end]

    dcr
}

plot_cnv <- function(x, dcr_paths, bins, outdir) {
    log_info(sprintf("plotting sample '%s' at variant '%s'", x$sample, x$variant))
    dcr <- tryCatch(
        get_dcr(x, dcr_paths, bins),
        error = function(cnd) {
            log_error(conditionMessage(cnd))
            NULL
        })

    if (is.null(dcr)) {
        log_warn("could not get sample dCR matrix")
        return()
    }

    if (nrow(dcr[!in_flank & !is.na(carrier), ]) == 0) {
        log_warn("dCR matrix does not have any non-missing intervals")
        return()
    }

    outfile <- file.path(outdir,
                         paste(x$sample, x$variant_name, x$svtype, sep = "__"))
    plot_dcr(x, dcr, outfile)
}

# Parse command-line arguments ------------------------------------------------
args <- parse_args()

suppressPackageStartupMessages(library(gelpers))
suppressPackageStartupMessages(library(GenomicRanges))

set.seed(42)

# Read inputs -----------------------------------------------------------------
log_info("reading CNV callset")
cnvs <- read_callset(args$CNVS)
cnvs[, `:=`(chr = as.character(chr),
            start = as.integer(start),
            end = as.integer(end))]

is_hg19 <- any(c(as.character(1:22), "X", "Y") %in% cnvs$chr)
log_info("reading bins")
bins <- read_gcnv_bins(args$BINS, reduce = is_hg19)
log_info("reading dCR paths")
dcr_paths <- read_dcr_list(args$DCRS)
setup_outdir(args$OUTDIR)

# Plot CNV calls --------------------------------------------------------------
log_info("plotting CNV evidence")
for (i in seq_len(nrow(cnvs))) {
    plot_cnv(as.list(cnvs[i, ]), dcr_paths, bins, args$OUTDIR)
}
log_info("completed CNV plotting")
