# de novo CNV plotting pipeline

CHILD_COLOR <- "#84A955"
FATHER_COLOR <-"#965DA7"
MOTHER_COLOR <- "#BC5D41"

usage <- function() {
'Usage: Rscript plot_cnv_evidence.R DENOVO BINS PED DCRS OUTDIR

DENOVO   de novo callset
BINS     Genomic intervals file used by the gCNV pipeline
PED      Pedigree file
DCRS     A file with either the paths, one per line, to the dCR matrices or
         two tab-separated columns with the batch ID in the first column and
         path to the dCR matrix for that batch in the second column
OUTDIR   Path to the output directory. The directory must not exist
'
}

parse_args <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    if (length(argv) != 5) {
        cat(usage(), sep = "", file = stderr())
        quit(save = "no", 2)
    }

    args <- as.list(argv)
    names(args) <- c("DENOVO", "BINS", "PED", "DCRS", "OUTDIR")

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

# Plot the case where the child has only one dCR interval overlapping its SV
plot_single_interval <- function(x, dcr, main) {
    dcr_point <- dcr[!is.na(child) & !in_flank, ]
    boxplot(as.double(dcr_point[, .SD, .SDcols = patterns("^bg_")]),
            horizontal = TRUE,
            outline = FALSE,
            xlab = "Denoised Coverage",
            ylab = "",
            ylim = c(0, 5),
            main = main)

    points(dcr_point$child,
           1 + runif(1, -0.07, 0.07),
           col = CHILD_COLOR,
           pch = if (x$phenotype == 2) 15 else 19,
           cex = 1.5)

    if (!is.na(dcr_point$father)) {
        points(dcr_point$father,
               1 + runif(1, -0.07, 0.07),
               col = FATHER_COLOR,
               pch = if (x$paternal_phenotype == 2) 15 else 19,
               cex = 1.5)
    }

    if (!is.na(dcr_point$mother)) {
        points(dcr_point$mother,
               1 + runif(1, -0.07, 0.07),
               col = MOTHER_COLOR,
               pch = if (x$maternal_phenotype == 2) 15 else 19,
               cex = 1.5)
    }

    legend("topleft",
           legend = c("offspring", "father", "mother"),
           col = c(CHILD_COLOR, FATHER_COLOR, MOTHER_COLOR),
           pch = 19,
           bty = "n")
    legend("topright",
           legend = c("Affected", "Unaffected"),
           pch = c(15, 19),
           bty = "n")
}

plot_sample_dcr_line <- function(x, color, phenotype) {
    lty = if (phenotype == 2) "solid" else "dashed"
    pch = if (phenotype == 2) 15 else 19
    lines(seq_along(x), x, lwd = 3, lty = lty, col = color)

    idx <- which(!is.na(x))
    na_after <- idx[-1] - idx[-length(idx)] > 1
    na_before <- which(na_after) + 1
    singleton_idx <- sort(intersect(idx[na_after], idx[na_before]))
    points(singleton_idx, x[singleton_idx], cex = 1, pch = pch, col = color)

    # Handle intervals at the start and end of vector
    # Only guarantee is that there are at least two intervals because this
    # function is only called if the child has at least two non-NA, non-flanking
    # intervals
    if (!is.na(x[[1]]) && is.na(x[[2]])) {
        points(1, x[[1]], cex = 1, pch = pch, col = color)
    }
    if (!is.na(x[[length(x)]]) && is.na(x[[length(x) - 1]])) {
        points(length(x), x[[length(x)]], cex = 1, pch = pch, col = color)
    }
}

plot_mult_interval <- function(x, dcr, main) {
    blocks <- split_dcr_blocks(dcr$in_flank)
    # If there are at least two consecutive intervals with values, then matplot
    # will be able to plot at least one line and because all the background
    # samples come from the same batch as the child, if the child has two
    # consecutive intervals, so will the background samples.
    child_has_lines <- has_consecutive_values(dcr$child)
    matplot(data.matrix(dcr[, .SD, .SDcols = patterns("^bg_")]),
            type = if (child_has_lines) "l" else "p",
            lty = 1,
            pch = 19,
            ylim = c(0, 5),
            col = "#77777777",
            xlab = "Interval",
            ylab = "Denoised Coverage",
            main = main)

    plot_sample_dcr_line(dcr$child, CHILD_COLOR, x$phenotype)
    plot_sample_dcr_line(dcr$father, FATHER_COLOR, x$paternal_phenotype)
    plot_sample_dcr_line(dcr$mother, MOTHER_COLOR, x$maternal_phenotype)

    if (length(blocks$left_flank) > 0) {
        rect(1, 0, blocks$left_flank[[2]] + 0.5, 5, col = "#33333311", border = NA)
        text(1, 0.3,
            labels = paste0(blocks$left_flank[[2]] - blocks$left_flank[[1]] + 1,
                            " intervals: ",
                            range2kb(dcr$coords[blocks$left_flank[[1]], "start"],
                                     dcr$coords[blocks$left_flank[[2]], "end"]),
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
                           range2kb(dcr$coords[blocks$right_flank[[1]], "start"],
                                    dcr$coords[blocks$right_flank[[2]], "end"]),
                           "kb"),
           pos = 2,
           offset = 0.3,
           cex = 0.9)
    }

    legend("topleft",
           legend = c("Offspring", "Father", "Mother"),
           col = c(CHILD_COLOR, FATHER_COLOR, MOTHER_COLOR),
           pch = 19,
           bty = "n")
    legend("topright",
           legend = c("Affected", "Unaffected"),
           lty = c("solid", "dashed"),
           pch = c(15, 19),
           bty = "n")
}

plot_dcr <- function(x, dcr, outfile) {
    png(filename = paste0(outfile, ".png"), width = 647, height = 400)
    main <- paste0(x$chr, ":", x$start, "-", x$end, "\n",
                   x$sample, ";", x$variant_name, ";", x$svtype)
    old_par <- par(font.lab = 2, mar = c(4, 4, 4, 2) + 0.1)
    if (nrow(dcr[!is.na(child) & !in_flank, ]) == 1) {
        plot_single_interval(x, dcr, main)
    } else {
        plot_mult_interval(x, dcr, main)
    }
    par(old_par)
    dev.off()
}

get_trio_dcr <- function(trio,
                         dcr_paths,
                         bins,
                         bg = 50,
                         pad = 0.2) {
    sv_region <- gregion(trio$chr, trio$start, trio$end)
    expanded_region <- expand_region_by_bins(sv_region, bins, pad)

    trio_ids <- c(trio$sample, trio$paternal_id, trio$maternal_id)
    batches <- c(trio$batch, trio$paternal_batch, trio$maternal_batch)
    coord_cols <- c("chr", "start", "end")

    batch_groups <- split(trio_ids, batches)
    trio_dcr <- vector(mode = "list", length = 3)
    names(trio_dcr) <- c("child", "father", "mother")
    bg_dcr <- NULL
    for (i in seq_along(batch_groups)) {
        batch <- names(batch_groups)[[i]]
        samples <- batch_groups[[i]]
        dcr <- get_samples_dcr(expanded_region,
                               gethash(dcr_paths, batch),
                               batch_groups[[i]],
                               include_bg = trio$sample %in% samples,
                               squeeze = TRUE,
                               reduce = TRUE)

        if (trio$sample %in% samples) {
            trio_dcr$child <- as.data.table(dcr[c(coord_cols, trio$sample)])
            setkey(trio_dcr$child, chr, start, end)
            setnames(trio_dcr$child, trio$sample, "child")
            bg_samples <- colnames(dcr)[!colnames(dcr) %in% c(coord_cols, samples)]
            if (length(bg_samples) == 0) {
                stop("child batch has no background samples", call. = FALSE)
            }
            bg_samples <- sample(bg_samples, min(length(bg_samples), bg))
            bg_dcr <- as.data.table(dcr[c(coord_cols, bg_samples)])
            setnames(bg_dcr, bg_samples, paste0("bg_", seq_along(bg_samples)))
            setkey(bg_dcr, chr, start, end)
        }

        if (trio$paternal_id %in% samples) {
            trio_dcr$father <- as.data.table(dcr[c(coord_cols, trio$paternal_id)])
            setkey(trio_dcr$father, chr, start, end)
            setnames(trio_dcr$father, trio$paternal_id, "father")
        }

        if (trio$maternal_id %in% samples) {
            trio_dcr$mother <- as.data.table(dcr[c(coord_cols, trio$maternal_id)])
            setkey(trio_dcr$mother, chr, start, end)
            setnames(trio_dcr$mother, trio$maternal_id, "mother")
        }
    }

    # Assume ranges are disjoint
    merged_dcr <- Reduce(\(x, y) merge(x, y, all = TRUE), trio_dcr)
    merged_dcr <- merge(merged_dcr, bg_dcr, all.x = TRUE)
    merged_dcr[, in_flank := end < sv_region$start | start > sv_region$end]

    merged_dcr
}

plot_denovo <- function(x, dcr_paths, bins, outdir) {
    log_info(sprintf("plotting sample '%s' at variant '%s'", x$sample, x$variant))
    dcr <- tryCatch(
        get_trio_dcr(x, dcr_paths, bins),
        error = function(cnd) {
            log_error(conditionMessage(cnd))
            NULL
        })

    if (is.null(dcr)) {
        log_warn("could not make trio dCR matrix")
        return()
    }

    if (nrow(dcr[!in_flank & !is.na(child), ]) == 0) {
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
log_info("reading de novo callset")
denovo <- read_callset(args$DENOVO)
denovo[,  `:=`(chr = as.character(chr),
                start = as.integer(start),
                end = as.integer(end))]
is_hg19 <- any(c(as.character(1:22), "X", "Y") %in% denovo$chr)
log_info("reading bins")
bins <- read_gcnv_bins(args$BINS, reduce = is_hg19)
log_info("reading pedigree")
pedigree <- read_pedigree(args$PED)
log_info("reading dCR paths")
dcr_paths <- read_dcr_list(args$DCRS)
setup_outdir(args$OUTDIR)

# Merge in phenotype ----------------------------------------------------------
phen <- pedigree[, list(sample = sample_id, phenotype = phenotype)]
# Just set all missing phenotypes to unaffected
phen[is.na(phenotype), phenotype := 1]
setkey(phen, sample)

denovo <- phen[denovo, on = "sample"]
denovo <- phen[, list(paternal_id = sample, paternal_phenotype = phenotype)][denovo, on = "paternal_id"]
denovo <- phen[, list(maternal_id = sample, maternal_phenotype = phenotype)][denovo, on = "maternal_id"]

# Plot de novo calls ----------------------------------------------------------
log_info("Plotting de novo calls")
for (i in seq_len(nrow(denovo))) {
    plot_denovo(as.list(denovo[i, ]), dcr_paths, bins, args$OUTDIR)
}
log_info("Completed de novo CNV plotting")
