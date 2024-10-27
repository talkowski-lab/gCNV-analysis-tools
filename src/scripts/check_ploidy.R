# Ploidy checker

TABIX_MAX_SEQLEN <- 536870912L

usage <- function() {
'Usage: Rscript make_ploidy_evidence.R SAMPLES CONTIGS DCRS OUTDIR

SAMPLES  A two-column, tab-delimited table with all the samples to process.
         The columns must be batch ID and sample ID.
CONTIGS  A list of contigs to check.
DCRS     A file with either the paths, one per line, to the dCR matrices
         or two tab-separated columns with the batch ID in the first
         column and path to the dCR matrix for that batch in the second
         column.
OUTDIR   Output directory. The directory must not exist.
'
}

parse_args <- function() {
    argv <- as.list(commandArgs(trailingOnly = TRUE))
    if (length(argv) != 4) {
        cat(usage(), sep = "", file = stderr())
        quit(save = "no", 2)
    }

    setNames(argv, c("SAMPLES", "CONTIGS", "DCRS", "OUTDIR"))
}

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

cleanup <- function(e) {
    done <- get0("done", envir = e, mode = "logical", inherits = FALSE)
    if (is.null(done) || isTRUE(done)) {
        return()
    }

    outdir <- get0("outdir", envir = e, mode = "character", inherits = FALSE)
    if (is.null(outdir) || !dir.exists(outdir)) {
        return()
    }

    ap_outdir <- get0("aneuploidies_outdir", envir = e, mode = "character", inherits = FALSE)
    if (!is.null(ap_outdir) && dir.exists(ap_outdir)) {
        ap_files <- list.files(ap_outdir)
        if (length(ap_files) == 0) {
            unink(ap_outdir, recursive = TRUE)
        }
    }

    p_files <- list.files(outdir)
    if (length(p_files) == 0) {
        unlink(outdir, recursive = TRUE)
    }
}

bin_dt <- function(x, nbins = 36L) {
    if (nrow(x) < nbins) {
        x$bin <- seq_len(nrow(x))

        return(x)
    }

    q <- nrow(x) %/% nbins
    r <- nrow(x) %% nbins

    binwidths <- rep.int(q, nbins)
    binwidths[seq_len(r)] <- q + 1

    x$bin <- factor(rep(seq_len(nbins), times = binwidths),
                    levels = seq_len(nbins),
                    ordered = TRUE)

    x
}

batch_contig_dcr <- function(contig, samples, dcr_paths) {
    region <- gregion(contig, 1, TABIX_MAX_SEQLEN)
    dcr <- as.data.table(get_samples_dcr(region, dcr_paths, samples))
    setkey(dcr, chr, start, end)

    dcr
}

cohort_contig_dcr <- function(contig, batches, sample_map, dcr_path_map) {
    dcrs <- lapply(batches,
                   \(x) batch_contig_dcr(contig,
                                         gethash(sample_map, x),
                                         gethash(dcr_path_map, x)))

    cohort_dcr <- Reduce(\(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE), dcrs)
    setkey(cohort_dcr, chr, start, end)

    cohort_dcr
}

plot_cohort_contig <- function(x,
                               contig,
                               samples,
                               outdir,
                               baseline = 2,
                               del_thresh = 1.5,
                               dup_thresh = 2.5) {
    binned_dcr <- x[, lapply(.SD, \(x) mean(x, na.rm = TRUE)), .SDcols = samples, by = "bin"] |>
        melt(measure.vars = samples, variable.name = "sample", value.name = "bin_dcr")
    setkey(binned_dcr, bin, sample)
    nbins <- uniqueN(binned_dcr, by = "bin")

    png(file.path(outdir, sprintf("%s--ploidy.png", contig)),
        res = 300,
        width = 2880,
        height = 1620)
    old_par <- par(mar = c(4.1, 4.1, 1.1, 1.1), cex.lab = 1.3, cex.axis = 1.3)
    plot(NULL,
         ylim = c(0, 5),
         xlim = c(0.5, nbins + 0.5),
         type = "n",
         axes = FALSE,
         xaxs = "i",
         xlab = "",
         ylab = "")
    if (nbins > 1) {
        for (i in seq(1, nbins, by = 2)) {
            rect(i - 0.5, par("usr")[[3]],
                 i + 0.5, par("usr")[[4]],
                 border = NA, col =  "#77777722")
        }
    }
    abline(h = baseline, lwd = 2)
    boxplot(split(binned_dcr$bin_dcr, binned_dcr$bin),
            las = 2,
            add = TRUE,
            xaxt = "n",
            pars = list(outpch = NA,
                        whisklty = "solid",
                        boxfill = "white"))

    ap_samples <- binned_dcr[bin_dcr < del_thresh | bin_dcr > dup_thresh, ]
    ap_cnts <- ap_samples[, list(count = .N), by = "sample"]
    ap <- ap_samples[ap_cnts[count > 1], on = "sample"]
    if (nrow(ap) > 0) {
        grps <- split(ap, by = "sample", drop = TRUE)
        for (g in grps) {
            bd <- rep.int(NA_real_, nbins)
            bd[g$bin] <- g$bin_dcr
            lines(seq_along(bd), bd, lwd = 0.7, col = "#44444477")
        }
    }

    dels <- binned_dcr[bin_dcr < del_thresh, ]
    if (nrow(dels) > 0) {
        points(dels$bin, dels$bin_dcr, col = "#D43925", cex = 0.5, pch = 19)
    }
    dups <- binned_dcr[bin_dcr > dup_thresh, ]
    if (nrow(dups) > 0) {
        points(dups$bin, dups$bin_dcr, col = "#2376B2", cex = 0.5, pch = 19)
    }

    text(1, y = 5,
         labels = paste0("N=", prettyNum(length(samples), big.mark = ", "), " samples"),
         adj = c(0, 1),
         cex = 1.3)
    title(xlab = sprintf("%s Binned (%d bins)", contig, nbins), line = 2)
    title(ylab = "Mean Denoised Coverage", line = 2.5)
    par(old_par)
    dev.off()

}

plot_sample_contig <- function(x,
                               sample_id,
                               contig,
                               outdir,
                               baseline = 2,
                               del_thresh = 1.5,
                               dup_thresh = 2.5) {
    png(file.path(outdir, sprintf("%s--%s--ploidy.png", contig, sample_id)),
        res = 300,
        width = 2880,
        height = 1620)
    old_par <- par(mar = c(4.1, 4.1, 3.1, 1.1), cex.lab = 1.3, cex.axis = 1.3)
    plot(seq_along(x),
         x,
         type = "l",
         lwd = 1,
         ylim = c(0, 5),
         las = 2,
         xaxt = "n",
         xaxs = "i",
         xlab = "",
         ylab = "",
         main = sample_id)
    del_idx <- which(x < del_thresh)
    if (length(del_idx) > 0) {
        points(del_idx, x[del_idx], col = "#D43925", cex = 0.5, pch = 19)
    }
    dup_idx <- which(x > dup_thresh)
    if (length(dup_idx) > 0) {
        points(dup_idx, x[dup_idx], col = "#2376B2", cex = 0.5, pch = 19)
    }
    abline(h = baseline, lwd = 2)

    title(xlab = sprintf("%s Binned (%d bins)", contig, length(x)), line = 2)
    title(ylab = "Mean Denoised Coverage", line = 2.5)
    par(old_par)
    dev.off()
}


args <- parse_args()

suppressPackageStartupMessages(library(gelpers))

log_info("reading samples")
samples_dt <- fread(args$SAMPLES,
                   header = FALSE,
                   sep = "\t",
                   col.names = c("batch", "sample"),
                   key = "batch")
samples_dt <- unique(samples_dt)
if (anyDuplicated(samples_dt, by ="sample") > 0) {
    log_error("samples table contains duplicate sample IDs")
    quit(save = "no", status = 1)
}
batches <- unique(samples_dt$batch)
sample_map <- hashtab(size = length(batches))
batch_splits <- split(samples_dt$sample, samples_dt$batch)
mapply(\(x, y) sethash(sample_map, x, y),
       names(batch_splits),
       batch_splits) |>
    invisible()

log_info("reading contig list")
contigs <- readLines(args$CONTIGS)

log_info("reading dCR paths")
dcr_path_map <- read_dcr_list(args$DCRS)

outdir <- args$OUTDIR

done <- FALSE
setup_outdir(outdir)
invisible(reg.finalizer(.GlobalEnv, cleanup, onexit = TRUE))
aneuploidies_outdir <- file.path(outdir, "aneuploidies")
dir.create(aneuploidies_outdir)

ploidies <- vector(mode = "list", length = length(contigs))
names(ploidies) <- contigs
aneuploidies <- vector(mode = "list", length = length(contigs))
names(aneuploidies) <- contigs

log_info(sprintf("checking ploidy on %d contigs", length(contigs)))
for (co in contigs) {
    log_info(co)
    contig_dcr <- cohort_contig_dcr(co, batches, sample_map, dcr_path_map) |>
        bin_dt()
    samples <- setdiff(colnames(contig_dcr), c("chr", "start", "end", "bin"))

    plot_cohort_contig(contig_dcr, co, samples, outdir) |>
       invisible()

    mean_dcr <- contig_dcr[, lapply(.SD, \(x) mean(x, na.rm = TRUE)), .SDcols = samples] |>
        melt(measure.vars = samples,
             variable.name = "sample",
             value.name = co,
             variable.factor = FALSE)
    setkey(mean_dcr, sample)
    ploidies[[co]] <- mean_dcr

    contig_aneuploidies <- mean_dcr[get(co) > 2.5 | get(co) < 1.5, ]
    if (nrow(contig_aneuploidies) == 0) {
        next
    }

    aneuploidies[[co]] <- melt(contig_aneuploidies,
                               measure.vars = co,
                               variable.name = "chr",
                               value.name = "mean_dCR",
                               variable.factor = FALSE)
    for (s in contig_aneuploidies$sample) {
        s_dcr <- bin_dt(contig_dcr[, .SD, .SDcols = s], nbins = 101)
        s_dcr <- s_dcr[, lapply(.SD, \(x) mean(x, na.rm = TRUE)), by = "bin", .SDcols = s]
        invisible(plot_sample_contig(s_dcr[[s]], s, co, aneuploidies_outdir))
    }
}

ploidy_matrix <- Reduce(\(x, y) merge(x, y, by = "sample"), ploidies)

ploidy_matrix[, c(contigs) := lapply(.SD, \(x) round(x, 3)), .SDcols = contigs]
fwrite(ploidy_matrix,
       file.path(outdir, "ploidy_matrix.tsv"),
       sep = "\t",
       col.names = TRUE,
       quote = FALSE)

aneuploidies <- Filter(Negate(is.null), aneuploidies)
if (length(aneuploidies) > 0) {
    aneuploidies <- rbindlist(aneuploidies)
    aneuploidies[, mean_dCR := round(mean_dCR, 3)]
    fwrite(aneuploidies,
           file.path(aneuploidies_outdir, "aneuploidies.tsv"),
           quote = FALSE,
           sep = "\t")
} else {
    unlink(aneuploidies_outdir, recursive = TRUE)
}

done <- TRUE
