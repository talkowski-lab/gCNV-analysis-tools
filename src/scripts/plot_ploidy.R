# Ploidy plotting workflow

usage <- function() {
'Usage: Rscript plot_ploidy.R DCRS BINS OUTPUT

DCRS    A file with either the paths, one per line, to the dCR matrices or
        two tab-separated columns with the batch ID in the first column and
        path to the dCR matrix for that batch in the second column
BINS    Genomic intervals file used by gCNV pipeline. The intervals can be a
        subset of the contigs.
OUTDIR  Path to the output directory. The directory must not exist
'
}

ncpus <- function() {
    n <- detectCores()
    max(n, 1, na.rm = TRUE)
}

parse_args <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    if (length(argv) != 3) {
        cat(usage(), sep = "", file = stderr())
        quit(save = "no", 2)
    }

    args <- as.list(argv)
    names(args) <- c("DCRS", "BINS", "OUTDIR")

    args
}

#' Create the output directory for the plots. Signal error if it already exists
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

bin_chr <- function(x, nbins = 36L) {
    if (nrow(x) < nbins) {
        x$bin <- seq_len(nrow(x))
        ranges <- x[, list(rstart = min(start), rend = max(end), binsize = .N), by = bin]

        return(ranges)
    }

    q <- nrow(x) %/% nbins
    r <- nrow(x) %% nbins

    binwidths <- rep.int(q, nbins)
    binwidths[seq_len(r)] <- q + 1

    bins <- rep(seq_len(nbins), times = binwidths)
    x$bin <- bins
    ranges <- x[, list(rstart = min(start), rend = max(end), binsize = .N), by = bin]

    ranges
}

all_dcr_paths <- function(h) {
    paths <- vector("list", numhash(h))
    idx <- 0
    maphash(h,
            function(k, v) {
                idx <<- idx + 1
                paths[[idx]] <<- v
            })

    paths
}

region_dcr <- function(region, dcr_paths) {
    dcr <- as.data.table(get_batch_dcr(region, dcr_paths))
    mean_dcr <- dcr[, lapply(.SD, mean, na.rm = TRUE), .SDcols = !c("chr", "start", "end")]
    mean_dcr <- melt(mean_dcr,
                   measure.vars = colnames(mean_dcr),
                   variable.name = "sample",
                   value.name = "mean_dcr",
                   variable.factor = FALSE)
    mean_dcr[, `:=`(chr = region$chr, start = region$start, end = region$end)]
    mean_dcr[is.nan(mean_dcr), mean_dcr := NA_real_]

    mean_dcr[, list(chr, start, end, sample, mean_dcr)]
}

batch_dcr <- function(bins, chr, batch) {
    apply(bins, 1, \(x) region_dcr(gregion(chr, x[[2]], x[[3]]), batch)) |>
        rbindlist()
}

chr_dcr <- function(bins, chr, batches) {
    x <- lapply(batches, \(y) batch_dcr(bins, chr, y)) |>
        rbindlist()
    setkey(x, chr, start)
    x
}

plot_chr <- function(x, chr, outdir) {
    groups <- split(x, by = c("chr", "start", "end")) |>
        lapply(\(y) y$mean_dcr)
    nsamples <- uniqueN(x, by = "sample")
    png(file.path(outdir, sprintf("%s_ploidy.png", chr)),
        res = 300,
        width = 2880,
        height = 1620)
    old_par <- par(mar = c(4.1, 4.1, 1.1, 1.1), cex.lab = 1.3, cex.axis = 1.3)
    plot(NULL,
         ylim = c(0, 5),
         xlim = c(0.5, length(groups) + 0.5),
         type = "n",
         axes = FALSE,
         xaxs = "i",
         xlab = "",
         ylab = "")
    if (length(groups) > 1) {
        for (i in seq(1, length(groups), by = 2)) {
            rect(i - 0.5, par("usr")[[3]],
                 i + 0.5, par("usr")[[4]],
                 border = NA, col =  "#77777722")
        }
    }
    lines(par("usr")[c(1,2)], c(2, 2), lwd = 1.2)
    boxplot(groups,
            las = 2,
            add = TRUE,
            xaxt = "n",
            pars = list(outcex = 0.4,
                        outpch = 19,
                        whisklty = "solid",
                        outcol = "#77777799",
                        boxfill = "white"))
    text(1, y = 5,
         labels = paste0("N=", prettyNum(nsamples, big.mark = ", "), " samples"),
         adj = c(0, 1),
         cex = 1.3)
    title(xlab = sprintf("%s Binned", chr), line = 2)
    title(ylab = "Mean Denoised Coverage", line = 2.5)
    dev.off()
    par(old_par)
}

args <- parse_args()

suppressPackageStartupMessages(library(gelpers))
suppressPackageStartupMessages(library(GenomicRanges))

log_info("reading dCR paths")
dcrs <- read_dcr_list(args$DCRS)

log_info("reading bins")
bins <- read_gcnv_bins(args$BINS, reduce = TRUE)
bins <- data.table(
    chr = as.character(seqnames(bins)),
    start = start(bins),
    end = end(bins)
)
setkey(bins, chr, start)
bins <- bins[grepl("[1-9]|1[0-9]|2[0-2]|[XY]", chr), ]
contigs <- unique(bins$chr)

plotting_bins <- lapply(split(bins, by = "chr"), bin_chr)

setup_outdir(args$OUTDIR)

dcr_paths_list <- all_dcr_paths(dcrs)

nworkers <- max(ncpus() - 1L, 1L)

ploidy <- mcmapply(
    chr_dcr,
    plotting_bins,
    names(plotting_bins),
    MoreArgs = list(batches = dcr_paths_list),
    SIMPLIFY = FALSE,
    mc.cores = nworkers
)

invisible(mcmapply(
    plot_chr,
    ploidy,
    names(ploidy),
    MoreArgs = list(outdir = args$OUTDIR),
    mc.cores = nworkers
))
