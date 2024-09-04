# dCR plotting pipeline
#
# Generate dCR plots for possible de novo CNVs or a group of clustered CNVs.
# Usage: Rscript plot_dcr.R CALLSET [DENOVO] BINS PED DCRS OUTDIR
# * CALLSET - If DENOVO is given, the full callset produced by the gCNV
#             pipeline. Otherwise, the subset of the callset for the variants
#             to plot
# * DENOVO  - Optional. de novo calls to plot
# * BINS    - Genomic bins file used by the gCNV pipeline
# * PED     - Pedigree file
# * DCRS    - List of paths, one per line, to the dCR matrices
# * OUTDIR  - Output directory. Must not exist

# Functions -------------------------------------------------------------------
ALPHANUM <- c(LETTERS, letters, as.character(0:9))

# Parse the command line arguments
parse_args <- function(argv) {
    argc <- length(argv)
    if (argc == 5) {
        args <- vector(mode = "list", 6)
        args[[1]] <- argv[[1]]
        args[3:6] <- argv[2:5]
    } else {
        args <- as.list(argv)
    }
    names(args) <- c(
        "callset", "denovo", "bins", "pedigree", "dcrs", "outdir"
    )

    args
}

# Generate n random alphanumeric strings of length k
random_str <- function(n, k = 7) {
    replicate(n,
              paste0(sample(ALPHANUM, k, replace = TRUE), collapse = ""),
              simplify = TRUE)
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

# From a data.table of CNV calls, return the largest variant as a gregion
# object
get_largest_variant <- function(x) {
    m <- x[which.max(end - start)]

    gregion(m$chr, m$start, m$end)
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
    trunc((end - start + 1) / 1000)
}

# Get the plotting colors for CNV types
svtype_color <- function(x) {
  if (is.na(x)) {
      return("#000000")
  } else if (x == "DEL") {
      return("#D43925")
  } else if (x == "DUP") {
      return("#2376B2")
  }

  "#000000"
}

# Take a data.table of CNV calls, the dCR set corresponding to those calls, a
# output path and a plot title, write the dCR plot to the output
plot_dcr <- function(x, dcr, outfile, main = "") {
    blocks <- split_dcr_blocks(dcr$coords$in_flank)
    png(filename = paste0(outfile, ".png"), width = 647, height = 400)
    old_par <- par(font.lab = 2, mar = c(4, 4, 4, 2) + 0.1)
    matplot(
        data.matrix(dcr$bg_dcr),
        type = "l",
        lty = 1,
        ylim = c(0, 5),
        col = "#77777777",
        xlab = "Interval",
        ylab = "Denoised Coverage",
        main = main
    )
    for (i in seq_len(nrow(x))) {
        sample_id <- x[i, ]$sample
        phen <- x[i, ]$phenotype
        svtype <- x[i, ]$svtype
        lines(
            dcr$sample_dcr[[sample_id]],
            lwd = 3,
            lty = if (!is.na(phen) && phen == 2) "solid" else "dashed",
            col = svtype_color(svtype)
        )
    }
    if (length(blocks$left_flank) > 0) {
        rect(1, 0, blocks$left_flank[[2]] + 0.5, 5, col = "#33333311", border = NA)
        text(
            1, 0.3,
            labels = paste0(
                blocks$left_flank[[2]] - blocks$left_flank[[1]] + 1,
                " intervals: ",
                range2kb(
                    dcr$coords[blocks$left_flank[[1]], "start"],
                    dcr$coords[blocks$left_flank[[2]], "end"]
                ),
                "kb"
            ),
            pos = 4,
            offset = 0.3,
            cex = 0.9
        )
    }
    if (length(blocks$right_flank) > 0) {
        rect(
            blocks$right_flank[[1]] - 0.5, 0,
            blocks$right_flank[[2]], 5,
            col = "#33333311",
            border = NA
        )
        text(
            blocks$right_flank[[2]], 0.3,
            labels = paste0(
            blocks$right_flank[[2]] - blocks$right_flank[[1]] + 1,
                " intervals: ",
                range2kb(
                    dcr$coords[blocks$right_flank[[1]], "start"],
                    dcr$coords[blocks$right_flank[[2]], "end"]
                ),
                "kb"
            ),
            pos = 2,
            offset = 0.3,
            cex = 0.9
        )
    }
    legend(
        "topleft",
        legend = c("DUP", "DEL", "Unknown"),
        col = c("#2376B2", "#D43925", "#000000"),
        pch = 19,
        bty = "n"
    )
    legend(
        "topright",
        legend = c("Affected", "Unaffected"),
        lty = c("solid", "dashed"),
        bty = "n"
    )
    par(old_par)
    dev.off()
}

# Get the dCR matrix for a group of samples. Samples that are missing are
# ignored.
get_group_dcr <- function(samples,
                          batches,
                          dcrs,
                          region,
                          bins,
                          bg = 50,
                          pad = 0.2) {
    bg_per_batch <- ceiling(bg / length(unique(batches)))
    expanded_region <- expand_region_by_bins(region, bins, pad)

    coord_cols <- c("chr", "start", "end")
    batch_groups <- split(samples, batches)
    dcr_list <- vector(mode = "list", length = length(batch_groups))
    bg_dcr_list <- vector(mode = "list", length = length(batch_groups))
    for (i in seq_along(batch_groups)) {
        batch <- names(batch_groups)[[i]]
        dcr <- tryCatch(
            get_samples_dcr(expanded_region,
                            gethash(dcrs, batch),
                            batch_groups[[i]],
                            include_bg = TRUE,
                            squeeze = TRUE,
                            reduce = TRUE),
            error = function(cnd) NULL
        )

        if (is.null(dcr) || nrow(dcr) == 0) {
            next
        }

        dcr_list[[i]] <- as.data.table(dcr[c(coord_cols, batch_groups[[i]])])
        setkey(dcr_list[[i]], chr, start, end)

        bg_samples <- colnames(dcr)[!colnames(dcr) %in% c(coord_cols, samples)]
        if (length(bg_samples) > 0) {
            bg_samples <- sample(bg_samples, min(length(bg_samples), bg_per_batch))
            bg_dcr <- dcr[c(coord_cols, bg_samples)]
            # It is possible that a single sample is in different batches which
            # creates the problem of duplicate column names when the background
            # dCR data.frames are joined so we change the column names to be
            # unique as we don't need the sample IDs anyways.
            colnames(bg_dcr) <- c(
                coord_cols, paste0(bg_samples, "_", i, random_str(length(bg_samples)))
            )
            bg_dcr_list[[i]] <- as.data.table(bg_dcr)
            setkey(bg_dcr_list[[i]], chr, start, end)
        }
    }

    dcr_list <- Filter(Negate(is.null), dcr_list)
    bg_dcr_list <- Filter(Negate(is.null), bg_dcr_list)

    if (length(dcr_list) == 0) {
        return(NULL)
    }

    merged_dcr <- Reduce(merge, c(dcr_list, bg_dcr_list))
    setkey(merged_dcr, chr, start)

    if (nrow(merged_dcr) == 0) {
        return(NULL)
    }

    # Assume ranges are disjoint
    merged_dcr[, in_flank := end < region$start | start > region$end]

    bg_samples <- colnames(merged_dcr)[!colnames(merged_dcr) %in% c(coord_cols, samples, "in_flank")]
    fg_samples <- samples[samples %in% colnames(merged_dcr)]
    list(
        coords = merged_dcr[, .SD, .SDcols = c(coord_cols, "in_flank")],
        sample_dcr = merged_dcr[, .SD, .SDcols = fg_samples],
        bg_dcr = merged_dcr[, .SD, .SDcols = bg_samples]
    )
}

plot_denovo_group <- function(x, dcrs, bins, outdir) {
    calls <- unique(x, by = "sample")
    region <- get_largest_variant(calls)
    dcr <- get_group_dcr(calls$sample, calls$batch, dcrs, region, bins)
    if (is.null(dcr)) {
        log_warn(paste0("could not plot de novo '",
                        x[1, ]$variant_name,
                        "' from family '",
                        x[1, ]$family_id, "'"))
        return()
    }

    outfile <- file.path(
        outdir,
        paste(x[1, ]$family_id, x[1, ]$variant_name, sep = "__")
    )
    main <- paste0(
        to_string(region), "\n",
        "Family: ", x[1, ]$family_id, " Variant: ", x[1, ]$variant_name
    )

    plot_dcr(x, dcr, outfile, main = main)
}

plot_variant_group <- function(x, dcrs, bins, outdir) {
  region <- get_largest_variant(x)
  dcr <- get_group_dcr(x$sample, x$batch, dcrs, region, bins)
  if (is.null(dcr)) {
      log_warn(paste0("Could not plot variant ", x[1, ]$variant_name))
      return()
  }

  outfile <- file.path(
      outdir,
      paste(x[1, ]$variant_name, region$chr, region$start, region$end, sep = "__")
  )
  main <- paste0(region$chr, ":", region$start, "-", region$end, "\n", x[1, ]$variant_name)

  plot_dcr(x, dcr, outfile, main = main)
}

# Parse command-line arguments ------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 5 && length(argv) != 6) {
    stop("Incorrect number of arguments to script", call. = FALSE)
}
args <- parse_args(argv)

suppressPackageStartupMessages(library(gelpers))

set.seed(42)

callset <- read_callset(args$callset)
pedigree <- read_pedigree(args$pedigree)
dcrs <- read_dcr_list(args$dcrs)
is_hg19 <- any(c(as.character(1:22), "X", "Y") %in% callset$chr)
bins <- read_gcnv_bins(args$bins, reduce = is_hg19)
setup_outdir(args$outdir)

# Run de novo plotting workflow ----------------------------------------------
if (!is.null(args$denovo)) {
    log_info("Running de novo plotting workflow")
    batches <- unique(callset[, c("sample", "batch")])
    callset <- callset[, list(chr, start, end, svtype, variant_name, sample)]
    denovo <- read_callset(args$denovo)
    denovo[, family_id := as.character(family_id)]
    setnames(denovo, "sample", "offspring")
    denovo <- pedigree[, c("family_id", "sample_id", "phenotype")][denovo, on = "family_id"]
    denovo <- denovo[, list(family_id, offspring, sample_id, phenotype, svtype, variant_name)]
    denovo <- callset[denovo, on = c(sample = "sample_id", "svtype", "variant_name"), mult = "first"]
    denovo[is.na(chr), svtype := NA]
    denovo <- batches[denovo, on = "sample", mult = "first"]
    denovo <- denovo[!is.na(batch)]

    # Sometimes there are families in which multiple offspring possess a denovo CNV
    # and causes there to be fewer family groups than denovo calls.
    grps <- split(denovo, by = c("family_id", "variant_name"))
    for (g in grps) {
        log_info(paste0("Plotting de novo '",
                        g[1, ]$variant_name,
                        "' for family '",
                        g[1, ]$family_id,
                        "'"))
        plot_denovo_group(g, dcrs, bins, args$outdir)
    }
    log_info("Completed de novo CNV plotting")
    quit(save = "no")
}

# Run variant group plotting workflow ----------------------------------------
log_info("Running variant plotting workflow")
pedigree <- pedigree[, list(sample_id, phenotype)]
setkey(pedigree, sample_id)

setkey(callset, sample)
callset <- pedigree[callset, on = c(sample_id = "sample"), mult = "first"]
setnames(callset, "sample_id", "sample")
setkey(callset, sample, variant_name)
callset <- unique(callset, by = c("sample", "variant_name"))

# Split callset by variant ID
var_grps <- split(callset, by = "variant_name")
for (vg in var_grps) {
    log_info(paste0("Plotting variant ", vg[1, ]$variant_name))
    plot_variant_group(vg, dcrs, bins, args$outdir)
}
log_info("Completed variant group plotting")
