# Parse the command line arguments.
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
        "callset", "denovo", "pedigree", "dcr_list", "bins", "outdir"
    )

    args
}

# Create the output directory for the plots. Signal error if it already exists.
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
# object.
get_largest_variant <- function(x) {
    r <- as.list(x[which.max(end - start + 1), list(chr, start, end)])

    gregion(r$chr, r$start, r$end)
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

# Convert a genomic range to kilobases.
range2kb <- function(start, end) {
  out <- trunc((end - start + 1) / 1000)
  return(out)
}

# Get the plotting colors to CNV types.
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
# output path and a plot title, write the dCR plot to the output.
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
    sample_id <- x[i, sample]
    phen <- x[i, phenotype]
    svtype <- x[i, svtype]
    lines(
      dcr$sample_dcr[, sample_id, with = FALSE],
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
          dcr$coords[blocks$left_flank[[1]], start],
          dcr$coords[blocks$left_flank[[2]], end]
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

# Get the dCR matrix for a group of samples. If this function fails to get the
# dCR matrix, it will either signal an error or return NULL.
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
        dcr <- get_samples_dcr(
            gethash(dcrs, batch),
            expanded_region,
            batch_groups[[i]],
            include_bg = TRUE
        )
        cols <- c(coord_cols, batch_groups[[i]])
        dcr_list[[i]] <- dcr[, ..cols]

        bg_samples <- colnames(dcr)[!colnames(dcr) %in% cols]
        if (length(bg_samples) > 0) {
            bg_samples <- sample(bg_samples, min(length(bg_samples), bg_per_batch))
            bg_cols <- c(coord_cols, bg_samples)
            bg_dcr_list[[i]] <- dcr[, ..bg_cols]
        }
    }

    shared_coords <- Reduce(
        \(x, y) x[y, ..coord_cols, on = coord_cols, nomatch = NULL],
        c(dcr_list, bg_dcr_list)
    )
    if (nrow(shared_coords) == 0) {
        warning("No shared dCR ranges among samples", call. = FALSE, immediate. = TRUE)
        return(NULL)
    }

    dcr_subset <- lapply(
        dcr_list,
        \(x) x[shared_coords, !..coord_cols, on = coord_cols, mult = "first"]
    )
    sample_dcr <- do.call(cbind, dcr_subset)
    bg_dcr_subset <- lapply(
        bg_dcr_list,
        \(x) x[shared_coords, !..coord_cols, on = coord_cols, mult = "first"]
    )
    bg_dcr <- do.call(cbind, bg_dcr_subset)

    shared_coords[, in_flank := end < region$start | start > region$end]

    list(
        coords = shared_coords,
        sample_dcr = sample_dcr,
        bg_dcr = bg_dcr
    )
}

plot_denovo_group <- function(x, dcrs, bins, outdir) {
    calls <- unique(x, by = "sample")
    region <- get_largest_variant(calls)
    dcr <- tryCatch(
        get_group_dcr(calls$sample, calls$batch, dcrs, region, bins),
        error = function(cnd) NULL
    )
    if (is.null(dcr)) {
        warning(
            paste0("could not plot de novo ", calls[1, family_id], " ", calls[1, variant_name]),
            immediate. = TRUE,
            call. = FALSE
        )
    }

    outfile <- file.path(
        outdir,
        paste(df[1, family_id], df[1, variant_name], sep = "__")
    )
    main <- paste0(
        to_string(region), "\n",
        "Family: ", x[1, family_id], " Variant: ", x[1, variant_name]
    )

    plot_dcr(x, dcr, outfile, main = main)
}

plot_variant_group <- function(x, dcrs, bins, outdir) {
  region <- get_largest_variant(x)
  dcr <- tryCatch(
    get_group_dcr(x$sample, x$batch, dcrs, region, bins),
    error = function(cnd) NULL
  )
  if (is.null(dcr)) {
    warning(
      paste0("Could not plot variant ", x[1, variant_name]),
      immediate. = TRUE,
      call. = FALSE
    )
  }

  outfile <- file.path(
    outdir,
    paste(x[1, variant_name], chr, start, end, sep = "__")
  )
  main <- paste0(chr, ":", start, "-", end, "\n", x[1, variant_name])

  plot_dcr(x, dcr, outfile, main = main)
}

# Parse command-line arguments ------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 5 || length(argv) != 6) {
    stop("Incorrect number of arguments to script", call. = FALSE)
}
args <- parse_args(argv)

library(data.table)
library(gelpers)

set.seed(42)

callset <- read_callset(args$callset)
pedigree <- read_pedigree(args$pedigree)
dcrs <- read_dcr_list(args$dcr_list)
bins <- read_gcnv_bins(args$bins)
setup_outdir(args$outdir)

if (!is.null(argv$denovo)) {
    message("Running de novo plotting workflow")
    denovo <- read_callset(argv$callset)
    denovo <- pedigree[, c("sample_id", "family_id")][denovo, on = c(sample_id = "sample")]
    setnames(denovo, "sample_id", "offspring")
    denovo <- pedigree[
        denovo,
        list(family_id, offspring, sample_id, phenotype, svtype, variant_name),
        on = "family_id"]
    denovo <- merge(
        callset[, list(chr, start, end, svtype, variant_name, sample)],
        denovo,
        by.x = c("sample", "svtype", "variant_name"),
        by.y = c("sample_id", "svtype", "variant_name"),
        all.y = TRUE
    )
    denovo[is.na(chr), svtype := NA_character_]
    batches <- unique(callset[, list(sample, batch)])
    denovo <- batches[denovo, on = "sample", mult = "first"][!is.na(batch)]
  
    # Sometimes there are families in which multiple offspring possess a denovo CNV
    # and causes there to be fewer family groups than denovo calls.
    denovo_groups <- split(denovo, by = c("family_id", "variant_name"))
    for (dg in denovo_groups) {
        message(
            paste0("Plotting de novo ", dg[1, variant_name], " for family ", dg[1, family_id])
        )
        plot_denovo_group(dg, dcr_paths, intervals, args$outdir)

    }
    quit(save = "no")
}

message("Running variant plotting workflow")
callset <- pedigree[, list(sample_id, phenotype, family_id)][
  callset, on = c(sample_id = "sample"), mult = "first"]
setnames(callset, "sample_id", "sample")

# Split callset by variant ID
variant_groups <- split(callset, callset$variant_name)
for (vg in variant_groups) {
  message(
    paste0("Plotting variant ", dg[1, variant_name])
  )
  plot_variant_group(vg, dcr_file_paths, intervals, args$outdir)
}
