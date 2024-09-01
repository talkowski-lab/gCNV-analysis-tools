# Test converting a valid data.frame
gdf <- data.frame(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90))
gr <- df_to_gr(gdf)
expect_inherits(gr, "GRanges")
expect_equal(as.character(GenomicRanges::seqnames(gr)),
             c("chr1", "chr2", "chr3"))
expect_equal(GenomicRanges::start(gr),
             c(87, 101, 18))
expect_equal(GenomicRanges::end(gr),
             c(111, 201, 90))
expect_equal(as.character(GenomicRanges::strand(gr)),
             rep("*", 3))

# Test converting a valid CNV data.frame
gdf <- data.frame(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  svtype = c("DEL", "DUP", "DUP"))
gr <- df_to_gr(gdf, cnv = TRUE)
expect_inherits(gr, "GRanges")
expect_equal(as.character(GenomicRanges::seqnames(gr)),
             c("chr1", "chr2", "chr3"))
expect_equal(GenomicRanges::start(gr),
             c(87, 101, 18))
expect_equal(GenomicRanges::end(gr),
             c(111, 201, 90))
expect_equal(as.character(GenomicRanges::strand(gr)),
             c("-", "+", "+"))

# Test converting a valid data.frame with metadata columns
gdf <- data.frame(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  foo = c("a", "b", "c"))
gr <- df_to_gr(gdf)
expect_inherits(gr, "GRanges")
expect_equal(as.character(GenomicRanges::seqnames(gr)),
             c("chr1", "chr2", "chr3"))
expect_equal(GenomicRanges::start(gr),
             c(87, 101, 18))
expect_equal(GenomicRanges::end(gr),
             c(111, 201, 90))
expect_equal(GenomicRanges::mcols(gr)$foo,
             c("a", "b", "c"))

# Test converting a valid CNV data.frame with metadata columns
gdf <- data.frame(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  svtype = c("DEL", "DUP", "DUP"),
                  foo = c("a", "b", "c"))
gr <- df_to_gr(gdf, cnv = TRUE)
expect_inherits(gr, "GRanges")
expect_equal(as.character(GenomicRanges::seqnames(gr)),
             c("chr1", "chr2", "chr3"))
expect_equal(GenomicRanges::start(gr),
             c(87, 101, 18))
expect_equal(GenomicRanges::end(gr),
             c(111, 201, 90))
expect_equal(as.character(GenomicRanges::strand(gr)),
             c("-", "+", "+"))
expect_equal(GenomicRanges::mcols(gr)$foo,
             c("a", "b", "c"))

# Test converting invalid CNV data.frame
gdf <- data.frame(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90))
expect_error(df_to_gr(gdf, cnv = TRUE))

gdf <- data.frame(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  svtype = c("DUP", "DEL", "CPX"))
expect_error(df_to_gr(gdf, cnv = TRUE))

# Test converting data.frame with prohibited metadata
gdf <- data.frame(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  svtype = c("DUP", "DEL", "DEL"),
                  strand = c("+", "+", "-"),
                  seqnames = c("1", "2", "3"))
expect_warning(df_to_gr(gdf, cnv = TRUE))

# Test converting a valid data.table
gdt <- data.table(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90))
gr <- df_to_gr(gdt)
expect_inherits(gr, "GRanges")
expect_equal(as.character(GenomicRanges::seqnames(gr)),
             c("chr1", "chr2", "chr3"))
expect_equal(GenomicRanges::start(gr),
             c(87, 101, 18))
expect_equal(GenomicRanges::end(gr),
             c(111, 201, 90))
expect_equal(as.character(GenomicRanges::strand(gr)),
             rep("*", 3))

# Test converting a valid CNV data.table
gdt <- data.table(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  svtype = c("DEL", "DUP", "DUP"))
gr <- df_to_gr(gdt, cnv = TRUE)
expect_inherits(gr, "GRanges")
expect_equal(as.character(GenomicRanges::seqnames(gr)),
             c("chr1", "chr2", "chr3"))
expect_equal(GenomicRanges::start(gr),
             c(87, 101, 18))
expect_equal(GenomicRanges::end(gr),
             c(111, 201, 90))
expect_equal(as.character(GenomicRanges::strand(gr)),
             c("-", "+", "+"))

# Test converting a valid data.table with metadata columns
gdt <- data.table(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  foo = c("a", "b", "c"))
gr <- df_to_gr(gdt)
expect_inherits(gr, "GRanges")
expect_equal(as.character(GenomicRanges::seqnames(gr)),
             c("chr1", "chr2", "chr3"))
expect_equal(GenomicRanges::start(gr),
             c(87, 101, 18))
expect_equal(GenomicRanges::end(gr),
             c(111, 201, 90))
expect_equal(GenomicRanges::mcols(gr)$foo,
             c("a", "b", "c"))

# Test converting a valid CNV data.table with metadata columns
gdt <- data.table(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  svtype = c("DEL", "DUP", "DUP"),
                  foo = c("a", "b", "c"))
gr <- df_to_gr(gdt, cnv = TRUE)
expect_inherits(gr, "GRanges")
expect_equal(as.character(GenomicRanges::seqnames(gr)),
             c("chr1", "chr2", "chr3"))
expect_equal(GenomicRanges::start(gr),
             c(87, 101, 18))
expect_equal(GenomicRanges::end(gr),
             c(111, 201, 90))
expect_equal(as.character(GenomicRanges::strand(gr)),
             c("-", "+", "+"))
expect_equal(GenomicRanges::mcols(gr)$foo,
             c("a", "b", "c"))

# Test converting invalid CNV data.table
gdt <- data.table(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90))
expect_error(df_to_gr(gdt, cnv = TRUE))

gdt <- data.table(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  svtype = c("DUP", "DEL", "CPX"))
expect_error(df_to_gr(gdt, cnv = TRUE))

# Test converting data.table with prohibited metadata
gdt <- data.table(chr = c("chr1", "chr2", "chr3"),
                  start = c(87, 101, 18),
                  end = c(111, 201, 90),
                  svtype = c("DUP", "DEL", "DEL"),
                  strand = c("+", "+", "-"),
                  seqnames = c("1", "2", "3"))
expect_warning(df_to_gr(gdt, cnv = TRUE))
