gr <- GenomicRanges::GRanges(c("chr1", "chr2", "19"),
                             ranges = IRanges::IRanges(1:3, 17:19),
                             strand = c("+", "-", "*"))
GenomicRanges::mcols(gr)$metadata <- c("foo", "bar", "baz")

# Test removing the same prefix
gr_p <- rmPrefix(gr, "chr")
expect_equal(as.character(GenomicRanges::seqnames(gr_p)),
             c("1", "2", "19"))
expect_identical(GenomicRanges::ranges(gr_p), GenomicRanges::ranges(gr))
expect_identical(GenomicRanges::strand(gr_p), GenomicRanges::strand(gr))
expect_identical(GenomicRanges::mcols(gr_p), GenomicRanges::mcols(gr))

# Test removing different prefixes
gr <- GenomicRanges::GRanges(c("foochr1", "barchr2", "chr19"),
                             ranges = IRanges::IRanges(1:3, 17:19),
                             strand = c("+", "-", "*"))
GenomicRanges::mcols(gr)$metadata <- c("foo", "bar", "baz")
gr_p <- rmPrefix(gr, c("foo", "bar", "baz"))
expect_equal(as.character(GenomicRanges::seqnames(gr_p)),
             paste0(c("chr1", "chr2", "chr19")))
expect_identical(GenomicRanges::ranges(gr_p), GenomicRanges::ranges(gr))
expect_identical(GenomicRanges::strand(gr_p), GenomicRanges::strand(gr))
expect_identical(GenomicRanges::mcols(gr_p), GenomicRanges::mcols(gr))

# Test removing empty string prefix
gr_p <- rmPrefix(gr, c("foo", "bar", ""))
expect_equal(as.character(GenomicRanges::seqnames(gr_p)),
             c("chr1", "chr2", "chr19"))
expect_identical(GenomicRanges::ranges(gr_p), GenomicRanges::ranges(gr))
expect_identical(GenomicRanges::strand(gr_p), GenomicRanges::strand(gr))
expect_identical(GenomicRanges::mcols(gr_p), GenomicRanges::mcols(gr))

# Test that incorrect length of prefix vector causes error
expect_error(rmPrefix(gr, c("foo", "bar")))

# Test that NA prefix causes error
expect_error(rmPrefix(gr, NA_character_))

# Test that removing a prefix identical to a sequence name causes error
expect_error(rmPrefix(gr, "chr19"))
