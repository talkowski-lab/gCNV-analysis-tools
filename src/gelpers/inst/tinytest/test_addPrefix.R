gr <- GenomicRanges::GRanges(c("chr1", "chr2", "chr19"),
                             ranges = IRanges::IRanges(1:3, 17:19),
                             strand = c("+", "-", "*"))
GenomicRanges::mcols(gr)$metadata <- c("foo", "bar", "baz")

# Test adding the same prefix
gr_p <- addPrefix(gr, "prefix")
expect_equal(as.character(GenomicRanges::seqnames(gr_p)),
             paste0("prefix", c("chr1", "chr2", "chr19")))
expect_identical(GenomicRanges::ranges(gr_p), GenomicRanges::ranges(gr))
expect_identical(GenomicRanges::strand(gr_p), GenomicRanges::strand(gr))
expect_identical(GenomicRanges::mcols(gr_p), GenomicRanges::mcols(gr))

# Test adding different prefixes
gr_p <- addPrefix(gr, c("bi", "bing", "bingo"))
expect_equal(as.character(GenomicRanges::seqnames(gr_p)),
             paste0(c("bi", "bing", "bingo"), c("chr1", "chr2", "chr19")))
expect_identical(GenomicRanges::ranges(gr_p), GenomicRanges::ranges(gr))
expect_identical(GenomicRanges::strand(gr_p), GenomicRanges::strand(gr))
expect_identical(GenomicRanges::mcols(gr_p), GenomicRanges::mcols(gr))

# Test adding empty string prefix
gr_p <- addPrefix(gr, c("foo", "bar", ""))
expect_equal(as.character(GenomicRanges::seqnames(gr_p)),
             paste0(c("foo", "bar", ""), c("chr1", "chr2", "chr19")))
expect_identical(GenomicRanges::ranges(gr_p), GenomicRanges::ranges(gr))
expect_identical(GenomicRanges::strand(gr_p), GenomicRanges::strand(gr))
expect_identical(GenomicRanges::mcols(gr_p), GenomicRanges::mcols(gr))

# Test that incorrect length of prefix vector causes error
expect_error(addPrefix(gr, c("foo", "bar")))
