# Test function handles 0-length objects
gr0 <- GenomicRanges::GRanges()
expect_equal(cnvCoverage(gr0, gr0), double(0))

# Test warning is emitted when there is a range with no strand
gr1 <- GenomicRanges::GRanges(c("chr1", "chr2"),
                              IRanges::IRanges(c(1, 13), c(7, 27)),
                              c("+", "*"))
expect_warning(cnvCoverage(gr0, gr1))

# Test when every range in x overlaps a range in y
gr2 <- GenomicRanges::GRanges(c("chr1", "chr2", "chr3"),
                              IRanges::IRanges(c(1, 1, 1), c(10, 10, 10)),
                              c("+", "+", "-"))
gr3 <- GenomicRanges::GRanges(c("chr1", "chr2", "chr3"),
                              IRanges::IRanges(c(6, 6, 6), c(10, 10, 10)),
                              c("+", "+", "-"))
expect_equal(cnvCoverage(gr2, gr3), c(0.5, 0.5, 0.5))

# Test that overlap computation respects strand
GenomicRanges::strand(gr3)[3] <- "+"
expect_equal(cnvCoverage(gr2, gr3), c(0.5, 0.5, 0))

# Test when range in x overlaps multiple ranges in y
gr4 <- GenomicRanges::GRanges(c("chr1", "chr1", "chr3"),
                              IRanges::IRanges(c(1, 7, 6), c(4, 10, 10)),
                              c("+", "+", "-"))
expect_equal(cnvCoverage(gr2, gr4), c(0.8, 0, 0.5))

# Test when there is no overlap between x and y
gr5 <- GenomicRanges::GRanges(c("chr1", "chr2", "chr3"),
                              IRanges::IRanges(c(50, 50, 50), c(70, 70, 70)),
                              c("+", "+", "-"))
expect_equal(cnvCoverage(gr2, gr5), rep(0, 3))
