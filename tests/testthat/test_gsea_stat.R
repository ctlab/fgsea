context("GSEA stat calculation")

test_that("calcGseaStats works", {
    stats <- -10:10

    sample <- 1:5
    expect_equal(calcGseaStat(stats, selectedStats = sample), 1)
    sample <- 16:21
    expect_equal(calcGseaStat(stats, selectedStats = sample), -1)

    sample <- (1:5)*2
    expect_equal(calcGseaStat(stats, selectedStats = sample), 0.71)
})

test_that("calcGseaStats(returnLeadingEdge=TRUE) works", {
    stats <- -10:10

    sample <- c(10, 1:5)
    gseaRes <- calcGseaStat(stats, selectedStats = sample,
                            returnLeadingEdge = TRUE,
                            returnAllExtremes = TRUE)
    expect_equal(gseaRes$res, calcGseaStat(stats, selectedStats = sample))
    expect_true(1 %in% gseaRes$leadingEdge)
    expect_false(10 %in% gseaRes$leadingEdge)

    sample <- c(10, 16:21)
    gseaRes <- calcGseaStat(stats, selectedStats = sample,
                            returnLeadingEdge = TRUE,
                            returnAllExtremes = TRUE)
    expect_equal(gseaRes$res, calcGseaStat(stats, selectedStats = sample))
    expect_true(21 %in% gseaRes$leadingEdge)
    expect_false(10 %in% gseaRes$leadingEdge)

    sample <- 10:12
    gseaRes <- calcGseaStat(stats, selectedStats = sample,
                            returnLeadingEdge = TRUE,
                            returnAllExtremes = TRUE)
    expect_equal(gseaRes$res, calcGseaStat(stats, selectedStats = sample))
    expect_false(12 %in% gseaRes$leadingEdge)
    expect_false(10 %in% gseaRes$leadingEdge)
})

#test_that("calcGseaStats returns zero when both sides are equally enriched", {
#    stats <- -10:10
#    sample <- 10:12
#    expect_equal(calcGseaStat(stats, selectedStats = sample), 0)
#
#    ess <- calcGseaStatCumulative(stats, selectedStats = sample, gseaParam = 1)
#    for (i in seq_along(sample)) {
#        expect_equal(ess[i], calcGseaStat(stats, sample[seq_len(i)]))
#    }
#})
#
#test_that("calcGseaStats* work with zero gene-level stat", {
#
#    expect_equal(
#        calcGseaStat(c(10:1, 0, -1:-20), selectedStats = 11, 1),
#        calcGseaStat(c(10:1, -1e-9, -1:-20), selectedStats = 11, 1))
#
#    expect_equal(
#        calcGseaStat(c(10:1, 0, 0, -1:-29), selectedStats = 11:12, 1),
#        calcGseaStat(c(10:1, 0.1, -0.2, -1:-29), selectedStats = 11:12, 1))
#
#    stats <- c(10:1, 0, 0, -1:-29)
#    sample <- 11:13
#
#    ess <- calcGseaStatCumulative(stats, selectedStats = sample, gseaParam = 1)
#    for (i in seq_along(sample)) {
#        expect_equal(ess[i], calcGseaStat(stats, sample[seq_len(i)]))
#    }
#})
#
#test_that("calcGseaStatsCumulative works", {
#    set.seed(42)
#    stats <- rnorm(100)
#    sample <- sample(seq_along(stats), 10)
#    ess <- calcGseaStatCumulative(stats, selectedStats = sample, gseaParam = 1)
#    for (i in seq_along(sample)) {
#        expect_equal(ess[i], calcGseaStat(stats, sample[seq_len(i)]))
#    }
#})

