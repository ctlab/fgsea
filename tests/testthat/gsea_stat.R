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

test_that("calcGseaStatsCumulative works", {
    set.seed(42)
    stats <- rnorm(100)
    sample <- sample(seq_along(stats), 10)
    ess <- calcGseaStatCumulative(stats, selectedStats = sample, gseaParam = 1)
    for (i in seq_along(sample)) {
        expect_equal(ess[i], calcGseaStat(stats, sample[seq_len(i)]))
    }
})
