context("Gene Set Co-Regulation Analysis")

test_that("GESECA: is reproducible", {
    data("exampleExpressionMatrix")
    data("examplePathways")

    set.seed(1)
    gr1 <- geseca(pathways=examplePathways, E=exampleExpressionMatrix)

    set.seed(1)
    gr2 <- geseca(pathways=examplePathways, E=exampleExpressionMatrix)

    expect_equal(gr1$pval, gr2$pval)
})


test_that("GESECA: works with zero pathways", {
    data("exampleExpressionMatrix")
    data("examplePathways")
    set.seed(42)
    sampleSize <- 11

    gr1 <- geseca(exampleExpressionMatrix, examplePathways,
                  sampleSize=sampleSize, minSize=50, maxSize=10)
    expect_equal(nrow(gr1), 0)

    gr2 <- geseca(pathways=examplePathways[5],
                  E=exampleRanks, sampleSize=sampleSize, minSize=5)

    expect_equal(colnames(gr2), colnames(gr1))
})

test_that("GESECA: throws a warning when there are duplicate gene names", {
    data("exampleExpressionMatrix")
    data("examplePathways")
    E <- exampleExpressionMatrix
    rownames(E)[1] <- rownames(E)[2]

    expect_warning(geseca(E, examplePathways, minSize=15))
})


test_that("GESECA: The eps parameter works correctly", {
    data("exampleExpressionMatrix")
    data("examplePathways")

    set.seed(42)
    expect_warning(gr <- geseca(E=exampleExpressionMatrix, pathways = examplePathways, eps = 1e-10))

    expect_true(all(gr$pval >= 1e-10))
    expect_true(any(is.na(gr$log2err)))
})



test_that("GESECA: throws a warning when sampleSize is less than 3", {
    data("exampleExpressionMatrix")
    data("examplePathways")

    set.seed(42)
    expect_silent(geseca(pathways   = examplePathways,
                         E          = exampleExpressionMatrix,
                         sampleSize = 5,
                         eps        = 0.0))
    expect_warning(geseca(pathways   = examplePathways,
                          E          = exampleExpressionMatrix,
                          sampleSize = 1,
                          eps        = 0.0))

})

test_that("GESECA: throws a warning when reaching eps", {
    data("exampleExpressionMatrix")
    data("examplePathways")

    set.seed(42)
    expect_warning(gr <- geseca(E=exampleExpressionMatrix, pathways = examplePathways, eps = 1e-10))
})
