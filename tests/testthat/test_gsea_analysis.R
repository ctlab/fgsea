context("GSEA analysis")

test_that("fgsea works", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    nperm <- 100
    fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=nperm, maxSize=500)


    expect_equal(fgseaRes[23, ES], 0.5788464)
    expect_equal(fgseaRes[23, nMoreExtreme], 0)
    expect_gt(fgseaRes[1237, nMoreExtreme], 50 * nperm / 1000)

    expect_true("70385" %in% fgseaRes[grep("5991851", pathway), leadingEdge][[1]])
    expect_true(!"68549" %in% fgseaRes[grep("5991851", pathway), leadingEdge][[1]])

    expect_true(!"11909" %in% fgseaRes[grep("5992314", pathway), leadingEdge][[1]])
    expect_true("69386" %in% fgseaRes[grep("5992314", pathway), leadingEdge][[1]])

    # analyzing one pathway is done in a different way
    fgsea1Res <- fgsea(examplePathways[1237], exampleRanks, nperm=nperm, maxSize=500)
    expect_gt(fgseaRes[1, nMoreExtreme], 50 * nperm / 1000)

    # specifying number of threads
    fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=2000, maxSize=100, nproc=2)

    # all nMoreExtreme being even is a sign of invalid parallelization
    expect_false(all(fgseaRes$nMoreExtreme %% 2 == 0))
})

test_that("fgsea works with zero pathways", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    nperm <- 100
    fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=nperm,
                      minSize=50, maxSize=10)
    expect_equal(nrow(fgseaRes), 0)
    fgseaRes1 <- fgsea(examplePathways[1], exampleRanks, nperm=nperm)
    expect_equal(colnames(fgseaRes), colnames(fgseaRes1))
})

test_that("fgseaLabel works", {
    mat <- matrix(rnorm(1000*20),1000,20)
    labels <- rep(1:2, c(10, 10))

    rownames(mat) <- as.character(seq_len(nrow(mat)))
    pathways <- list(sample(rownames(mat), 100),
                     sample(rownames(mat), 200))

    fgseaRes <- fgseaLabel(pathways, mat, labels, nperm = 1000, minSize = 15, maxSize = 500)

})
