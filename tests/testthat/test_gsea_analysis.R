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
