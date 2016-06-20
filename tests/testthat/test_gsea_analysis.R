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

    # analyzing one pathway is done in a different way
    fgsea1Res <- fgsea(examplePathways[1237], exampleRanks, nperm=nperm, maxSize=500)
    expect_gt(fgseaRes[1, nMoreExtreme], 50 * nperm / 1000)
})

test_that("fgseaReactome works", {
    if (!requireNamespace("reactome.db")) {
        skip("No reactome.db")
    }
    data(exampleRanks)
    set.seed(42)
    fgseaRes <- fgseaReactome(exampleRanks, nperm=100, minSize=15, maxSize=100)
    expect_equal(fgseaRes[match("5990981", pathwayId), nMoreExtreme], 0)
    expect_gt(fgseaRes[match("5991238", pathwayId), nMoreExtreme], 5)
})
