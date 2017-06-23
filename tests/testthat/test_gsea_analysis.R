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

test_that("fgseaL works", {
    skip_on_bioc()
    set.seed(42)

    skip_if_not(require("limma"))
    skip_if_not(require("GEOquery"))

    library(limma)
    library(GEOquery)
    es <- getGEO("GSE19429", AnnotGPL = TRUE)[[1]]
    exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")
    es <- es[!grepl("///", fData(es)$`Gene ID`), ]
    es <- es[fData(es)$`Gene ID` != "", ]
    es <- es[order(apply(exprs(es), 1, mean), decreasing=T), ]
    es <- es[!duplicated(fData(es)$`Gene ID`), ]
    rownames(es) <- fData(es)$`Gene ID`

    pathways <- reactomePathways(rownames(es))
    mat <- exprs(es)
    labels <- as.numeric(as.factor(gsub(" .*", "", es$title)))
    fgseaRes <- fgseaL(pathways, mat, labels, nperm = 1000, minSize = 15, maxSize = 500)
    expect_true(abs(fgseaRes[1, ES] - -0.1754936) < 1e-5)
    expect_true(fgseaRes[1, pval] > 0.5)
    expect_true(fgseaRes[1, pval] < 0.95)

})
