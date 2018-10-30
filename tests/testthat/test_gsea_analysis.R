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

test_that("fgsea is reproducable independent of bpparam settings", {

    data(examplePathways)
    data(exampleRanks)
    nperm <- 2000

    set.seed(42)
    fr <- fgsea(examplePathways[1:2], exampleRanks, nperm=nperm, maxSize=500, nproc=1)


    set.seed(42)
    fr1 <- fgsea(examplePathways[1:2], exampleRanks, nperm=nperm, maxSize=500)
    expect_equal(fr1$nMoreExtreme, fr$nMoreExtreme)

    set.seed(42)
    fr2 <- fgsea(examplePathways[1:2], exampleRanks, nperm=nperm, maxSize=500, nproc=0)
    expect_equal(fr2$nMoreExtreme, fr$nMoreExtreme)

    set.seed(42)
    fr3 <- fgsea(examplePathways[1:2], exampleRanks, nperm=nperm, maxSize=500, nproc=2)
    expect_equal(fr3$nMoreExtreme, fr$nMoreExtreme)
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
    expect_true(!is.null(fgseaRes))
})

test_that("fgseaLabel example works", {
    skip_on_bioc()
    skip_on_travis() # not working on travis for some reason, check later (2018-04-10)
    example("fgseaLabel", run.donttest = TRUE, local = FALSE)
    expect_true(!is.null(fgseaRes))
})

test_that("Ties detection in ranking works", {
    data(examplePathways)
    data(exampleRanks)
    exampleRanks.ties <- exampleRanks
    exampleRanks.ties[41] <- exampleRanks.ties[42]
    exampleRanks.ties.zero <- exampleRanks.ties
    exampleRanks.ties.zero[41] <- exampleRanks.ties.zero[42] <- 0

    expect_silent(fgsea(examplePathways, exampleRanks, nperm=100, minSize=10, maxSize=50, nproc=1))

    expect_warning(fgsea(examplePathways, exampleRanks.ties, nperm=100, minSize=10, maxSize=50, nproc=1))

    expect_silent(fgsea(examplePathways, exampleRanks.ties.zero, nperm=100, minSize=10, maxSize=50, nproc=1))
})

test_that("fgsea throws a warning when there are duplicate gene names", {
    data(examplePathways)
    data(exampleRanks)
    exampleRanks.dupNames <- exampleRanks
    names(exampleRanks.dupNames)[41] <- names(exampleRanks.dupNames)[42]

    expect_warning(fgsea(examplePathways, exampleRanks.dupNames, nperm=100, minSize=10, maxSize=50, nproc=1))

})

test_that("fgsea returns leading edge ordered by decreasing of absolute statistic value", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    nperm <- 100
    fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=nperm, maxSize=50)

    expect_true(abs(exampleRanks[fgseaRes$leadingEdge[[1]][1]]) >
                abs(exampleRanks[fgseaRes$leadingEdge[[1]][2]]))

    expect_true(abs(exampleRanks[fgseaRes$leadingEdge[[2]][1]]) >
                abs(exampleRanks[fgseaRes$leadingEdge[[2]][2]]))
})

test_that("collapsePathways work", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    nperm <- 100
    pathways <- list(p1=examplePathways$`5991851_Mitotic_Prometaphase`)
    pathways <- c(pathways, list(p2=unique(c(pathways$p1, sample(names(exampleRanks), 20)))))
    pathways <- c(pathways, list(p3=sample(pathways$p1, floor(length(pathways$p1) * 0.8))))
    fgseaRes <- fgsea(pathways, exampleRanks, nperm=nperm, maxSize=500)
    collapsedPathways <- collapsePathways(fgseaRes[order(pval)],
                                          pathways, exampleRanks)
    collapsedPathways$mainPathways
    expect_identical("p1", collapsedPathways$mainPathways)
})
