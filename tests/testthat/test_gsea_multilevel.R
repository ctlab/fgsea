context("GSEA multilevel analysis")

test_that("fgseaMultilevel works", {
	data(examplePathways)
	data(exampleRanks)
	set.seed(42)
	sampleSize <- 50
	fgseaMultilevelRes <- fgseaMultilevel(examplePathways, exampleRanks, sampleSize=sampleSize,
										  maxSize=500)
	expect_equal(fgseaMultilevelRes[23, ES], 0.5788464)

	expect_true("70385" %in% fgseaMultilevelRes[grep("5991851", pathway), leadingEdge][[1]])
	expect_true(!"68549" %in% fgseaMultilevelRes[grep("5991851", pathway), leadingEdge][[1]])

	expect_true(!"11909" %in% fgseaMultilevelRes[grep("5992314", pathway), leadingEdge][[1]])
    expect_true("69386" %in% fgseaMultilevelRes[grep("5992314", pathway), leadingEdge][[1]])

    # specifying number of threads
    fgseaMultilevelRes <- fgseaMultilevel(examplePathways, exampleRanks, sampleSize=100, maxSize=100, nproc=2)
	})


test_that("fgseaMultilevel is reproducable independent of bpparam settings", {

    data(examplePathways)
    data(exampleRanks)
    sampleSize <- 50

    set.seed(42)
    fr <- fgseaMultilevel(examplePathways[1:2], exampleRanks, sampleSize=sampleSize, maxSize=500, nproc=1)


    set.seed(42)
    fr1 <- fgseaMultilevel(examplePathways[1:2], exampleRanks, sampleSize=sampleSize, maxSize=500)
    expect_equal(fr1$nMoreExtreme, fr$nMoreExtreme)

    set.seed(42)
    fr2 <- fgseaMultilevel(examplePathways[1:2], exampleRanks, sampleSize=sampleSize, maxSize=500, nproc=0)
    expect_equal(fr2$nMoreExtreme, fr$nMoreExtreme)

    set.seed(42)
    fr3 <- fgseaMultilevel(examplePathways[1:2], exampleRanks, sampleSize=sampleSize, maxSize=500, nproc=2)
    expect_equal(fr3$nMoreExtreme, fr$nMoreExtreme)
})


test_that("fgseaMultilevel works with zero pathways", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    sampleSize <- 50
    fgseaMultilevelRes <- fgseaMultilevel(examplePathways, exampleRanks, sampleSize=sampleSize,
                      minSize=50, maxSize=10)
    expect_equal(nrow(fgseaMultilevelRes), 0)
    fgseaMultilevelRes1 <- fgseaMultilevel(examplePathways[1], exampleRanks, sampleSize=sampleSize)
    expect_equal(colnames(fgseaMultilevelRes), colnames(fgseaMultilevelRes1))
})


test_that("fgseaMultilevel throws a warning when there are duplicate gene names", {
    data(examplePathways)
    data(exampleRanks)
    exampleRanks.dupNames <- exampleRanks
    names(exampleRanks.dupNames)[41] <- names(exampleRanks.dupNames)[42]

    expect_warning(fgseaMultilevel(examplePathways, exampleRanks.dupNames, sampleSize=100, minSize=10, maxSize=50, nproc=1))

})

test_that("fgseaMultilevel: Ties detection in ranking works", {
    data(examplePathways)
    data(exampleRanks)
    exampleRanks.ties <- exampleRanks
    exampleRanks.ties[41] <- exampleRanks.ties[42]
    exampleRanks.ties.zero <- exampleRanks.ties
    exampleRanks.ties.zero[41] <- exampleRanks.ties.zero[42] <- 0

    expect_silent(fgseaMultilevel(examplePathways, exampleRanks,
                                  minSize=10, maxSize=50, nproc=1))

    expect_warning(fgseaMultilevel(examplePathways, exampleRanks.ties,
                                   minSize=10, maxSize=50, nproc=1))

    expect_silent(fgseaMultilevel(examplePathways, exampleRanks.ties.zero,
                                  minSize=10, maxSize=50, nproc=1))
})

test_that("fgseaMultilevel gives valid P-value for 5990980_Cell_Cycle", {
    data(examplePathways)
    data(exampleRanks)
    ranks <- sort(exampleRanks, decreasing = TRUE)
    example.pathway <- examplePathways["5990980_Cell_Cycle"]

    set.seed(42)
    fgseaMRes <- fgseaMultilevel(example.pathway, ranks)
    pval <- fgseaMRes$pval
    expect_true(1e-28 <= pval && pval <= 1e-26)
})

test_that("The absEps parameter works correctly with 5990980_Cell_Cycle", {
    data(examplePathways)
    data(exampleRanks)
    ranks <- sort(exampleRanks, decreasing = TRUE)
    example.pathway <- examplePathways["5990980_Cell_Cycle"]

    set.seed(42)
    fgseaMRes <- fgseaMultilevel(example.pathway, ranks, absEps = 1e-10)

    pval <- fgseaMRes$pval
    expect_true(pval == 1e-10)
    expect_true(is.na(fgseaMRes$log2err))
})

test_that("The absEps parameter works correctly with 5991504_Extension_of_Telomeres", {
    data(examplePathways)
    data(exampleRanks)
    example.pathway <- examplePathways["5991504_Extension_of_Telomeres"]

    set.seed(42)
    pvals <- replicate(fgseaMultilevel(example.pathway,
                                       exampleRanks, absEps = 1e-5)$pval,
                       n = 20)

    expect_true(all(pvals >= 1e-5))

})


test_that("fgseaMultilevel works correctly with zeros in the tail of the conditional probability vector", {
    data(exampleRanks)
    ranks <- sort(exampleRanks, decreasing = TRUE)

    expect_silent(
        fgsea:::fgseaMultilevelCpp(enrichmentScores = 0.001, ranks = ranks,
                               pathwaySize = 100, sampleSize = 3,
                               seed = 333, absEps = 0.0, sign = FALSE)
    )
})
