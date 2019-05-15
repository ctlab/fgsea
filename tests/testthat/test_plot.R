context("Plots")

test_that("plotGseaTable works", {
    data(examplePathways)
    data(exampleRanks)
    fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=1000,
                      minSize=15, maxSize=100)
    tf <- tempfile("plot", fileext = ".png")
    topPathways <- fgseaRes[head(order(pval), n=15)][order(NES), pathway]
    png(filename = tf, width=2000, height=1600, res = 300)
    plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, gseaParam=0.5)
    dev.off()
    expect_true(TRUE) # check that didn't fail before
})

test_that("plotEnrichment works", {
    data(examplePathways)
    data(exampleRanks)
    g <- plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
                        exampleRanks)
    tf <- tempfile("plot", fileext = ".png")
    ggsave(tf, plot=g)
    expect_true(TRUE) # check that didn't fail before
})

test_that("plotGseaTable ignores empty pathways", {
    data(examplePathways)
    data(exampleRanks)
    fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=1000,
                      minSize=15, maxSize=100)

    expect_equal(length(intersect(examplePathways[477], names(exampleRanks))), 0)
    p <- plotGseaTable(examplePathways[477], exampleRanks, fgseaRes, gseaParam=0.5,
                       render = FALSE)
    expect_true(is(p, "gtable"))
})
