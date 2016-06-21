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
})

test_that("plotEnrichment works", {
    data(examplePathways)
    data(exampleRanks)
    g <- plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
                        exampleRanks)
    tf <- tempfile("plot", fileext = ".png")
    ggsave(tf, plot=g)
})
