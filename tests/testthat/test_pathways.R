context("Pathways")

test_that("reactomePathways works", {
    if (!requireNamespace("reactome.db")) {
        skip("No reactome.db")
    }
    data(exampleRanks)
    pathways <- reactomePathways(names(exampleRanks))

    expect_true("12504" %in% pathways$`TCR signaling`)
    expect_false(any(duplicated(names(pathways))))
})

test_that("gmtPathways works", {
    pathways <- gmtPathways(system.file("extdata", "mouse.reactome.gmt", package="fgsea"))

    expect_true("12504" %in% pathways$`5991751_TCR_signaling`)
})
