context("Pathways")

test_that("reactomePathways works", {
    if (!requireNamespace("reactome.db")) {
        skip("No reactome.db")
    }
    data(exampleRanks)
    pathways <- reactomePathways(names(exampleRanks))

    expect_true("11461" %in% pathways$`Chromatin organization`)
})

test_that("gmtPathways works", {
    pathways <- gmtPathways(system.file("extdata", "mouse.reactome.gmt", package="fgsea"))

    expect_true("11461" %in% pathways$`5992314_Chromatin_organization`)
})
