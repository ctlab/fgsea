context("fora analysis")

test_that("fora works if pathways consists of one element", {
    pthw <- c("1", "2", "10")
    genes <- paste(1:5)
    universe <- paste(1:10)


    res <- fora(list("test_pathway" = pthw), genes, universe)
    expect_true(nrow(res) == 1)
    expect_equal(res$pval, phyper(1, 3, 7, 5, lower.tail = FALSE))
})


test_that("fora works if pathways consist of multiple elements", {
    pthw1 <- c("1", "2", "10")
    pthw2 <- paste(1:4)
    pthw3 <- paste(6:8)

    genes <- paste(1:5)
    universe <- paste(1:10)


    res <- fora(list("test_pathway_1" = pthw1,
                     "test_pathway_2" = pthw2,
                     "test_pathway_3" = pthw3),
                genes, universe)

    expect_true(nrow(res) == 3)

    expect_equal(res$pval, sort(c(phyper(1, 3, 7, 5, lower.tail = FALSE),
                                  phyper(3, 4, 6, 5, lower.tail = FALSE),
                                  phyper(-1, 3, 7, 5, lower.tail = FALSE))))
})

test_that("fora correctly handles the case when pathways don't intersect with universe", {
    pthw1 <- c("11", "12", "13")
    genes <- paste(1:5)
    universe <- paste(1:10)

    res1 <- fora(list("p1" = pthw1), genes, universe)

    expect_true(nrow(res1) == 0)

    pthw2 <- paste(1:4)
    pthw3 <- paste(6:8)

    res2 <- fora(list("p1" = pthw1,
                      "p2" = pthw2,
                      "p3" = pthw3), genes, universe)
    expect_true(nrow(res2) == 2)
})
