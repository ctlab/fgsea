gmt.file <- "./inst/extdata/mouse.reactome.gmt"

pathwayLines <- strsplit(readLines(gmt.file), "\t")
examplePathways <- lapply(pathwayLines, tail, -2)
names(examplePathways) <- sapply(pathwayLines, head, 1)
devtools::use_data(examplePathways)
