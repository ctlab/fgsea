library(fgsea)
library(data.table)
library(BiocParallel)

data(exampleRanks)
data(examplePathways)

ranks <- sort(exampleRanks, decreasing = TRUE)

ranks <- abs(floor(ranks * 0.5 + 0.5))


write.table(c(length(ranks), ranks), file="roundRanks.txt", row.names = FALSE, col.names = FALSE)



fgseaRes <- fgsea(examplePathways, ranks, nperm=100,
                  minSize = 15, maxSize=500)

inpData <- fgseaRes[, .(size, ES)]


filePath <- "inpPathways.txt"

write.table(nrow(inpData), file = filePath,
            col.names = FALSE,
            row.names = FALSE)
write.table(inpData,
            file = filePath,
            append = TRUE,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t")


set.seed(42)
inpData[, `seed` := sample.int(1e9, size = nrow(inpData))]

inpData <- split(inpData, seq(nrow(inpData)))

sampleSize <- 1001
multilevelPvals <- unlist(bplapply(inpData, function(x) fgsea:::fgseaMultilevelCpp(enrichmentScores = x$ES,ranks = ranks,
                                                                            pathwaySize = x$size, sampleSize = sampleSize,
                                                                            seed = x$seed,
                                                                            eps = 0.0,
                                                                            sign = TRUE)))

write.table(multilevelPvals,
            file="multilevelResults.tsv",
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t")
