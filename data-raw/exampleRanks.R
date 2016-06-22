rnk.file <- "./inst/extdata/naive.vs.th1.rnk"

exampleRanks <- read.table(rnk.file,
                    header=T, colClasses = c("character", "numeric"))
exampleRanks <- setNames(exampleRanks$t, exampleRanks$ID)
devtools::use_data(exampleRanks)
