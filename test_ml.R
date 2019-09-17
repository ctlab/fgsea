devtools::load_all()
library(BiocParallel)
register(SerialParam())

data(exampleRanks)
data(examplePathways)

set.seed(42)
system.time({
fgseaMultilevelRes <- fgseaMultilevel(pathways = examplePathways["5991454_M_Phase"],
#fgseaMultilevelRes <- fgseaMultilevel(pathways = examplePathways["5990980_Cell_Cycle"],
#fgseaMultilevelRes <- fgseaMultilevel(pathways = examplePathways,
                                      stats = exampleRanks,
                                      minSize=15,
                                      maxSize=500,
                                      sampleSize = 101
                                      #, absEps = 1e-10
                                      )
})
head(fgseaMultilevelRes[order(pval), ])

