devtools::load_all()
data("examplePathways")
data("exampleRanks")

set.seed(42)

start_profiler("prof.out")

system.time(
fgseaMultilevelRes <- fgseaMultilevel(pathways = examplePathways["5991454_M_Phase"],
#fgseaMultilevelRes <- fgseaMultilevel(pathways = examplePathways,
                                      stats = exampleRanks,
                                      minSize=15,
                                      maxSize=500,
                                      sampleSize = 1001,
                                      nproc=1
                                      )
)
stop_profiler()


