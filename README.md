# fgsea
An R-package for preranked gene set enrichment analysis (GSEA).

See https://github.com/ctlab/fgsea-paper/ for more details on the algorithm.

## Installation

```{r}
library(devtools)
install_github("ctlab/fgsea")
```

## Quick run

Loading libraries

```{r}
library(data.table)
library(fgsea)
```

Files with example data:

```{r}
rnk.file <- system.file("naive.vs.th1.rnk", package="fgsea")
gmt.file <- system.file("mouse.reactome.gmt", package="fgsea")
```

Loading ranks.

```{r}
ranks <- read.table(rnk.file,
                    header=T, colClasses = c("character", "numeric"))
ranks <- structure(ranks$t, names=ranks$ID)
str(ranks)
```

Loading pathways.

```{r}
pathwayLines <- strsplit(readLines(gmt.file), "\t")
pathways <- lapply(pathwayLines, tail, -2)
names(pathways) <- sapply(pathwayLines, head, 1)
str(head(pathways))
```

Running fgsea (should take about 10 seconds):
```{r}
fgseaRes <- fgsea(pathways = pathways, 
                  stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000,
                  nproc=1)
```
