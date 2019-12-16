[![Travis-CI Build Status](https://travis-ci.org/ctlab/fgsea.svg?branch=master)](https://travis-ci.org/ctlab/fgsea)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ctlab/fgsea?branch=master&svg=true)](https://ci.appveyor.com/project/assaron/fgsea)
[![codecov](https://codecov.io/gh/ctlab/fgsea/branch/master/graph/badge.svg)](https://codecov.io/gh/ctlab/fgsea)


# fgsea
An R-package for fast preranked gene set enrichment analysis (GSEA). The package 
implements a special algorithm to calculate the empirical enrichment score null distributions simultaneously
for all the gene set sizes, which allows up to **several hundred times faster** execution time compared to original
Broad implementation. See [the preprint](http://biorxiv.org/content/early/2016/06/20/060012) for algorithmic details.

Full vignette can be found here: http://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

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
library(ggplot2)
```

Loading example pathways and gene-level statistics:
```{r}
data(examplePathways)
data(exampleRanks)
```

Running fgsea (should take about 10 seconds):
```{r}
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats = exampleRanks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
```

The head of resulting table sorted by p-value:
```
                                    pathway         pval        padj        ES      NES nMoreExtreme size
 1:                      5990980_Cell_Cycle 0.0001221001 0.002356366 0.5388497 2.682641            0  369
 2:             5990979_Cell_Cycle,_Mitotic 0.0001252976 0.002356366 0.5594755 2.738270            0  317
 3:        5991210_Signaling_by_Rho_GTPases 0.0001310444 0.002356366 0.4238512 2.009963            0  231
 4:                         5991454_M_Phase 0.0001372872 0.002356366 0.5576247 2.548794            0  173
 5:     5991023_Metabolism_of_carbohydrates 0.0001378930 0.002356366 0.4944766 2.240678            0  160
 6:            5991209_RHO_GTPase_Effectors 0.0001386001 0.002356366 0.5248796 2.369573            0  157
 7:  5991502_Mitotic_Metaphase_and_Anaphase 0.0001450326 0.002356366 0.6052907 2.633267            0  123
 8:                5991600_Mitotic_Anaphase 0.0001452222 0.002356366 0.6015521 2.613626            0  122
 9: 5991599_Separation_of_Sister_Chromatids 0.0001460707 0.002356366 0.6164600 2.661433            0  116
10:           5991130_Programmed_Cell_Death 0.0001465201 0.002356366 0.4537506 1.941445            0  108
```

One can make an enrichment plot for a pathway:
```{r}
plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death")

```

![enrichment.png](https://www.dropbox.com/s/zusn9pju7f608sn/enrichment.png?raw=1)

Or make a table plot for a bunch of selected pathways:
```{r}
topPathwaysUp <- fgseaRes[ES > 0, ][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0, ][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)
```

<img src="https://www.dropbox.com/s/uthtzn8wgo176f6/enrichmentTable.png?raw=1" width="600">
