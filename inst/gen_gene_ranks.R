library(GEOquery)
library(limma)
library(org.Mm.eg.db)
library(data.table)
# for collapseBy
source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")

gse14308 <- getGEO("GSE14308")[[1]]

pData(gse14308)$condition <- sub("-.*$", "", gse14308$title)

es <- collapseBy(gse14308, fData(gse14308)$ENTREZ_GENE_ID, FUN=median)
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]

exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")

es.design <- model.matrix(~0+condition, data=pData(es))

fit <- lmFit(es, es.design)

fit2 <- contrasts.fit(fit, makeContrasts(conditionTh1-conditionNaive,
                                         levels=es.design))
fit2 <- eBayes(fit2)
de <- data.table(topTable(fit2, adjust.method="BH", number=12000, sort.by = "B"), keep.rownames = TRUE)


ranks <- de[order(t), list(rn, t)]
write.table(ranks, "inst/extdata/naive.vs.th1.rnk", sep="\t", quote = FALSE, row.names = FALSE)
