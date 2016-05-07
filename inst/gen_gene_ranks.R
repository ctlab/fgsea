library(GEOquery)
library(limma)
library(org.Mm.eg.db)
library(data.table)
# for collapseBy and normalizeGeneDE
source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")

gse14308 <- getGEO("GSE14308")[[1]]

pData(gse14308)$condition <- sub("-.*$", "", gse14308$title)

es <- collapseBy(gse14308, fData(gse14308)$ENTREZ_GENE_ID, FUN=median)
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]

# our symbols are better
fData(es) <- data.frame(row.names = rownames(es))
fData(es)$symbol <- sapply(mget(rownames(es),
                                org.Mm.egSYMBOL,
                                ifnotfound = NA),
                           head, n=1)

exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")

es.design <- model.matrix(~0+condition, data=pData(es))

fit <- lmFit(es, es.design)

fit2 <- contrasts.fit(fit, makeContrasts(conditionTh1-conditionNaive,
                                         levels=es.design))
fit2 <- eBayes(fit2)
de <- normalizeGeneDE(topTable(fit2, adjust.method="BH", number=12000, sort.by = "B"))


ranks <- de[order(t), list(ID, t)]
write.table(ranks, "inst/naive.vs.th1.rnk", sep="\t", quote = F, row.names = F)
