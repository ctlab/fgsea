#' Calculates GSEA statistics for a given query gene set
#' @param stats Named numeric vector with gene-level statistics
#' @param selectedStats indexes of selected genes in a 'stats' array
#' @param gseaParam GSEA weight parameter (0 is unweighted, suggested value is 1)
#' @return value of GSEA statistic
#' @export
calcGseaStat <- function(stats, selectedStats, gseaParam=1) {
    S <- selectedStats
    r <- stats
    p <- gseaParam

    S <- sort(S)

    m <- length(S)
    N <- length(r)
    NR <- (sum(abs(r[S])^p))
    rAdj <- abs(r[S])^p
    if (NR == 0) {
        # this is equivalent to rAdj being rep(eps, m)
        rCumSum <- seq_along(rAdj) / length(rAdj)
    } else {
        rCumSum <- cumsum(rAdj) / NR
    }


    tops <- rCumSum - (S - seq_along(S)) / (N - m)
    if (NR == 0) {
        # this is equivalent to rAdj being rep(eps, m)
        bottoms <- tops - 1 / m
    } else {
        bottoms <- tops - rAdj / NR
    }

    maxP <- max(tops)
    minP <- min(bottoms)

    if(maxP > -minP) {
        geneSetStatistic <- maxP
    } else if (maxP < -minP) {
        geneSetStatistic <- minP
    } else {
        geneSetStatistic <- 0
    }
    geneSetStatistic
}

#' Runs preranked gene set enrichment analysis.
#' @param pathways List of gene sets to check.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param nproc Number of parallel processes to use.
#' @param gseaParam GSEA parameter value.
#' @return A table with GSEA results.
#' @export
#' @import data.table
#' @import parallel
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' \dontrun{
#'  fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=10000, maxSize=500)
#' }
fgsea <- function(pathways, stats, nperm,
                  minSize=1, maxSize=Inf,
                  nproc=1,
                  gseaParam=1) {
    minSize <- max(minSize, 1)
    stats <- sort(stats, decreasing=T)
    pathwaysFiltered <- lapply(pathways, intersect, names(stats))
    pathwaysSizes <- sapply(pathwaysFiltered, length)
    pathwaysFiltered <- pathwaysFiltered[
        minSize <= pathwaysSizes & pathwaysSizes <= maxSize]

    pathwaysSizes <- sapply(pathwaysFiltered, length)
    K <- max(pathwaysSizes)
    m <- length(pathwaysFiltered)
    npermActual <- nperm
#     npermActual <- if (npermIsActual) nperm else nperm * m
#     message(sprintf("%s pathways left", m))
#     message(sprintf("Setting actual permutations number to %s", npermActual))

    pathwaysIndexes <- lapply(pathwaysFiltered, match, names(stats))

    pathwayScores <- sapply(pathwaysIndexes, calcGseaStat, stats=stats)


    permPerProc <- rep(npermActual %/% nproc, nproc) +
        c(rep(1, npermActual %% nproc), rep(0, nproc - npermActual %% nproc))
    counts <- mclapply(permPerProc, function(nperm1) {
        leEs <- rep(0, m)
        geEs <- rep(0, m)
        leZero <- rep(0, m)
        geZero <- rep(0, m)
        leZeroSum <- rep(0, m)
        geZeroSum <- rep(0, m)
        for (i in seq_len(nperm1)) {
            randSample <- sample(seq_along(stats), size=K)
            if (m == 1) {
                randEsP <- calcGseaStat(
                    stats = stats,
                    selectedStats = randSample,
                    gseaParam = gseaParam)
            } else {
                randEs <- calcGseaStatCumulative(
                    stats = stats,
                    selectedStats = randSample,
                    gseaParam = gseaParam)
                randEsP <- randEs[pathwaysSizes]
            }
            leEs <- leEs + (randEsP <= pathwayScores)
            geEs <- geEs + (randEsP >= pathwayScores)
            leZero <- leZero + (randEsP <= 0)
            geZero <- geZero + (randEsP >= 0)
            leZeroSum <- leZeroSum + pmin(randEsP, 0)
            geZeroSum <- geZeroSum + pmax(randEsP, 0)
        }
        data.table(pathway=seq_len(m),
                   leEs=leEs, geEs=geEs,
                   leZero=leZero, geZero=geZero,
                   leZeroSum=leZeroSum, geZeroSum=geZeroSum
                   )
    }, mc.cores=nproc)

    counts <- rbindlist(counts)

    # Getting rid of check NOTEs
    leEs=leZero=geEs=geZero=leZeroSum=geZeroSum=NULL
    pathway=padj=pval=ES=NES=geZeroMean=leZeroMean=NULL
    nMoreExtreme=nGeEs=nLeEs=size=NULL
    .="damn notes"


    pvals <- counts[,
        list(pval=min((1+sum(leEs)) / (1 + sum(leZero)),
                 (1+sum(geEs)) / (1 + sum(geZero))),
             leZeroMean = sum(leZeroSum) / sum(leZero),
             geZeroMean = sum(geZeroSum) / sum(geZero),
             nLeEs=sum(leEs),
             nGeEs=sum(geEs)
             )
        ,
        by=.(pathway)]
    pvals[, padj := p.adjust(pval, method="BH")]
    pvals[, ES := pathwayScores[pathway]]
    pvals[, NES := ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean))]
    pvals[, leZeroMean := NULL]
    pvals[, geZeroMean := NULL]

    pvals[, nMoreExtreme :=  ifelse(ES > 0, nGeEs, nLeEs)]
    pvals[, nLeEs := NULL]
    pvals[, nGeEs := NULL]

    pvals[, size := pathwaysSizes[pathway]]
    pvals[, pathway := names(pathwaysFiltered)[pathway]]

    pvals
}


#' Runs preranked gene set enrichment analysis for reactome pathways.
#' @param stats Named vector of gene-level stats. Names should be Entrez IDs.
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param nproc Number of parallel processes to use.
#' @param gseaParam GSEA parameter value.
#' @return A table with GSEA results.
#' @export
#' @examples
#' data(exampleRanks)
#' \dontrun{
#'  fgseaRes <- fgseaReactome(exampleRanks, nperm=10000, maxSize=500)
#' }
fgseaReactome <- function(stats, nperm,
                  minSize=1, maxSize=Inf,
                  nproc=1,
                  gseaParam=1) {
    stopifnot(requireNamespace("reactome.db"))
    stopifnot(requireNamespace("AnnotationDbi"))
    pathways <- na.omit(AnnotationDbi::select(reactome.db::reactome.db,
                                              keys=names(stats),
                                              c("PATHID"),
                                              keytype = 'ENTREZID'))

    pathways <- split(pathways$ENTREZID, pathways$PATHID)

    pathway2name <- as.data.table(na.omit(AnnotationDbi::select(
        reactome.db::reactome.db,
        names(pathways),
        c("PATHNAME"),
        'PATHID')))

    # Hack to get rid of check NOTEs
    pathwayId=pathway=PATHNAME=NULL

    # Remove organism prefix
    pathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]
    pathway2name <- structure(pathway2name$PATHNAME, names=pathway2name$PATHID)
    date()

    res <- fgsea(pathways, stats, nperm,
                 minSize, maxSize,
                 nproc, gseaParam)
    res[, pathwayId := pathway]
    res[, pathway := pathway2name[pathwayId]]
    res
}
