#' @useDynLib fgsea
NULL

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
    rCumSum <- cumsum(rAdj) / NR

    tops <- rCumSum - (S - seq_along(S)) / (N - m)
    bottoms <- tops - rAdj / NR
    maxP <- max(tops)
    minP <- min(bottoms)

    if(maxP > -minP) {
        geneSetStatistic <- maxP
    } else {
        geneSetStatistic <- minP
    }
    geneSetStatistic
}

#' @export
#' @import data.table
#' @import parallel
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
    npermActual <- nperm * m
    message(sprintf("%s pathways left", m))
    message(sprintf("Setting actual permutations number to %s", npermActual))

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
