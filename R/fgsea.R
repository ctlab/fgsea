#' Calculates GSEA statistics for a given query gene set
#'
#' Takes \emph{O(k log k)} time, where \emph{k} is a size of `selectedSize`.
#' @param stats Named numeric vector with gene-level statistics
#'  sorted in decreasing order (order is not checked).
#' @param selectedStats Indexes of selected genes in the `stats` array.
#' @param gseaParam GSEA weight parameter (0 is unweighted, suggested value is 1).
#' @param returnAllExtremes If TRUE return not only the most extreme point, but all of them. Can be used for enrichment plot
#' @param returnLeadingEdge If TRUE return also leading edge genes.
#' @return Value of GSEA statistic if both returnAllExtremes and returnLeadingEdge are FALSE.
#' Otherwise returns list with the folowing elements:
#' \itemize{
#' \item res -- value of GSEA statistic
#' \item tops -- vector of top peak values of cumulative enrichment statistic for each gene;
#' \item bottoms -- vector of bottom peak values of cumulative enrichment statistic for each gene;
#' \item leadingGene -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#' @export
#' @examples
#' data(exampleRanks)
#' data(examplePathways)
#' ranks <- sort(exampleRanks, decreasing=TRUE)
#' es <- calcGseaStat(ranks, na.omit(match(examplePathways[[1]], names(ranks))))
calcGseaStat <- function(stats, selectedStats, gseaParam=1,
                         returnAllExtremes=FALSE,
                         returnLeadingEdge=FALSE) {
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

    if (!returnAllExtremes && !returnLeadingEdge) {
        return(geneSetStatistic)
    }

    res <- list(res=geneSetStatistic)
    if (returnAllExtremes) {
        res <- c(res, list(tops=tops, bottoms=bottoms))
    }
    if (returnLeadingEdge) {
        leadingEdge <- if (maxP > -minP) {
            S[seq_along(S) <= which.max(bottoms)]
        } else if (maxP < -minP) {
            S[seq_along(S) >= which.min(bottoms)]
        } else {
            NULL
        }

        res <- c(res, list(leadingEdge=leadingEdge))
    }
    res
}

#' Runs preranked gene set enrichment analysis.
#'
#' The function takes about \emph{O(nk^\{3/2\})} time,
#' where \emph{n} is number of permutations and \emph{k} is a maximal
#' size of the pathways. That means that setting `maxSize` parameter with a value of ~500
#' is strongly recommended.
#' @param pathways List of gene sets to check.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param nproc Number of parallel processes to use.
#' @param gseaParam GSEA parameter value.
#' @return A table with GSEA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'  \item pathway -- name of the pathway as in `names(pathway)`;
#'  \item pval -- an enrichment p-value;
#'  \item padj -- a BH-adjusted p-value;
#'  \item ES -- enrichment score, same as in Broad GSEA implementation;
#'  \item NES -- enrichment score normalized to mean enrichment of random samples of the same size;
#'  \item nMoreExtreme` -- a number of times a random gene set had a more
#'      extreme enrichment score value;
#'  \item size -- size of the pathway after removing genes not present in `names(stats)`.
#' }
#'
#' @export
#' @import data.table
#' @import parallel
#' @import stats
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=10000, maxSize=500)
#' # Testing only one pathway is implemented in a more efficient manner
#' fgseaRes1 <- fgsea(examplePathways[1], exampleRanks, nperm=10000)
fgsea <- function(pathways, stats, nperm,
                  minSize=1, maxSize=Inf,
                  nproc=1,
                  gseaParam=1) {
    minSize <- max(minSize, 1)
    stats <- sort(stats, decreasing=TRUE)
    pathwaysFiltered <- lapply(pathways, function(p) { as.vector(na.omit(match(p, names(stats)))) })
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

    gseaStatRes <- do.call(rbind,
                lapply(pathwaysFiltered, calcGseaStat,
                       stats=stats,
                       returnLeadingEdge=TRUE))


    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(gseaStatRes[, "res"])



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

    pvals[, leadingEdge := .(leadingEdges)]


    # Makes pvals object printable immediatly
    pvals <- pvals[]

    pvals
}
