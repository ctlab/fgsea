fetchSubDataTables <- function(pDt){
    subDtPos <- pDt[pDt$ES >= 0, ]
    if (nrow(subDtPos) > 0){
        subDtPos <- setorder(subDtPos, -size, -ES)
        subDtPos <- split(subDtPos, subDtPos$size)
        subDtPos <- subDtPos[order(sample(1:length(subDtPos)))]
    } else{
        subDtPos <- NULL
    }
    subDtNeg <- pDt[pDt$ES < 0, ]
    if (nrow(subDtNeg) > 0){
        subDtNeg <- setorder(subDtNeg, -size, ES)
        subDtNeg <- split(subDtNeg, subDtNeg$size)
        subDtNeg <- subDtNeg[order(sample(1:length(subDtNeg)))]
    } else{
        subDtNeg <- NULL
    }
    return(list(pos=subDtPos, neg=subDtNeg))
}

lcalcPvals <- function(subDt, stats, samplesSize, seed, absEps, sign=FALSE, BPPARAM=NULL){
    if (is.null(subDt[["pos"]])){
        pvals <- bplapply(
            subDt[["neg"]],
            function(x) fgseaMultilevelCpp(as.numeric(x[, size]), as.numeric(x[, ES]),
                                         stats, samplesSize, seed, absEps, sign),
            BPPARAM=BPPARAM
        )
        return(list(pos=NULL, neg=pvals))
    } else if (is.null(subDt[["neg"]])){
        pvals <- bplapply(
            subDt[["pos"]],
            function(x) fgseaMultilevelCpp(as.numeric(x[, size]), as.numeric(x[, ES]),
                                         stats, samplesSize, seed, absEps, sign),
            BPPARAM=BPPARAM
        )
        return(list(pos=pvals, neg=NULL))
    }
    else{
        resPos <- bplapply(
            subDt[["pos"]],
            function(x) fgseaMultilevelCpp(as.numeric(x[, size]), as.numeric(x[, ES]),
                                         stats, samplesSize, seed, absEps, sign),
            BPPARAM=BPPARAM
        )
        resNeg <- bplapply(
            subDt[["neg"]],
            function(x) fgseaMultilevelCpp(as.numeric(x[, size]), as.numeric(x[, ES]),
                                         stats, samplesSize, seed, absEps, sign),
            BPPARAM=BPPARAM
        )
        return(list(posEsPval=resPos, negEsPval=resNeg))
    }

}

#'Runs preranked gene set enrichment analysis.
#'
#' This approach is based on multilevel splitting approach,
#' that server for estimating the probability of rare events.
#' @param pathways List of gene sets to check.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param samplesSize The size of a random set of genes which in turn has size = pathwaySize
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param absEps This parameter sets the boundary for calculating the p value.
#' @param signed If TRUE returns two probabilities with a sign dependency.
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param BPPARAM Parallelization parameter used in bplapply.
#'@export
#' @import fastmatch
#' @import data.table
#' @import BiocParallel
#'@return A table with GSEA results. Each row corresponds to a tested pathway.
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' fgseaMultilevelRes <- fgseaMultilevel(examplePathways, exampleRanks, maxSize=500)
fgseaMultilevel <- function(pathways, stats, samplesSize=101,
                            minSize=1, maxSize=Inf, absEps=0,
                            signed = FALSE, nproc=0, BPPARAM=NULL)
{
    minSize <- max(minSize, 1)
    stats <- sort(stats, decreasing = TRUE)
    if (samplesSize %% 2 == 0){
        samplesSize <-  samplesSize + 1
    }
    pathwaysFiltered <- lapply(pathways, function(p) { as.vector(na.omit(fmatch(p, names(stats)))) })
    pathwaysSizes <- sapply(pathwaysFiltered, length)

    toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
    m <- length(toKeep)

    if (m == 0) {
        return(data.table(pathway=character(),
                          pval=numeric(),
                          padj=numeric(),
                          ES=numeric(),
                          size=integer(),
                          leadingEdge=list()))
    }
    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]

    K <- max(pathwaysSizes)

    gseaStatRes <- do.call(rbind,
                           lapply(pathwaysFiltered, calcGseaStat,
                                  stats=stats,
                                  returnLeadingEdge=TRUE))

    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(gseaStatRes[, "res"])


    nperm <- 1000
    granularity <- 1000
    permPerProc <- rep(granularity, floor(nperm / granularity))
    if (nperm - sum(permPerProc) > 0) {
        permPerProc <- c(permPerProc, nperm - sum(permPerProc))
    }

    seeds <- sample.int(10^9, 1)

    if (is.null(BPPARAM)) {
        if (nproc != 0) {
            if (.Platform$OS.type == "windows") {
                # windows doesn't support multicore, using snow instead
                BPPARAM <- SnowParam(workers = nproc)
            } else {
                BPPARAM <- MulticoreParam(workers = nproc)
            }
        } else {
            BPPARAM <- bpparam()
        }
    }

    universe <- seq_along(stats)

    counts <- bplapply(seq_along(permPerProc), function(i) {
        nperm1 <- permPerProc[i]
        leEs <- rep(0, m)
        geEs <- rep(0, m)
        leZero <- rep(0, m)
        geZero <- rep(0, m)
        leZeroSum <- rep(0, m)
        geZeroSum <- rep(0, m)
        if (m == 1) {
            for (i in seq_len(nperm1)) {
                randSample <- sample.int(length(universe), K)
                randEsP <- calcGseaStat(
                    stats = stats,
                    selectedStats = randSample,
                    gseaParam = 1)
                leZero <- leZero + (randEsP <= 0)
                geZero <- geZero + (randEsP >= 0)
                leZeroSum <- leZeroSum + pmin(randEsP, 0)
                geZeroSum <- geZeroSum + pmax(randEsP, 0)
            }
        } else {
            aux <- calcGseaStatCumulativeBatch(
                stats = stats,
                gseaParam = 1,
                pathwayScores = pathwayScores,
                pathwaysSizes = pathwaysSizes,
                iterations = nperm1,
                seed = seeds[i])
            leEs = get("leEs", aux)
            geEs = get("geEs", aux)
            leZero = get("leZero", aux)
            geZero = get("geZero", aux)
            leZeroSum = get("leZeroSum", aux)
            geZeroSum = get("geZeroSum", aux)
        }
        data.table(pathway=seq_len(m),
                   leEs=leEs, geEs=geEs,
                   leZero=leZero, geZero=geZero,
                   leZeroSum=leZeroSum, geZeroSum=geZeroSum
        )
    }, BPPARAM = BPPARAM)

    counts <- rbindlist(counts)
    dtES <- counts[, list(leZeroMean = sum(leZeroSum) / sum(leZero),
                          geZeroMean = sum(geZeroSum) / sum(geZero)), by=.(pathway)]
    dtES[, ES := pathwayScores[pathway]]
    dtES[, NES := ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean))]


    pathwaysDt <-  data.table(pathway = names(pathwaysFiltered), size = pathwaysSizes,
                         ES = pathwayScores, NES = dtES$NES,
                         leadingEdge = I(leadingEdges))
    subDt <- fetchSubDataTables(pathwaysDt)

    if (is.null(subDt[["pos"]])){
        result <- rbindlist(subDt[["neg"]])
    } else if (is.null(subDt[["neg"]])){
        result <- rbindlist(subDt[["pos"]])
    } else{
        result <- rbindlist(c(subDt[["pos"]], subDt[["neg"]]))
    }

    seed=sample.int(1e9, size=1)
    if (signed){
        pvals_plus <- lcalcPvals(subDt, stats, samplesSize,
                            seed, absEps, sign=TRUE, BPPARAM=BPPARAM)
        stats <- sort(stats)
        pvals_minus <- lcalcPvals(subDt, stats, samplesSize,
                                 seed, absEps, sign=TRUE, BPPARAM=BPPARAM)
    } else{
        pvals <- lcalcPvals(subDt, stats, samplesSize,
                            seed, absEps, sign=FALSE, BPPARAM=BPPARAM)
    }
    if (signed){
        result[, p_plus := unlist(pvals_plus)]
        result[, p_minus := unlist(pvals_minus)]
        result[, p_minus_adj := p.adjust(p_minus, method = "BH")]
        result[, p_plus_adj := p.adjust(p_plus, method = "BH")]
        setcolorder(result, c("pathway", "p_plus", "p_plus_adj", "p_minus", "p_minus_adj", "ES", "NES", "size", "leadingEdge"))
        setorder(result, pathway)
    } else{
        result[, pval := unlist(pvals)]
        result[, padj := p.adjust(pval, method = "BH")]
        result[, log2err := sqrt(floor(-log2(pval) + 1) * (trigamma((samplesSize+1)/2) - trigamma(samplesSize+1))/log(2))]
        setcolorder(result, c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge"))
        setorder(result, pathway)
    }
    result <- result[]

    result
}

#'Calculates the expected error.
#'
#' @param pval P-value
#' @param C equivavlent to samplesSize
#'@export
#'@return The value of the expected error
#' @examples
#' expectedError <- multilevelError(pval=1e-10, C=1001)
multilevelError <- function(pval, C){
    return(sqrt(floor(-log2(pval) + 1) * (trigamma((C+1)/2) - trigamma(C+1))/log(2)))
}
