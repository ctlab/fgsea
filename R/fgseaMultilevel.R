#' Runs preranked gene set enrichment analysis.
#'
#' This feature is based on the adaptive multilevel splitting Monte Carlo approach.
#' This allows us to exceed the results of simple sampling and calculate arbitrarily small P-values.
#' @param pathways List of gene sets to check.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param sampleSize The size of a random set of genes which in turn has size = pathwaySize
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param absEps This parameter sets the boundary for calculating the p value.
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#'@export
#' @import fastmatch
#' @import data.table
#' @import BiocParallel
#' @return A table with GSEA results. Each row corresponds to a tested pathway. The columns are the following
#' \itemize{
#' \item pathway -- name of the pathway as in `names(pathway)`;
#' \item pval -- an enrichment p-value;
#' \item padj -- a BH-adjusted p-value;
#' \item log2err -- the expected error for the standard deviation of the P-value logarithm.
#' \item ES -- enrichment score, same as in Broad GSEA implementation;
#' \item NES -- enrichment score normalized to mean enrichment of random samples of the same size;
#' \item size -- size of the pathway after removing genes not present in `names(stats)`.
#' \item leadingEdge -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' fgseaMultilevelRes <- fgseaMultilevel(examplePathways, exampleRanks, maxSize=500)
fgseaMultilevel <- function(pathways, stats, sampleSize=101,
                            minSize=1, maxSize=Inf, absEps=0,
                            nproc=0, BPPARAM=NULL)
{
    # Error if pathways is not a list
    if (!is.list(pathways)) {
        stop("pathways should be a list with each element containing names of the stats argument")
    }

    # Error if stats is not named
    if (is.null(names(stats))) {
        stop("stats should be named")
    }

    # Warning message for ties in stats
    ties <- sum(duplicated(stats[stats != 0]))
    if (ties != 0) {
        warning("There are ties in the preranked stats (",
                paste(round(ties * 100 / length(stats), digits = 2)),
                "% of the list).\n",
                "The order of those tied genes will be arbitrary, which may produce unexpected results.")
    }
    # Warning message for duplicate gene names
    if (any(duplicated(names(stats)))) {
        warning("There are duplicate gene names, fgsea may produce unexpected results")
    }
    #To avoid warnings during the check
    log2err=nMoreExtreme=pathway=pval=padj=NULL
    ES=NES=size=leadingEdge=NULL
    .="damn notes"

    nPermSimple <- 1000 # number of samples for initial fgseaSimple run: fast and good enough
    minSize <- max(minSize, 1)
    stats <- sort(stats, decreasing = TRUE)
    if (sampleSize %% 2 == 0){
        sampleSize <-  sampleSize + 1
    }
    pathwaysFiltered <- lapply(pathways, function(p) { as.vector(na.omit(fmatch(p, names(stats)))) })
    pathwaysSizes <- sapply(pathwaysFiltered, length)

    toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
    m <- length(toKeep)
    if (m == 0) {
        return(data.table(pathway=character(),
                          pval=numeric(),
                          padj=numeric(),
                          log2err=numeric(),
                          ES=numeric(),
                          NES=numeric(),
                          size=integer(),
                          leadingEdge=list()))
    }
    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]

    gseaStatRes <- do.call(rbind,
                           lapply(pathwaysFiltered, calcGseaStat,
                                  stats=stats,
                                  returnLeadingEdge=TRUE))

    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(gseaStatRes[, "res"])

    universe <- seq_along(stats)

    seeds <- sample.int(10^9, 1)
    BPPARAM <- setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)

    simpleFgseaRes <- fgseaSimpleImpl(pathwayScores=pathwayScores, pathwaysSizes=pathwaysSizes,
                                      pathwaysFiltered=pathwaysFiltered, leadingEdges=leadingEdges,
                                      permPerProc=nPermSimple, seeds=seeds, toKeepLength=m,
                                      stats=stats, BPPARAM=SerialParam())
    simpleError <- 1/log(2)*(trigamma(simpleFgseaRes$nMoreExtreme) - trigamma(nPermSimple))
    multError <- sapply((simpleFgseaRes$nMoreExtreme + 1) / nPermSimple, multilevelError, sampleSize)


    if (all(multError > simpleError)){
        simpleFgseaRes[, log2err := 1/log(2)*(trigamma(nMoreExtreme) - trigamma((nPermSimple)))]
        setorder(simpleFgseaRes, pathway)
        simpleFgseaRes <- simpleFgseaRes[, .(pathway, pval, padj, log2err, ES, NES, size, leadingEdge)]
        simpleFgseaRes <- simpleFgseaRes[]
        return(simpleFgseaRes)
    }

    dtSimpleFgsea <- simpleFgseaRes[simpleError < multError]
    dtSimpleFgsea[, log2err := 1/log(2)*(trigamma(nMoreExtreme) - trigamma(nPermSimple))]
    dtMultilevel <- simpleFgseaRes[multError < simpleError]

    multilevelPathwaysList <- split(dtMultilevel, by="size")
    # In most cases, this gives a speed increase with parallel launches.
    indxs <- sample(1:length(multilevelPathwaysList))
    multilevelPathwaysList <- multilevelPathwaysList[indxs]

    seed=sample.int(1e9, size=1)
    pvals <- multilevelImpl(multilevelPathwaysList, stats, sampleSize,
                            seed, absEps, BPPARAM=BPPARAM)
    result <- rbindlist(multilevelPathwaysList)
    result[, pval := unlist(pvals)]
    result[, padj := p.adjust(pval, method = "BH")]
    result[, log2err := sqrt(floor(-log2(pval) + 1) * (trigamma((sampleSize+1)/2) - trigamma(sampleSize+1))/log(2))]
    result <- rbindlist(list(result, dtSimpleFgsea), use.names = TRUE)
    result <- result[, .(pathway, pval, padj, log2err, ES, NES, size, leadingEdge)]
    setorder(result, pathway)
    result <- result[]
    result
}

#'Calculates the expected error for the standard deviation of the P-value logarithm.
#'
#' @param pval P-value
#' @param sampleSize equivavlent to sampleSize in fgseaMultilevel
#'@export
#'@return The value of the expected error
#' @examples
#' expectedError <- multilevelError(pval=1e-10, sampleSize=1001)
multilevelError <- function(pval, sampleSize){
    return(sqrt(floor(-log2(pval) + 1) * (trigamma((sampleSize+1)/2) - trigamma(sampleSize+1))/log(2)))
}

#' Calculates P-values for preprocessed data.
#' @param multilevelPathwaysList List of pathways for which P-values will be calculated.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param sampleSize The size of a random set of genes which in turn has size = pathwaySize
#' @param seed `seed` parameter from `fgseaMultilevel`
#' @param absEps This parameter sets the boundary for calculating the p value.
#' @param sign This option will be used in future implementations.
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @return List of P-values.
multilevelImpl <- function(multilevelPathwaysList, stats, sampleSize,
                           seed, absEps, sign=FALSE, BPPARAM=NULL){
    #To avoid warnings during the check
    size=ES=NULL
    pvals <- bplapply(multilevelPathwaysList,
                      function(x) fgseaMultilevelCpp(x[, ES], stats, unique(x[, size]),
                                                     sampleSize, seed, absEps, sign),
                      BPPARAM=BPPARAM)
    return(pvals)
}
