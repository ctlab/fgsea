checkGesecaArgs <- function(E, pathways){
    # Error if pathways is not a list
    if (!is.list(pathways)) {
        stop("pathways should be a list with each element containing names of the stats argument")
    }

    # Error if stats is not named
    if (is.null(rownames(E))) {
        stop("E rows should be named")
    }

    # Error if E has non-finite values
    if (any(!is.finite(E))){
        stop("Not all E values are finite numbers")
    }

    # Warning message for duplicate gene names
    if (any(duplicated(rownames(E)))) {
        warning("There are duplicate gene names, geseca may produce unexpected results.")
    }
}

gesecaPreparePathways <- function(E, pathways, minSize, maxSize){
    minSize <- max(minSize, 2)
    pathwaysFiltered <- lapply(pathways, function(p) {unique(na.omit(fmatch(p, rownames(E))))})
    pathwaysSizes <- sapply(pathwaysFiltered, length)

    toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]

    return(list(filtered=pathwaysFiltered,
                sizes=pathwaysSizes))
}


calcGesecaScores <- function(indxs, E){
    return(var(colSums(E[indxs, ])))
}
