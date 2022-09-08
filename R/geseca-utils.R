checkGesecaArgs <- function(E, pathways){
    # Error if pathways is not a list
    if (!is.list(pathways)) {
        stop("pathways should be a list with each element containing names of the stats argument")
    }

    # Error if E matrix doesn't have rownames
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
    maxSize <- min(nrow(E) - 1, maxSize)

    pathwaysFiltered <- lapply(pathways, function(p) {unique(na.omit(fmatch(p, rownames(E))))})
    pathwaysSizes <- sapply(pathwaysFiltered, length)

    toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]

    return(list(filtered=pathwaysFiltered,
                sizes=pathwaysSizes))
}


# E should be row-centered
calcGesecaScores <- function(indxs, E){
    return(sum(colSums(E[indxs, ])**2))
}

# E should be row-centered
getGesecaLeadingEdge <- function(indxs, E, extend=FALSE, sort=TRUE) {
    profile <- colSums(E[indxs, ])
    profile <- profile / sqrt(sum(profile**2))
    Ey <- E[indxs,]
    if (extend) {
        Ey <- E
    }
    weights <- (Ey %*% profile)[,1]
    weights <- weights / sqrt(apply(Ey**2, 1, sum))
    if (sort) {
        weights <- base::sort(weights, decreasing = TRUE)
    }
    return(names(weights[weights >= 0.5]))
}
