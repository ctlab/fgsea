#' Returns a list of Reactome pathways for given Entrez gene IDs
#' @param genes Entrez IDs of query genes.
#' @return A list of vectors with gene sets.
#' @export
#' @examples
#' data(exampleRanks)
#' pathways <- reactomePathways(names(exampleRanks))
reactomePathways <- function(genes) {
    stopifnot(requireNamespace("reactome.db"))
    stopifnot(requireNamespace("AnnotationDbi"))
    pathways <- na.omit(AnnotationDbi::select(reactome.db::reactome.db,
                                              keys=genes,
                                              c("PATHID"),
                                              keytype = 'ENTREZID'))

    pathways <- split(pathways$ENTREZID, pathways$PATHID)

    pathway2name <- as.data.table(AnnotationDbi::select(
        reactome.db::reactome.db,
        names(pathways),
        c("PATHNAME"),
        'PATHID'))

    # workaround for reactome.db having one-to-many PATHID to PATHNAM mapping
    PATHID=NULL
    pathway2name <- pathway2name[!duplicated(PATHID)]



    # Hack to get rid of check NOTEs
    PATHNAME=NULL

    # Remove organism prefix
    pathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]

    # workaround for reactome.db having duplicated names for pathways
    # example: Immune system (5991156 and 6096039)
    name2pathways <- split(pathway2name$PATHID, pathway2name$PATHNAME)

    pathways <- lapply(name2pathways, function(x) unique(do.call(c, pathways[x])))

    pathways <- pathways[!is.na(names(pathways))]

    pathways
}

#' Returns a list of pathways from a GMT file.
#' @param gmt.file Path to a GMT file.
#' @return A list of vectors with gene sets.
#' @export
#' @examples
#' pathways <- gmtPathways(system.file(
#'      "extdata", "mouse.reactome.gmt", package="fgsea"))
#' @export
gmtPathways <- function(gmt.file) {
    pathwayLines <- strsplit(readLines(gmt.file), "\t")
    pathways <- lapply(pathwayLines, utils::tail, -2)
    names(pathways) <- sapply(pathwayLines, head, 1)
    pathways
}
