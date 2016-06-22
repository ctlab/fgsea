#' @useDynLib fgsea
#' @import Rcpp
NULL

#' Example vector of gene-level statistics obtained for Th1 polarization.
#'
#' The data were obtained by doing differential expression between
#' Naive and Th1-activated states for GEO dataset GSE14308. The
#' exact script is available as system.file("gen_gene_ranks.R", package="fgsea")
#' @docType data
#' @name exampleRanks
NULL

#' Example list of mouse Reactome pathways.
#'
#' The list was obtained by selecting all the pathways from `reactome.db`
#' package that contain mouse genes. The exact script is available as
#' system.file("gen_reactome_pathways.R", package="fgsea")
#' @docType data
#' @name examplePathways
NULL
