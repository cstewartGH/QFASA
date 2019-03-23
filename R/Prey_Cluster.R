#' This function performs a hierarchical cluster analysis of prey
#' fatty acid signatures using a matrix of dissimilarities for the n
#' objects being clustered.  Initially, each object is assigned as its
#' own cluster and then the algorithm proceeds iteratively, at each
#' stage joining the two most similar clusters, until there
#' is just a single cluster. 
#'
#' @export
#' @param prey.fa data frame of prey fatty acid signature
#'     samples. Species column is used to group samples. Other columns
#'     are assumed to be fatty acid proportions.
#' @param method the agglomeration method to be used.  This should be
#'     one of \code{'single'}, \code{'complete'}, \code{'average'}, \code{'median'}, \code{'centroid'}.
#' @param FUN distance function
#' @return an object of class \code{hclust} which describes the tree produced by
#'           the clustering process.
#'
prey.cluster <- function(prey.fa,
                         method,
                         FUN) {}
