#' This function performs a hierarchical cluster analysis of prey
#' fatty acid signatures using a set  of dissimilarities for the n
#' objects being clustered.  Initially, each object is assigned to its
#' own cluster and then the algorithm proceeds iteratively, at each
#' stage joining the two most similar clusters, continuing until there
#' is just a single cluster. 
#'
#' @export
#' @param prey.fa data frame of prey fatty acid signature
#'     samples. Species column is used to group samples. Other columns
#'     are assumed to be fatty acid proportions.
#' @param method the agglomeration method to be used.  This should be
#'     one of "single", "complete", "average", "median", "centroid".
#' @param FUN distance function
#' @return an object of class *hclust* which describes the tree produced by
#'           the clustering process.
#'
prey.cluster <- function(prey.fa,
                         method,
                         FUN) {}
