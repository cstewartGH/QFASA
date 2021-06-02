#' Produces a dendrogram using distances between the mean
#' FA signatures of the prey types.
#'
#'
#' Performs a hierarchical cluster analysis of mean prey
#' fatty acid signatures using function \code{hclust}.
#' @export
#' @param prey.fa data frame of prey fatty acid signature
#'     samples. First column must be species used to group samples. Other columns
#'     are assumed to be fatty acid proportions.
#' @param method the agglomeration method to be used.  This should be
#'     one of the possible methods in \code{hclust} such as \code{"single"},
#'     \code{"complete"} or  \code{"average"}.  Default is \code{"complete"}.
#' @param dist.meas distance measure to use for calculating dissimilarities:
#'    1=KL, 2=AIT or 3=CS.  Default is \code{AIT}.
#' @return Plot (dendrogram)
#'
#' @examples
#'
#' ## Fatty Acids
#' data(FAset)
#' fa.set = as.vector(unlist(FAset))
#'
#' ## prey.cluster requires full prey database.
#' data(preyFAs)
#' prey.sub=(preyFAs[,4:(ncol(preyFAs))])[fa.set]
#' prey.sub=prey.sub/apply(prey.sub,1,sum)
#' group=as.vector(preyFAs$Species)
#' prey.matrix=cbind(group,prey.sub)
#'
#' prey.cluster(prey.matrix,method="average",dist.meas=3)

prey.cluster <- function(prey.fa,
                         method="complete",
                         dist.meas=2) {

  mult.rep <- function(x.vec) {

    no.zero <- sum(x.vec == 0)
    x.vec[x.vec == 0] <- 1e-05
    x.vec[x.vec > 0] <- (1 - no.zero * 1e-05) * x.vec[x.vec > 0]

    return(x.vec)

  }

  prey.fa[,-1] <- prey.fa[,-1]/apply(prey.fa[,-1],1,sum)
  sort.preytype <- order(prey.fa[,1])
  prey.fa <- prey.fa[sort.preytype,]


  prey.fa[,-1] <- t(apply(prey.fa[,-1],1,mult.rep))


#Calculate distances

prey.means <- MEANmeth(prey.fa)
I <- nrow(prey.means)


dist.mat <- matrix(rep(NA,I*I),nrow=I)

if (dist.meas==1) {

  my.ylab <- "KL Distance"

  for (i in 2:I) {

    for (j in 1:(i-1)) {


          dist.mat[i,j] <- KL.dist(as.vector(prey.means[i,]),as.vector(prey.means[j,]))

    }

  }

} else if (dist.meas==2) {

  my.ylab <- "Aitchison Distance"

  for (i in 2:I) {

   for (j in 1:(i-1)) {


     dist.mat[i,j] <- AIT.dist(as.vector(prey.means[i,]),as.vector(prey.means[j,]))

    }

  }

} else {

  my.ylab <- "CS Distance"

  for (i in 2:I) {

    for (j in 1:(i-1)) {


      dist.mat[i,j] <- chisq.dist(as.vector(prey.means[i,]),as.vector(prey.means[j,]),1)

    }

  }

}


plot(stats::hclust(stats::as.dist(dist.mat),method=method),labels=rownames(prey.means),xlab="", ylab=my.ylab, sub=NA)


}









