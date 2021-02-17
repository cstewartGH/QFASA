#' Generate a pseudo predator by sampling with replacement from prey database.
#'
#' Generates a single pseudo predator by sampling with replacement from prey database. To generate a sample of
#' pseudo predators, please refer to example code.
#'
#' @export
#' @param diet the "true" of "desired" diet of the pseudo predator. A compositional vector of proportions that sum to one with
#'     length equal to the number of prey species.
#' @param preybase prey database from which to generate the pseudo predator.  First column must provide the species name.
#' @param cal.vec vector of calibration coefficients whose length is the same as the number of fatty acids in prey database.
#' @param fat.vec vector of fat content whose length is the same as the number of species.
#' @param preysize number of prey to sample from prey database.  If preysize=1, then one prey is selected from each species.
#'       Otherwise, a sample of n_k signatures (where n_k is sample size for species k) is obtained by sampling with replacement.
#' @return A simulated predator FA signature.
#' @details The default is to re-sample all of the prey signatures within each species
#' (that is, preysize=2).  Alternatively, one prey may be randomly selected from each species yielding potentially more variable
#' pseudo-predators.
#' For details on simulating realistic predators signatures, see
#' Bromaghin, J. (2015)
#' Simulating realistic predator signatures in quantitative fatty acid signature analysis,
#' Ecological Informatics, 30, 68-71.
#'
#'
#' @examples
#' data(preyFAs)
#'
#' # Generating a sample of 10 pseudo predators each with "true" diet being
#' # (1/11,1/11,...1/11), no calibration effect and no fat content.  The QFASA diet estimate
#' # is then computed for each pseudo predator.
#'
#' # Note: To incorporate calibration and fat content in a simulation study,
#' # one set of calibration and fat content is generally used to simulate the pseudo predator
#' # and another is used to estimate the diet.
#'
#' set.seed(11)
#' p.mat <- matrix(rep(NA,10*11),nrow=10)
#' for (i in 1: 10) {
#'     my.seal <- pseudo.pred(rep(1/11,11),
#'                             preyFAs[,-c(1,3)],
#'                             rep(1,ncol(preyFAs[,-c(1,3)])-1),
#'                             rep(1,11))
#'      p.mat[i,] <- p.QFASA(my.seal,
#'                           MEANmeth(preyFAs[,-c(1,3)]),
#'                           rep(1,length(my.seal)),
#'                           2,
#'                           ext.fa=colnames(preyFAs[,-c(1:3)]))$`Diet Estimates`
#'  }
#'
#' # Can verify that average diet estimate of the 10 pseudo predators is close to
#' # "true" diet.
#'
#' round(apply(p.mat,2,mean),3)
#'

pseudo.pred <- function(diet, preybase, cal.vec, fat.vec, preysize=2) {

  preybase[,-1] <- preybase[,-1]/apply(preybase[,-1],1,sum)
  sort.preytype <- order(preybase[,1])
  preybase <- preybase[sort.preytype,]


  diet <- fat.vec * diet
  diet <- diet/sum(diet)


  I <- length(unique(preybase[,1]))
  nk.vec <- tapply(preybase[,1],preybase[,1],length)
  pred <- matrix(rep(NA,(ncol(preybase)-1)*I),nrow=ncol(preybase)-1,ncol=I,byrow=T)
  diet <- unlist(diet)

  cum <- cumsum(nk.vec)

  preybase.mat <- as.matrix(preybase[,-1])

  if (preysize==1) {

    # SELECTING ONE PREY SIGNATURE

    for (k in 1:I) {

      if (k==1) {


        ind <- sample(seq(1,cum[k],1),size=1)
        pred[,k] <- diet[k] * preybase.mat[ind, ]

      } else {


        # SELECTING ONE PREY SIGNATURE
        ind <- sample(seq((cum[k-1]+1),cum[k],1),size=1)
        pred[,k] <- diet[k] * preybase.mat[ind,]

      }

    } # end loop

  } else {

    # SAMPLING PREY DATA BASE

    for (k in 1:I) {

      if (k==1) {

        inds <- sample(seq(1,cum[k],1),size=nk.vec[k],replace=TRUE)
        pred[,k] <- diet[k] * apply(preybase.mat[inds, ],2,mean)

      } else {

        inds <- sample(seq((cum[k-1]+1),cum[k],1),size=nk.vec[k], replace=TRUE)
        pred[,k] <- diet[k] * apply(preybase.mat[inds,],2,mean)
      }

    } # end loop

  } # end outer if statement

  y.hat <- apply(pred,1,sum)*cal.vec

  y.hat <- y.hat/sum(y.hat)

  y.hat <- matrix(y.hat,nrow = 1)
  y.hat <- as.data.frame(y.hat)
  dimnames(y.hat)[[2]] <- dimnames(preybase[,-1])[[2]]

  return(y.hat)
}
