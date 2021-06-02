#'Simultaneous estimation of diet composition and calibration coefficients
#'
#'Computes the diet estimate for each predator in \emph{pred.sig} as well as an
#'overall estimate of the calibration coefficient vector.
#'
#' @param pred.sig matrix containing the FA signatures of the predator
#' @param prey.mat matrix containing a representative FA signature from each
#'                 prey group (usually the mean).  The first column must index
#'                 the prey group.
#' @param FC vector of fat content of length equal to the number of prey groups
#'           (or species)
#'
#' @details Starting values for the diet estimates are equal proportions
#'           and a vector of ones is used for the calibration coefficients.
#'
#' @return A list with components:
#' \item{diet.est}{A matrix of the diet estimates for each predator where each
#' row corresponds to a predator and the columns to prey species. The estimates
#' are expressed as proportions summing to one.}
#' \item{cc.est}{Estimated vector of calibration coefficients}
#'
#' @references Bromaghin, Jeffrey F., Budge, Suzanne M.,  Thiemann, Gregory and
#' Rode, Karyn D. (2017) Simultaneous estimation of the diet composition and
#' calibration coefficients with fatty acid signature data.
#' Ecology and Evolution, 7(16), 6103-6113
#'
#' @export
#'
#' @examples
#' ## Fatty Acids
#'data(FAset)
#'fa.set = as.vector(unlist(FAset))
#'
#' ## Predators
#'data(predatorFAs)
#'tombstone.info = predatorFAs[,1:4]
#'predator.matrix = predatorFAs[,5:(ncol(predatorFAs))]
#'npredators = nrow(predator.matrix)
#'
#' ## Need predator and prey to have same length
#'
#'predator.ext <- predator.matrix[fa.set]
#'predator.ext <- predator.ext/rowSums(predator.ext)
#'
#' ## Prey
#'data(preyFAs)
#'prey.sub=(preyFAs[,4:(ncol(preyFAs))])[fa.set]
#'prey.sub=prey.sub/apply(prey.sub,1,sum)
#'group=as.vector(preyFAs$Species)
#'prey.matrix=cbind(group,prey.sub)
#'prey.matrix=MEANmeth(prey.matrix)
#'
#' ## Fat Content
#'FC = preyFAs[,c(2,3)]
#'FC = as.vector(tapply(FC$lipid,FC$Species,mean,na.rm=TRUE))
#'
#'Q.sim <-p.sim.QFASA(predator.ext,prey.matrix,FC)
#' ## Average Diet Estimate
#'round(colMeans(Q.sim[[1]]),3)
#' ## Calibration Coefficients
#'Q.sim[[2]]

p.sim.QFASA <- function(pred.sig,prey.mat,FC=rep(1,nrow(prey.mat))){


  mult.rep <- function(x.vec) {

    no.zero <- sum(x.vec == 0)
    x.vec[x.vec == 0] <- 1e-05
    x.vec[x.vec > 0] <- (1 - no.zero * 1e-05) * x.vec[x.vec > 0]

    return(x.vec)

  }

  mean.geometric <- function(x) {
    D <- length(x)
    return(prod(x)^(1/D))
  }

  # constraint equation
  con.eqn <- function(pars){


    dietpars <- pars[1:(J*I)]
    dietpars <- matrix(dietpars, nrow=J, ncol= I, byrow=TRUE)
    sumpar <- rowSums(dietpars)


    return( c(sumpar, sum(pars[-(1:(J*I))]) ))
  }

  AIT.objective <- function(pars){


    pred <- t(apply(pred,1,mult.rep))

    diet.est <- pars[1:(J*I)]
    diet.est <- matrix(diet.est, nrow=J, ncol=I, byrow=TRUE)

    cal.est <- pars[-(1:(J*I))]

    prey <- Compositional::perturbation(prey,cal.est,oper="*")

    pred.hat <- diet.est%*%as.matrix(prey)

    pred.hat <- t(apply(pred.hat,1,mult.rep))

    geo.mn.1 <- as.vector(apply(pred,1,mean.geometric))
    geo.mn.2 <- as.vector(apply(pred.hat,1,mean.geometric))

    opt.fun <- sum(sqrt(rowSums( (log(pred/geo.mn.1)-log(pred.hat/geo.mn.2) )^2 )))

  }

  # number of predators
  J <- nrow(pred.sig)
  # number of FA
  K <- ncol(prey.mat)
  # number of prey types
  I <- nrow(prey.mat)

  if (K<I)
    { cat("Number of fatty acids is less than number of prey species.")
    return()
    }

  if (J < (K-1)/(K-I)) {

    cat("Degrees of freedom must equal or exceed number of parameters.")
    cat("Need a larger sample of predators.")
    return()

  }

  pred <- pred.sig
  prey <- prey.mat


  #initial diet proportions = inverse of number of prey types = 1/I
  start.diet <- rep((1/I),(J*I))
  #initial calibration coefficients = 1
  start.cal <- rep(1,K)

  p.all <- Rsolnp::solnp(pars=c(start.diet,start.cal),fun=AIT.objective,
                         eqfun = con.eqn, eqB = c(rep(1,length=J),K),
                         LB = c(rep(0,(J*I)), rep(0.02,K)) , UB = c(rep(1,(J*I)), rep(K,K)),
                         control = list(tol=1e-8,delta=1e-8))


  pars.est <- p.all[["pars"]]

  diet.est <- matrix(pars.est[1:(J*I)],ncol=I,byrow = T)

  if (is.matrix(FC)) {
    FC.mat <- FC
  }
  else {
    FC.mat <- matrix(rep(FC, J), byrow = T, nrow=J, ncol=I)
  }

  diet.est <- diet.est/FC.mat
  diet.est <- diet.est/rowSums(diet.est)

  colnames(diet.est) <- rownames(prey)
  cc.est <- pars.est[-c(1:(J*I))]
  names(cc.est) <- colnames(prey)

 return(list(diet.est,cc.est))
}
