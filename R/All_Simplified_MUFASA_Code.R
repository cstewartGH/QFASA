#' Returns simplified MLE diet estimates corresponding to a sample of predators.
#'
#' Computes the diet estimate for each predator in \emph{pred.mat} using the
#' simplified MLE method, without the use of random effects.
#'
#' @export
#' @param pred.mat matrix containing the FA signatures of the predators.
#' @param prey.mat matrix containing the FA signatures of the individual prey.
#'                 The first column must index the prey group.
#'                 \emph{prey.mat} is the prey database.
#' @param cal.mat matrix of calibration factors where the \emph{i} th
#'     column is to be used with the \emph{i} th predator. If modelling is to be
#'     done without calibration coefficients, simply pass a vector or matrix of
#'     ones.
#' @param FC vector of fat content of length equal to the number of prey groups
#'           or species.
#' @param ext.fa subset of fatty acids to be used to obtain diet estimates.
#' @import dplyr compositions
#' @return A list with components:
#' \item{Diet_Estimates}{A matrix of the diet estimates for each predator where
#' each row corresponds to a predator and the columns to prey species.
#' The estimates are expressed as proportions summing to one.}
#' \item{Var_Epsilon}{Optimized values of error variance.  See reference.}
#' \item{nll}{Negative log likelihood values.  As per \emph{solnp} documentation,
#'            \emph{nll} is a "vector of function values during optimization
#'  with last one the value at the optimal".}
#' @details The assumed model is similar to the MUFASA model but the random
#' effects are replaced by the prey species' sample means to speed up
#' computations. Unlike \emph{p.MUFASA}, this function does not require
#' integration and hence is much faster.
#' @references Steeves, Holly (2020) Maximum likelihood approach to diet
#' estimation and inference based on fatty acid signatures. PhD thesis available
#' at https://dalspace.library.dal.ca/handle/10222/80034.
#'
#' @examples
#'
#'##  This example takes some time to run.
#'## Please uncomment the code below to run.
#'
#'
#'#library(dplyr)
#'#library(compositions)
#'
#'
#'## Fatty Acids
#'#data(FAset)
#'#ext.fa <- as.vector(unlist(FAset))
#'
#'## Predators
#'#data(predatorFAs)
#'#pred.mat <- predatorFAs[, -c(1:4)]
#'#n.pred <- nrow(pred.mat)
#'
#'## Prey
#'#data(preyFAs)
#'#prey.mat <- preyFAs[, -c(1,3)]
#'
#'#FC = preyFAs[,c(2,3)]
#'#FC = as.vector(tapply(FC$lipid,FC$Species,mean,na.rm=TRUE))
#'
#'## Calibration Coefficients
#'#data(CC)
#'#cal.vec = CC[,2]
#'#cal.mat = replicate(n.pred, cal.vec)
#'#rownames(cal.mat) <- CC$FA
#'
#'
#'## Diet Estimates
#'#mle.est <- p.MLE(pred.mat, prey.mat, cal.mat, FC, ext.fa)
#'#mle.est$"Diet Estimates"

p.MLE <- function(pred.mat,
                  prey.mat,
                  cal.mat,
                  FC,
                  ext.fa) {

  # Adjust prey input
  colnames(prey.mat)[1] <- "Species"
  prey.mat <- prey.mat[order(prey.mat$Species), ]
  species.vec <- prey.mat[, 1]
  prey.mat <- prey.mat[, ext.fa]
  prey.mat <- multiplicativeReplacement(prey.mat)
  prey.mat <- prey.mat / apply(prey.mat, 1, sum)
  prey.mat <- data.frame(Species = species.vec, prey.mat)


  # Adjust predator input
  pred.mat <- pred.mat[, ext.fa]
  pred.mat <- pred.mat / apply(pred.mat, 1, sum)
  pred.mat <- multiplicativeReplacement(pred.mat)

  # Select FAs from CCs

  cal.mat <- cal.mat[ext.fa,]

  # Predator and prey counts
  I <- length(unique(prey.mat$Species))
  n.pred <- nrow(pred.mat)

  # Mean and covariance of ilr transformed prey
  prey.mat.t <- compositions::ilr(prey.mat[, -1])
  prey.mat.t <- data.frame(Species = species.vec, prey.mat.t)
  prey.means <- MEANmeth(prey.mat)
  prey.means.t <- MEANmeth(prey.mat.t)
  prey.var.t <- POOLVARmeth(prey.mat.t)$Spool

  # QFASA estimates for use in error start values
  Q <- p.QFASA(pred.mat,
               prey.means,
               cal.mat,
               dist.meas = 2,
               as.vector(FC),
               ext.fa = ext.fa
  )
  Q.est <- Q$"Diet Estimates"

  # Adjust predators by calibration coefficients
  if ((is.vector(cal.mat)) || (nrow(cal.mat) == 1) || (ncol(cal.mat) == 1)) {
    pred.mat <- t(t(pred.mat) / as.vector(unlist(cal.mat)))
    pred.mat <- pred.mat / apply(pred.mat, 1, sum)
    pred.mat <- as.matrix(pred.mat)
  } else {
    pred.mat <- pred.mat / t(cal.mat)
    pred.mat <- pred.mat / apply(pred.mat, 1, sum)
    pred.mat <- as.matrix(pred.mat)
  }

  # Ilr transform predator matrix
  pred.mat.t <- Compositional::alfa(pred.mat, 0)$aff

  # Get starting values for errors
  ers <- matrix(NA,
                nrow = 100 * n.pred,
                ncol = (ncol(prey.mat) - 1)
  )

  for (j in 1:100) {
    lb <- (j - 1) * n.pred + 1
    ub <- j * n.pred

    y.est <- matrix(NA,
                    nrow = n.pred,
                    ncol = (ncol(prey.mat) - 1)
    )

    for (i in 1:nrow(y.est)) {
      y.est[i, ] <- pseudo.pred.norm(
        prey.means.t,
        prey.var.t,
        Q.est[i, ]
      )
    }

    y.est <- compositions::acomp(y.est)
    pred.acomp <- compositions::acomp(pred.mat)

    ers[lb:ub, ] <- pred.acomp - y.est
  }

  ers <- compositions::acomp(ers)
  V <- compositions::ilrBase(D = length(ext.fa))
  G <- V %*% t(V)

  # Error start values
  sep.start <- diag(-(1 / 2)
                    * t(V)
                    %*% G
                    %*% compositions::variation(ers)
                    %*% G
                    %*% V)

  # Separate into four quantiles
  quan <- stats::quantile(sep.start)
  quan.start <- numeric(4)
  groupind <- numeric(length(sep.start))
  for (j in 1:4) {
    quan.start[j] <- mean(sep.start[sep.start >= quan[j] &
                                      sep.start <= quan[j + 1]])
    groupind[sep.start >= quan[j] &
               sep.start <= quan[j + 1]] <- j
  }

  # Function to optimize:
  opt.nll <- function(par.vec) {
    alpha <- par.vec[1:(length(par.vec) - length(quan.start))]
    alpha <- matrix(alpha, nrow = n.pred, byrow = FALSE)
    quan.start <- par.vec[(length(par.vec) - length(quan.start) + 1):
                            length(par.vec)]
    alpha <- cbind(alpha, 1 - apply(alpha, 1, sum))

    # Covariance matrix for error
    sep <- rep(1, length(ext.fa) - 1)
    for (l in 1:length(groupind)) {
      if (groupind[l] == 1) {
        sep[l] <- quan.start[1]
      } else if (groupind[l] == 2) {
        sep[l] <- quan.start[2]
      } else if (groupind[l] == 3) {
        sep[l] <- quan.start[3]
      } else {
        sep[l] <- quan.start[4]
      }
    }

    eta <- alpha %*% prey.means
    eta <- Compositional::alfa(eta, 0)$aff

    # Negative log likelihood
    nll <- rep(0, n.pred)
    for (j in 1:n.pred) {
      nll[j] <- -mvtnorm::dmvnorm(as.vector(pred.mat.t[j, ]),
                                  as.vector(eta[j, ]),
                                  sigma = diag(sep),
                                  log = TRUE
      )
    }

    return(sum(nll))
  }


  # Alpha sum constraint
  al.sum <- function(par.vec) {
    npars <- length(par.vec)
    alpha <- par.vec[1:(length(par.vec) - length(quan.start))]
    alpha <- matrix(alpha, nrow = n.pred, byrow = FALSE)
    return(c(apply(alpha, 1, sum)))
  }

  # Start values
  alpha.start <- Q.est[, 1:ncol(Q.est) - 1]
  par.vec <- c(as.vector(alpha.start), quan.start)

  LB <- c(
    rep(0, n.pred * (I - 1)),
    rep(1e-06, length(quan.start))
  )
  UB <- c(
    rep(1, n.pred * (I - 1)),
    rep(Inf, length(quan.start))
  )

  # Optimize
  optnt <- Rsolnp::solnp(
    pars = par.vec,
    fun = opt.nll,
    ineqfun = al.sum,
    ineqLB = rep(0, n.pred),
    ineqUB = rep(1, n.pred),
    LB = LB, UB = UB
    #control = list(tol = 1e-05)
  )

  # Get outputs
  L <- optnt$values
  alpha <- matrix(optnt$pars[1:((I - 1) * (n.pred))], nrow = n.pred)
  alpha <- cbind(alpha, 1 - apply(alpha, 1, sum))
  seps <- optnt$pars[((I - 1) * (n.pred) + 1):length(optnt$pars)]

  # Adjust for fat content
  if (is.matrix(FC)) {
    FC.mat <- FC
  } else {
    FC.mat <- matrix(rep(FC, n.pred),
                     byrow = T,
                     nrow = n.pred,
                     ncol = I
    )
  }

  alpha <- alpha / FC.mat
  alpha <- alpha / rowSums(alpha)
  colnames(alpha) <- colnames(Q.est)

  return(list(
    "Diet Estimates" = alpha,
    "Var Epsilon" = seps,
    "nll Value" = L
  ))
}

#' Multiplicative replacement of zeroes
#' @export
#' @param compositional.mat matrix containing compositional data (for example, FA signatures)
#'                          that may contain zeros.
#' @param delta imputed value
#' @return The compositional matrix with zeros modified.
multiplicativeReplacement <- function(compositional.mat,
                                      delta = 1e-05) {
  for (i in 1:nrow(compositional.mat)) {
    # Count the zeroes in the row
    rep.row <- compositional.mat[i, ]
    zero.count <- sum(rep.row == 0)

    # Index the elements
    zero <- (rep.row == 0)
    nonzero <- (rep.row > 0)

    # Multiplicative replacement
    rep.row[nonzero] <- (1 - zero.count * delta) * rep.row[nonzero]
    rep.row[zero] <- delta
    compositional.mat[i, ] <- rep.row
  }

  return(compositional.mat)
}


