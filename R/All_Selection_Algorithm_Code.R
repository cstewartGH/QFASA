#' Returns diet estimates corresponding to a sample of predators based on a
#' forward selection algorithm that chooses the prey species to be included in
#' the modelling.
#'
#' @export
#' @param pred.mat matrix containing the FA signatures of the predators, where
#'                 each row corresponds to a predator and each column to a FA.
#' @param prey.mat data frame containing the FA signatures of the prey, where
#'                 each row corresponds to a single individual prey.
#'                 The first column must index the prey group, while the
#'                 remaining columns correspond to FAs.
#' @param cal.vec	 numeric vector of calibration coefficients corresponding to
#'                 the FAs contained in the modelling subset ext.fa. A vector of
#'                 ones in the length of ext.fa may be used for modelling
#'                 without calibration coefficients.
#' @param FC	     numeric vector of the average lipid contents for each prey
#'                 group, in the order of the \emph{alphabetized} prey groups.
#'                 vector of ones equal in length to the number of prey groups
#'                 may be used for modelling without adjustment for fat content.
#' @param ext.fa	 character vector containing the subset of FAs to be used in
#'                 modelling.
#' @param k        scaling factor to be used in calculating the value of the
#'                 information criterion (IC). The default value of 2
#'                 corresponds to the Akaike Information Criterion. For a sample
#'                 of size n, k = log(n) corresponds to the Bayesian Information
#'                 Criterion. As this factor is numeric, values corresponding to
#'                 other IC may be used freely.
#' @param min.spec optional integer value specifying the minimum final model
#'                 size for forward selection. By default, forward selection
#'                 will add species to the model until the value of the chosen
#'                 IC ceases to improve. If this parameter is increased, forward
#'                 selection will add the best available species up to the
#'                 specified minimum model size before continuing with the
#'                 default selection process.
#' @param starting.spec optional character vector specifying the starting
#'                      species for the forward selection algorithm. Where known,
#'                      two or more species may be specified to ensure their
#'                      inclusion in the final model, reducing computation times
#'                      for the algorithm.  The default is NULL.
#' @param silence if true, additional information is printed. Default is false.
#' @import dplyr compositions
#' @return A list with components:
#' \item{Diet_Estimates}{A matrix of the diet estimates for each predator where
#' each row corresponds to a predator and each column to a prey species. The
#' estimates are expressed as proportions summing to one.}
#' \item{Selection_Order}{A data frame summarizing each step of the algorithm,
#' giving the order of species selection and the corresponding IC values.}
#' \item{Selection_Tables}{A list containing a data frame for each step of the
#' selection process, providing the IC values associated with adding any one
#' candidate species at that step.}
#'
#' @details The function uses a forward selection algorithm and the simplified
#' MLE method to choose the prey species to be included in the model and then
#' returns the diet estimates corresponding to these species.
#' @seealso \emph{backward.elimination()}
#' @examples
#'
#'  ## This example takes some time to run.
#'  ## Please uncomment code below to run.
#'
#'#library(dplyr)
#'#library(compositions)
#'## Package data: FAs
#'#data(FAset)
#'#fa.set = as.vector(unlist(FAset))
#'
#'## Package data: Prey
#'#data(preyFAs)
#'#prey.sub=(preyFAs[,4:(ncol(preyFAs))])[fa.set]
#'#prey.sub=prey.sub/apply(prey.sub,1,sum)
#'#group=as.vector(preyFAs$Species)
#'#prey.sub = cbind(group,prey.sub)
#'#sort.preytype <- order(prey.sub[, 1])
#'#prey.matrix <- prey.sub[sort.preytype,]
#'
#'## Package data: Predators
#'#data(predatorFAs)
#'#tombstone.info = predatorFAs[,1:4]
#'#predator.matrix = predatorFAs[,5:(ncol(predatorFAs))]
#'#npredators = nrow(predator.matrix)
#'
#'## Package data: Fat content
#'#FC = preyFAs[,c(2,3)]
#'#FC = as.vector(tapply(FC$lipid,FC$Species,mean,na.rm=TRUE))
#'
#'## Package data: Calibration coefficients
#'#data(CC)
#'#cal.vec = CC[,2]
#'#cal.mat = replicate(npredators, cal.vec)
#'#rownames(cal.mat) <- CC$FA
#'#names(cal.vec) <- rownames(cal.mat)
#'
#'## QFASA (KL)
#'#sample.qfasa <- p.QFASA(predator.matrix,MEANmeth(prey.matrix),cal.mat,
#'#dist.meas = 1,gamma=1,FC,
#'#start.val = rep(1,nrow(MEANmeth(prey.matrix))),fa.set)
#'
#'## Forward Selection
#'#sample.fs <- forward.selection(predator.matrix,prey.matrix,cal.vec,FC,fa.set,
#'#min.spec = 5,starting.spec = c("capelin", "herring"))
#'## Output
#'#fs.estimates <- sample.fs$`Diet Estimates`

forward.selection <- function(pred.mat,
                              prey.mat,
                              cal.vec,
                              FC,
                              ext.fa,
                              k = 2,
                              min.spec = 5,
                              starting.spec = NULL,
                              silence = FALSE) {

  # Arrange prey matrix
  prey.mat <- cbind(prey.mat[, 1], prey.mat[, ext.fa])
  colnames(prey.mat) <- c("Species", ext.fa)
  prey.mat <- prey.mat[order(prey.mat$Species), ]

  # Calibration coefficients
  if (length(cal.vec) > length(ext.fa)) {
    if (identical(cal.vec, unname(cal.vec)) == FALSE) {
      cal.vec <- cal.vec[ext.fa]
    }
  }

  # Save full prey library
  species.list <- sort(unique(prey.mat$Species))
  full.prey.mat <- prey.mat

  # Initialize output object
  aic.tables <- list()

  # Handle starting species
  if (!is.null(starting.spec)) {
    prey.mat       <- prey.mat[prey.mat$Species%in%starting.spec, ]
    species.index  <- which(species.list%in%starting.spec)
    starting.model <- evaluateModel(pred.mat,
                                    prey.mat,
                                    cal.vec,
                                    FC[species.index],
                                    ext.fa,
                                    k)
    aic.old         <- starting.model$'AIC Value'
    it.index        <- 0
  } else {
    starting.model <- bestSubset(pred.mat,
                                 prey.mat,
                                 cal.vec,
                                 FC,
                                 ext.fa,
                                 k,
                                 silence)
    aic.old         <- starting.model$'AIC Value'
    prey.mat        <- starting.model$'Prey Mat'
    aic.tables[[1]] <- starting.model$'Combinations'
    it.index        <- 1
    species.index   <- which(species.list%in%starting.model$'Best Subset')
  }

  # Record selection order
  selection.order <- data.frame("Model Size" =
                                  c(2, 2),
                                "Species Added" =
                                  c(species.list[species.index[1]],
                                    species.list[species.index[2]]),
                                "AIC" =
                                  c(aic.old, aic.old),
                                check.names = FALSE)

  # Evaluate models of the next size up
  new.models <- bestForwardModel(pred.mat,
                                 full.prey.mat,
                                 prey.mat,
                                 cal.vec,
                                 FC,
                                 ext.fa,
                                 k,
                                 silence)

  # Update results
  it.index               <- it.index + 1
  aic.new                <- new.models[[1]]$'AIC Value'
  aic.tables[[it.index]] <- new.models$'AIC Values'
  estimates              <- new.models$'Best Model'[[1]]
  model.size             <- length(unique(prey.mat$Species)) + 1
  species.added     <- new.models[[1]]$'Species'
  species.3         <- which(species.list == species.added)
  selection.details <- data.frame("Model Size" =
                                    model.size,
                                  "Species Added" =
                                    species.added,
                                  "AIC" = aic.new,
                                  check.names = FALSE)
  selection.order   <- rbind(selection.order, selection.details)

  # Iterate until improvement stops
  while (aic.new < aic.old | model.size <= min.spec) {
    aic.old  <- aic.new
    prey.mat <- new.models$'Best Model'[[2]]

    # Evaluate models of the next size down
    new.models <- bestForwardModel(pred.mat,
                                   full.prey.mat,
                                   prey.mat,
                                   cal.vec,
                                   FC,
                                   ext.fa,
                                   k,
                                   silence)

    # Find the best model and update results
    aic.new   <- new.models[[1]]$'AIC Value'

    # Update results
    model.size             <- model.size + 1
    it.index               <- it.index + 1
    aic.tables[[it.index]] <- new.models$'AIC Values'

    if (aic.new < aic.old | model.size <= min.spec) {
      estimates     <- new.models$'Best Model'[[1]]
      species.added <- new.models[[1]]$'Species'
      species.index <- which(species.list == species.added)
      selection.details <- data.frame("Model Size" =
                                        model.size,
                                      "Species Added" =
                                        species.added,
                                      "AIC" = aic.new,
                                      check.names = FALSE)
      selection.order   <- rbind(selection.order, selection.details)
    }
  }

  # Gather results
  estimates         <- zeroEstimates(estimates, full.prey.mat)
  results           <- list("Diet Estimates"    = estimates,
                            "Selection Order"   = selection.order,
                            "Selection Tables"  = aic.tables)
  return(results)
}

#' Returns diet estimates corresponding to a sample of predators based on a
#' backward elimination algorithm that chooses the prey species to be included in
#' the modelling.
#'
#' @export
#' @param pred.mat matrix containing the FA signatures of the predators, where
#'                 each row corresponds to a predator and each column to a FA.
#' @param prey.mat data frame containing the FA signatures of the prey, where
#'                 each row corresponds to a single individual prey.
#'                 The first column must index the prey group, while the
#'                 remaining columns correspond to FAs.
#' @param cal.vec	 numeric vector of calibration coefficients corresponding to
#'                 the FAs contained in the modelling subset ext.fa. A vector of
#'                 ones in the length of ext.fa may be used for modelling
#'                 without calibration coefficients.
#' @param FC	     numeric vector of the average lipid contents for each prey
#'                 group, in the order of the \emph{alphabetized} prey groups.
#'                 vector of ones equal in length to the number of prey groups
#'                 may be used for modelling without adjustment for fat content.
#' @param ext.fa	 character vector containing the subset of FAs to be used in
#'                 modelling.
#' @param k        scaling factor to be used in calculating the value of the
#'                 information criterion (IC). The default value of 2
#'                 corresponds to the Akaike Information Criterion. For a sample
#'                 of size n, k = log(n) corresponds to the Bayesian Information
#'                 Criterion. As this factor is numeric, values corresponding to
#'                 other IC may be used freely.
#' @param cutoff   numeric proportion to be used as the threshold for the
#'                 candidacy for removal of a species in backward
#'                 elimination. If initial diet estimates for any
#'                 individual predator find a species to be present in
#'                 proportions greater than this threshold, that species
#'                 will not be considered for removal from the model. This
#'                 reduces computation times and safeguards against the
#'                 removal of species which may be present in the diets of
#'                 few predators. All species are considered candidates for
#'                 removal at a value of 1, while lower values are more
#'                 conservative. The default value is 0.1.
#' @param silence if true, additional information is printed. Default is false.
#' @return A list with components:
#' \item{Diet_Estimates}{A matrix of the diet estimates for each predator where
#' each row corresponds to a predator and each column to a prey species. The
#' estimates are expressed as proportions summing to one.}
#' \item{Selection_Order}{A data frame summarizing each step of the algorithm,
#' giving the order of species removal and the corresponding IC values.}
#' \item{Selection_Tables}{A list containing a data frame for each step of the
#' selection process, providing the IC values associated with removing any one
#' candidate species at that step.}
#'
#' @details The function uses a backward elimination algorithm and the simplified
#' MLE method to choose the prey species to be included in the model and then
#' returns the diet estimates corresponding to these species.
#' @seealso \emph{forward.selection()}
#' @examples
#'
#'  ## This example takes some time to run.
#'  ## Please uncomment code below to run.
#'
#'#library(dplyr)
#'#library(compositions)
#'## Package data: FAs
#'#data(FAset)
#'#fa.set = as.vector(unlist(FAset))
#'
#'## Package data: Prey
#'#data(preyFAs)
#'#prey.sub=(preyFAs[,4:(ncol(preyFAs))])[fa.set]
#'#prey.sub=prey.sub/apply(prey.sub,1,sum)
#'#group=as.vector(preyFAs$Species)
#'#prey.sub = cbind(group,prey.sub)
#'#sort.preytype <- order(prey.sub[, 1])
#'#prey.matrix <- prey.sub[sort.preytype,]
#'
#'## Package data: Predators
#'#data(predatorFAs)
#'#tombstone.info = predatorFAs[,1:4]
#'#predator.matrix = predatorFAs[,5:(ncol(predatorFAs))]
#'#npredators = nrow(predator.matrix)
#'
#'## Package data: Fat content
#'#FC = preyFAs[,c(2,3)]
#'#FC = as.vector(tapply(FC$lipid,FC$Species,mean,na.rm=TRUE))
#'
#'## Package data: Calibration coefficients
#'#data(CC)
#'#cal.vec = CC[,2]
#'#cal.mat = replicate(npredators, cal.vec)
#'#rownames(cal.mat) <- CC$FA
#'#names(cal.vec) <- rownames(cal.mat)
#'
#'## QFASA (KL)
#'#sample.qfasa <- p.QFASA(predator.matrix,MEANmeth(prey.matrix),cal.mat,
#'#dist.meas = 1,gamma=1,FC,
#'#start.val = rep(1,nrow(MEANmeth(prey.matrix))),fa.set)
#'
#'## Backward Elimination
#'#sample.be <- backward.elimination(predator.matrix, prey.matrix, cal.vec,FC,
#'#fa.set,cutoff = 0.01)
#'
#'## Output
#'#be.estimates <- sample.be$`Diet Estimates`

backward.elimination <- function(pred.mat,
                                 prey.mat,
                                 cal.vec,
                                 FC,
                                 ext.fa,
                                 k = 2,
                                 cutoff = 0.1,
                                 silence = FALSE) {

  # Sort prey mat! (new)
  prey.mat <- cbind(prey.mat[, 1], prey.mat[, ext.fa])
  colnames(prey.mat) <- c("Species", ext.fa)
  prey.mat <- prey.mat[order(prey.mat$Species), ]

  # Check calibration vector
  if (length(cal.vec) > length(ext.fa)) {
    cal.vec <- cal.vec[ext.fa]
  }

  # Set up
  species.list  <- sort(unique(prey.mat$Species))
  full.prey.mat <- prey.mat
  aic.tables    <- list()

  # Evaluate full model
  full.model <- evaluateModel(pred.mat,
                              prey.mat,
                              cal.vec,
                              FC,
                              ext.fa,
                              k)
  aic.old <- full.model$'AIC Value'

  # Record selection order
  selection.order <- data.frame("Model Size" =
                                  length(unique(prey.mat$Species)),
                                "Species Dropped" = NA,
                                "AIC" = aic.old,
                                check.names = FALSE)

  # Evaluate models of the next size down
  new.models <- bestBackwardModel(pred.mat,
                                  prey.mat,
                                  cal.vec,
                                  FC,
                                  ext.fa,
                                  full.model$Estimates,
                                  k,
                                  cutoff,
                                  silence)

  # Update results
  if (is.na(new.models[1])) {
    estimates <- full.model$Estimates
    results   <- list("Diet Estimates" = estimates)
  } else {
    # Find lowest AIC of new models
    aic.new   <- new.models[[1]]$'AIC Value'

    # Assign previous estimates if AIC hasn't improved
    if (aic.old < aic.new) {
      estimates <- full.model$Estimates
    }

    # Iterate until improvement stops
    while (aic.new < aic.old) {
      # Update model information
      estimates <- new.models$'Best Model'[[1]]
      prey.mat  <- new.models$'Best Model'[[2]]
      FC   <- new.models$'Best Model'[[3]]

      # Update output
      aic.tables <- c(aic.tables,
                      list(stats::na.omit(new.models$'AIC Values')))
      model.size        <- length(unique(prey.mat$Species))
      dropped.species   <- new.models[[1]]$'Species'
      species.index     <- which(species.list == dropped.species)
      selection.details <- data.frame("Model Size" =
                                        model.size,
                                      "Species Dropped" =
                                        dropped.species,
                                      "AIC" = aic.new,
                                      check.names = FALSE)
      selection.order   <- rbind(selection.order, selection.details)

      # Evaluate models of the next size down
      aic.old <- aic.new
      new.models <- bestBackwardModel(pred.mat,
                                      prey.mat,
                                      cal.vec,
                                      FC,
                                      ext.fa,
                                      estimates,
                                      k,
                                      cutoff,
                                      silence)

      # Find lowest AIC of new models
      if (!is.na(new.models[1])) {
        aic.new   <- new.models[[1]]$'AIC Value'
      } else {
        break
      }
    }

    # Gather results
    estimates         <- zeroEstimates(estimates, full.prey.mat)
    results           <- list("Diet Estimates"    = estimates,
                              "Selection Order" = selection.details,
                              "Selection Tables" = aic.tables)
  }

  return(results)
}

#' Get simplified MLE estimates
#' @keywords internal
likelihoodEstimates <- function(pred.mat,
                                prey.mat,
                                cal.mat,
                                fat.vec,
                                ext.fa) {

  # Adjust prey input
  colnames(prey.mat)[1] <- "Species"
  prey.mat              <- prey.mat[order(prey.mat$Species), ]
  species.vec           <- prey.mat[, 1]
  prey.mat              <- prey.mat[, ext.fa]
  prey.mat              <- multiplicativeReplacement(prey.mat)
  prey.mat              <- prey.mat / apply(prey.mat, 1, sum)
  prey.mat              <- data.frame(Species = species.vec, prey.mat)

  # Adjust predator input
  pred.mat <- pred.mat[, ext.fa]
  pred.mat <- pred.mat / apply(pred.mat, 1, sum)
  pred.mat <- multiplicativeReplacement(pred.mat)

  # Predator and prey counts
  I      <- length(unique(prey.mat$Species))
  n.pred <- nrow(pred.mat)

  # Mean and covariance of ilr transformed prey
  prey.mat.t   <- compositions::ilr(prey.mat[, -1])
  prey.mat.t   <- data.frame(Species = species.vec, prey.mat.t)
  prey.means   <- MEANmeth(prey.mat)
  prey.means.t <- MEANmeth(prey.mat.t)
  prey.var.t   <- POOLVARmeth(prey.mat.t)$Spool

  # QFASA estimates for use in error start values
  Q     <- p.QFASA(pred.mat,
                   prey.means,
                   cal.mat,
                   dist.meas = 2,
                   as.vector(fat.vec),
                   ext.fa = ext.fa)
  Q.est <- Q$'Diet Estimates'

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
  pred.mat.t <- Compositional::alfa(pred.mat,0)$aff

  # Get starting values for errors
  ers <- matrix(NA,
                nrow = 50 * n.pred,
                ncol = (ncol(prey.mat) - 1))

  for (j in 1:50) {
    lb <- (j - 1) * n.pred + 1
    ub <- j * n.pred

    y.est <- matrix(NA,
                    nrow = n.pred,
                    ncol = (ncol(prey.mat) - 1))

    for (i in 1:nrow(y.est)) {
      y.est[i, ] <- pseudo.pred.norm(prey.means.t,
                                     prey.var.t,
                                     Q.est[i, ])
    }

    y.est      <- compositions::acomp(y.est)
    pred.acomp <- compositions::acomp(pred.mat)

    ers[lb:ub,] = pred.acomp - y.est
  }

  ers <- compositions::acomp(ers)
  V   <- compositions::ilrBase(D = length(ext.fa))
  G   <- V %*% t(V)

  # Error start values
  sep.start <- diag(- (1 / 2)
                    * t(V)
                    %*% G
                    %*% compositions::variation(ers)
                    %*% G
                    %*% V)

  # Separate into four quantiles
  quan       <- stats::quantile(sep.start)
  quan.start <- numeric(4)
  groupind   <- numeric(length(sep.start))
  for (j in 1:4) {
    quan.start[j] <- mean(sep.start[sep.start >= quan[j]
                                    & sep.start <= quan[j + 1]])
    groupind[sep.start >= quan[j]
             & sep.start <= quan[j + 1]] <- j
  }

  # Function to optimize:
  opt.nll <- function(par.vec) {
    alpha      <- par.vec[1:(length(par.vec) - length(quan.start))]
    alpha      <- matrix(alpha, nrow = n.pred, byrow = FALSE)
    quan.start <- par.vec[(length(par.vec) - length(quan.start)+1):
                            length(par.vec)]
    alpha      <- cbind(alpha, 1 - apply(alpha, 1, sum))

    # Covariance matrix for error
    sep <- rep(1, length(ext.fa) - 1)
    for (l in 1:length(groupind)) {
      if        (groupind[l] == 1) {
        sep[l] = quan.start[1]
      } else if (groupind[l] == 2) {
        sep[l] = quan.start[2]
      } else if (groupind[l] ==3) {
        sep[l] = quan.start[3]
      } else {
        sep[l] = quan.start[4]
      }
    }

    eta <- alpha %*% prey.means
    eta <- Compositional::alfa(eta, 0)$aff

    # Negative log likelihood
    nll <- rep(0, n.pred)
    for(j in 1:n.pred) {
      nll[j] <- -mvtnorm::dmvnorm(as.vector(pred.mat.t[j, ]),
                                  as.vector(eta[j, ]),
                                  sigma = diag(sep),
                                  log = TRUE)
    }

    return(sum(nll))
  }


  # Alpha sum constraint
  al.sum <- function(par.vec)
  {
    npars <- length(par.vec)
    alpha <- par.vec[1:(length(par.vec) - length(quan.start))]
    alpha <- matrix(alpha, nrow = n.pred, byrow = FALSE)
    return(c(apply(alpha, 1, sum)))
  }

  # Start values
  alpha.start <- Q.est[, 1:ncol(Q.est) - 1]
  par.vec     <- c(as.vector(alpha.start), quan.start)

  LB = c(rep(0, n.pred * (I-1)),
         rep(1e-06, length(quan.start)))
  UB = c(rep(1, n.pred * (I-1)),
         rep(Inf, length(quan.start)))

  # Optimize
  optnt <- Rsolnp::solnp(pars = par.vec,
                         fun = opt.nll,
                         ineqfun = al.sum,
                         ineqLB = rep(0, n.pred),
                         ineqUB = rep(1, n.pred),
                         LB = LB, UB = UB,
                         control = list(tol = 1e-05))

  # Get outputs
  L     <- optnt$values
  alpha <- matrix(optnt$pars[1:((I - 1) * (n.pred))], nrow = n.pred)
  alpha <- cbind(alpha, 1 - apply(alpha, 1, sum))
  seps  <- optnt$pars[((I - 1) * (n.pred) + 1):length(optnt$pars)]

  # Adjust for fat content
  if (is.matrix(fat.vec)) {
    fat.vec.mat <- fat.vec
  } else {
    fat.vec.mat <- matrix(rep(fat.vec, n.pred),
                          byrow = T,
                          nrow = n.pred,
                          ncol = I)
  }

  alpha           <- alpha/fat.vec.mat
  alpha           <- alpha/rowSums(alpha)
  colnames(alpha) <- colnames(Q.est)

  return(list("Diet Estimates" = alpha,
              "Var Epsilon" = seps,
              "nll Value" = L))
}

#' Get nll and IC values for a given model (note that the model is determined by
#' the species included in prey.mat - prey.mat is adjusted to add or remove
#' species before entering this function
#'
#' @keywords internal

evaluateModel <- function(pred.mat,
                          prey.mat,
                          cal.vec,
                          fat.vec,
                          ext.fa,
                          k = 2) {

  # Parameters
  I       <- length(unique(prey.mat$Species))
  n.pred  <- nrow(pred.mat)
  cal.mat <- replicate(n.pred, cal.vec)

  # Evaluate model
  model.output <- likelihoodEstimates(pred.mat,
                                      prey.mat,
                                      cal.mat,
                                      fat.vec,
                                      ext.fa)
  estimates    <- model.output$'Diet Estimates'

  # Calculate AIC value
  nll.value    <- utils::tail(model.output$'nll Value', n=1)
  penalty.term <- k * I * n.pred
  aic.value    <- (2 * nll.value) + penalty.term

  # Output
  return(list('nll Value' = nll.value,
              'AIC Value' = aic.value,
              'Estimates' = estimates))
}

#' Reintroduce excluded species to diet estimates as forced zeroes
#' @keywords internal

zeroEstimates <- function(estimates,
                          prey.mat) {

  '%notin%'      <- Negate('%in%')

  # Subset the prey species
  species.list   <- sort(unique(prey.mat$Species))
  used.species   <- sort(colnames(estimates))
  used.index     <- which(species.list%in%used.species)
  unused.species <- sort(species.list[species.list%notin%used.species])
  unused.index   <- which(species.list%in%unused.species)

  # Add zero estimates for unused species
  zero.vec <- rep(0, nrow(estimates))
  new.ind  <- vector()
  if (length(unused.species) > 0) {
    for (i in 1:length(unused.species)) {
      ind                 <- unused.index[i]
      estimates           <- cbind(estimates, zero.vec)
      new.ind             <- c(new.ind, species.list[ind])
      colnames(estimates) <- c(used.species, new.ind)
    }
  }

  estimates <- estimates[, species.list]

  return(estimates)
}

#' Find the pair of starting species with the best (lowest) IC value among all
#' possible pairs of species in prey.mat Used in forward selection when no
#' starting species are specified
#' @keywords internal

bestSubset <- function(pred.mat,
                       prey.mat,
                       cal.vec,
                       fat.vec,
                       ext.fa,
                       k = 2,
                       silence = FALSE) {

  # Index combinations
  species.list       <- sort(unique(prey.mat$Species))
  I                  <- length(species.list)
  combinations.index <- utils::combn(1:I, 2)
  n.comb             <- choose(I, 2)

  # Print update
  if (silence == FALSE) {
    start.time <- Sys.time()
    print(paste0(n.comb, " models of size 2 to be evaluated."))
  }

  # Evaluate the models
  model.output <- apply(combinations.index,
                        2,
                        evaluateCombination,
                        species.list,
                        pred.mat,
                        prey.mat,
                        cal.vec,
                        fat.vec,
                        ext.fa,
                        k)

  # Table of results
  results.table             <- as.data.frame(t(model.output))
  colnames(results.table)   <- c("Species 1",
                                 "Species 2",
                                 "nll Value",
                                 "AIC Value")
  results.table$'nll Value' <- as.numeric(results.table$'nll Value')
  results.table$'AIC Value' <- as.numeric(results.table$'AIC Value')

  # Find the best subset
  aic.value <- min(results.table$'AIC Value', na.rm = TRUE)
  aic.index <- which(results.table$'AIC Value' == aic.value)
  species   <- results.table[aic.index, 1:2]
  prey.mat  <- prey.mat[prey.mat$Species%in%species, ]

  # Print update
  if (silence == FALSE) {
    runtime <- Sys.time() - start.time
    units(runtime) <- "mins"
    print(paste0("Finished evaluating ", n.comb, " combinations in ",
                 runtime, " minutes."))
    print(paste0("The best combination is ", tolower(species[1]),
                 " and ", tolower(species[2]), " with an AIC value of ",
                 aic.value, "."))
  }

  # Gather results
  results <- list("Combinations" = results.table,
                  "Best Subset"  = species,
                  "AIC Value"    = aic.value,
                  "Prey Mat"     = prey.mat)

  return(results)
}

#' Get nll and IC values for a pair of species (this is used in forward
#' selection if the starting species are not specified; used on a list of all
#' possible combinations of two species)
#' @keywords internal

evaluateCombination <- function(species.index,
                                species.list,
                                pred.mat,
                                prey.mat,
                                cal.vec,
                                fat.vec,
                                ext.fa,
                                k = 2) {

  # Identify species in combination
  species.1 <- species.list[species.index[1]]
  species.2 <- species.list[species.index[2]]
  species   <- c(species.1, species.2)

  # Prey mat and fat content
  prey.mat <- prey.mat[prey.mat$Species%in%species, ]
  fat.vec  <- fat.vec[species.index]

  # Evaluate the model
  model.output <- evaluateModel(pred.mat,
                                prey.mat,
                                cal.vec,
                                fat.vec,
                                ext.fa,
                                k)

  # Gather results
  results <- c(species.1,
               species.2,
               model.output$'nll Value',
               model.output$'AIC Value')

  return(results)
}

#' Finds the species combination with the best (lowest) IC value achievable by
#' dropping a single species from the prey.mat. Used for each iteration of
#' backward elimination.
#' @keywords internal

bestBackwardModel <- function(pred.mat,
                              prey.mat,
                              cal.vec,
                              fat.vec,
                              ext.fa,
                              starting.estimates,
                              k = 2,
                              cutoff = 0.01,
                              silence = FALSE) {

  # Count and list prey species
  species.list <- sort(unique(prey.mat$Species))
  I            <- length(species.list)

  # Find number of models to be evaluated
  max.proportion <- apply(starting.estimates,
                          2,
                          max,
                          na.rm = TRUE)
  model.count    <- sum(max.proportion < cutoff)

  # Print update
  if (silence == FALSE) {
    start.time   <- Sys.time()
    model.number <- 0
    print(paste0(model.count, " models of size ",
                 I-1, " to be assessed."))
  }

  if (model.count > 0) {
    # Initialize outputs
    estimates     <- list()
    results.table <- data.frame(matrix(NA,
                                       nrow = I,
                                       ncol = 3))
    colnames(results.table) <- c("Dropped Species",
                                 "nll Value",
                                 "AIC Value")

    # Evaluate all candidate models
    for (i in 1:I) {
      if (max.proportion[i] < cutoff) {
        # Update arguments
        fat.vec.new  <- fat.vec[-i]
        prey.mat.new <- prey.mat[prey.mat$Species != species.list[i], ]

        # Evaluate the model omitting species i
        model.output <- evaluateModel(pred.mat,
                                      prey.mat.new,
                                      cal.vec,
                                      fat.vec.new,
                                      ext.fa,
                                      k)

        # Update outputs
        results.table[i, 1] <- species.list[i]
        results.table[i, 2] <- model.output$'nll Value'
        results.table[i, 3] <- model.output$'AIC Value'
        estimates[[i]]      <- model.output$'Estimates'

        # Print update
        if (silence == FALSE) {
          model.number <- model.number + 1
          print(paste0("Finished evaluating model ",
                       model.number, " of ", model.count,
                       " models of size ", I-1, "."))
        }
      }
    }

    # Find the best model
    aic.value     <- min(results.table$'AIC Value', na.rm = TRUE)
    aic.index     <- which(results.table$'AIC Value' == aic.value)
    best.species  <- results.table[aic.index, 1]
    species.index <- which(species.list == best.species)

    # Update arguments
    prey.mat <- prey.mat[prey.mat$Species != best.species, ]
    fat.vec  <- fat.vec[-species.index]

    # Gather results
    best.model <- list("Estimates" = estimates[[aic.index]],
                       "Prey Mat"  = prey.mat,
                       "Fat Vec"   = fat.vec,
                       "Species"   = best.species,
                       "AIC Value" = aic.value)
    results <- list("Best Model"   = best.model,
                    "AIC Values"   = stats::na.omit(results.table),
                    "Estimates"    = estimates)

    # Print update
    if (silence == FALSE) {
      runtime <- Sys.time() - start.time
      units(runtime) <- "mins"
      print(stats::na.omit(results.table))
      print(paste0("Finished evaluating ", model.count,
                   " models of size ", I-1, " in ",
                   runtime, " minutes."))
    }

  } else {
    results <- NA
    print("There are no candidates for removal.")
  }

  return(results)
}

#'Finds the species combination with the best# (lowest) IC value achievable by
#'adding a single species from full.prey.mat. Used for each iteration of
#'forward selection.
#'@keywords internal

bestForwardModel <- function(pred.mat,
                             full.prey.mat,
                             prey.mat,
                             cal.vec,
                             fat.vec,
                             ext.fa,
                             k = 2,
                             silence = FALSE) {

  '%notin%' <- Negate('%in%')

  # Subset prey species
  species.list   <- sort(unique(full.prey.mat$Species))
  I              <- length(species.list)
  used.species   <- sort(unique(prey.mat$Species))
  used.index     <- which(species.list%in%used.species)
  unused.species <- species.list[species.list%notin%used.species]
  unused.index   <- which(species.list%in%unused.species)

  # Initialize results
  estimates               <- list()
  results.table           <- data.frame(matrix(NA,
                                               nrow = length(unused.species),
                                               ncol = 3))
  colnames(results.table) <- c("Species Added",
                               "nll Value",
                               "AIC Value")

  # Print update
  if (silence == FALSE) {
    start.time   <- Sys.time()
    model.number <- 0
    print(paste0(length(unused.species), " models of size ",
                 length(used.species) + 1, " to be evaluated."))
  }

  # Check each unused species
  for (i in 1:length(unused.species)) {
    species.index <- sort(c(used.index, unused.index[i]))
    species       <- species.list[species.index]
    new.fat.vec   <- fat.vec[species.index]
    new.prey.mat  <- full.prey.mat[full.prey.mat$Species%in%species, ]

    # Evaluate the model
    model.output  <- evaluateModel(pred.mat,
                                   new.prey.mat,
                                   cal.vec,
                                   new.fat.vec,
                                   ext.fa,
                                   k)

    # Update results
    results.table[i, ] <- c(unused.species[i],
                            model.output$'nll Value',
                            model.output$'AIC Value')
    estimates[[i]]     <- model.output$'Estimates'

    # Print update
    if (silence == FALSE) {
      model.number <- model.number + 1
      print(paste0("Finished evaluating model ", model.number,
                   " of ", length(unused.species),
                   " models of size ", length(used.species) + 1,
                   "."))
    }
  }

  # Find the best model
  results.table$'AIC Value' <- as.numeric(results.table$'AIC Value')
  aic.value     <- min(results.table$'AIC Value', na.rm = TRUE)
  aic.index     <- which(results.table$'AIC Value' == aic.value)
  species       <- sort(c(used.species, results.table[aic.index, 1]))
  species.index <- which(species.list%in%species)
  new.prey.mat  <- full.prey.mat[full.prey.mat$Species%in%species, ]
  new.fat.vec   <- fat.vec[species.index]

  # Print update
  if (silence == FALSE) {
    runtime        <- Sys.time() - start.time
    units(runtime) <- "mins"
    print(results.table)
    print(paste0("Finished evaluating ", length(unused.species),
                 " models of size ", length(used.species) + 1,
                 " in ", runtime, " minutes."))
    print(paste0("The best species to add is ",
                 tolower(results.table[aic.index, 1]),
                 " with an AIC value of ",
                 aic.value, "."))
  }

  # Gather results
  best.model <- list("Estimates"     = estimates[[aic.index]],
                     "Prey Mat"      = new.prey.mat,
                     "Fat Vec"       = new.fat.vec,
                     "Species"       = results.table[aic.index, 1],
                     "AIC Value"     = aic.value)
  results    <- list("Best Model"    = best.model,
                     "AIC Values"    = results.table,
                     "Estimates"     = estimates)

  return(results)
}

