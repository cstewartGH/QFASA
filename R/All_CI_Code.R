#' Confidence Intervals for Diet Proportions
#'
#' Returns simultaneous confidence intervals for the diet of each species in the
#' prey database.
#'
#' @export
#' @param predator.mat matrix containing fatty acid signatures of the predators
#'     with fatty acids summing to one.
#' @param prey.mat prey database. A data frame with first column a
#'     Species label and other columns fatty acid proportions summing to one..
#' @param p.mat matrix of previously computed predator diet estimates needed for
#'     confidence interval calculation.
#' @param cal.mat matrix or vector of calibration coefficients of
#'     predators. Each COLUMN corresponds to a different predator. Default is a
#'     vector of ones. The number of fatty acids should be the same as the number
#'     of predator and prey fatty acids.
#' @param dist.meas distance measure to use for estimation: 1=KL,
#'     2=AIT or 3=CS
#' @param FC vector of prey fat content, one for each individual in prey database.
#'     Note that this vector is
#'     passed to the \code{\link{gen.pseudo.seals}} which expects fat
#'     content values for individual prey samples while
#'     \code{\link{pseudo.seal}} and  \code{\link{p.QFASA}}
#'     expect a species average.  Default is a vector of ones.
#' @param alpha 1-alpha is the family-wise or overall confidence level.
#'     Default is 0.05 for an overall confidence level of 0.95.
#' @param nprey number of prey to sample from the prey database
#'     when generating pseudo predators for the nuisance parameter
#'     estimation using original QFASA simulating code. Default is 30.
#' @param noise proportion of noise to include in the generation of pseudo predators
#'     using original QFASA simulating code.
#' @param R.p number of times to re-sample data.  Due to algorithm being slow,
#'     the default parameter is 1.
#' @param R.ps number of pseudo predators to generate when estimating
#'     nuisance parameters.  Default is 100.
#' @param R number of bootstrap replicates to use when generating
#'     p-values for confidence interval estimation.  Default is 100.
#' @param R.bias number of replicates for bias computation.  Default is 100.
#' @param ext.fa subset of fatty acids to be used.  These should be the same as
#'     those in predator.mat, prey.mat and cal.mat.
#' @return Simultaneous (1-alpha)*100% confidence intervals based on the
#'    zero-inflated beta distribution.
#' @details Intervals are biased corrected as recommended in Stewart, C. (2013).
#'    Intervals are slow to obtain, particularly if there are many prey types.
#'    See vignette on parallel execution to speed up calculations.
#' @references Stewart, C. (2013) Zero-inflated beta distribution for
#'     modeling the proportions in quantitative fatty acid signature
#'     analysis. Journal of Applied Statistics, 40(5), 985-992.
#'
#'
#'@examples
#'## Reducing prey database to three species so that code below will run more quickly.
#'## Please uncomment code to run.
#'
#'#set.seed(1234)
#'## Fatty Acids
#'#data(FAset)
#'#fa.set = as.vector(unlist(FAset))
#'
#'## Sample of Predators
#'#data(predatorFAs)
#'#predator.matrix = predatorFAs[, -c(1:4)]
#'#predator.matrix.ext = predatorFAs[,fa.set]
#'#predator.matrix.ext = predator.matrix.ext/rowSums(predator.matrix.ext)
#'
#'# Prey Database
#'#prey.red =
#'#preyFAs[preyFAs$Species=="capelin"|preyFAs$Species=="herring"|preyFAs$Species=="sandlance", ]
#'#prey.red = prey.red[,-c(1,3)]
#'#prey.red.ext = prey.red[,c("Species",fa.set)]
#'#prey.red.ext[,-1] <- prey.red.ext[,-1]/rowSums(prey.red.ext[,-1])
#'#prey.red.ext.means = MEANmeth(prey.red.ext)
#'
#'## Calibration Coefficients
#'
#'#data(CC)
#'#cal.vec = CC[CC$FA %in% fa.set, 2]
#'
#'#diet.est <- p.QFASA(predator.mat = predator.matrix.ext,
#'#                    prey.mat = prey.red.ext.means,
#'#                    cal.mat = cal.vec,
#'#                   dist.meas = 2,
#'#                    start.val = rep(1,nrow(prey.red.ext.means)),
#'#                    ext.fa = fa.set)[['Diet Estimates']]
#'
#'## conf.meth needs the full prey matrix unlike in p.QFASA
#'#ci <- conf.meth(predator.mat = predator.matrix.ext, prey.mat = prey.red.ext, cal.mat = cal.vec,
#'#          p.mat = diet.est, dist.meas = 2, ext.fa = fa.set)

conf.meth <- function(predator.mat,prey.mat,p.mat,cal.mat = rep(1, length(ext.fa)), dist.meas, FC = rep(1, nrow(prey.mat)),
                      alpha = 0.05,
                      nprey = 30,
                      R.p=1,
                      R.ps=100,
                      R=100,
                      R.bias=100,
                      noise=0,
                      ext.fa) {


  ## Number of prey species
  I <- length(unique(prey.mat[, 1]))

  ## Number of predator samples
  ns <- nrow(predator.mat)

  ## JUST TO AVOID PROBLEMS IN OTHER PROGRAMS THAT ARE EXPECTING A
  ## LIST OF LISTS (see parmeth.2.CI) ???
  beta.list <- vector("list", 1)
  data.star.seals <- predator.mat

  ## Generate estimates of beta distribution nuisance parameters
  ## (phi, theta) for R.ps pseudo predators.
  ## These are subsequently used in the confidence interval calculation where we
  ## bootstrap multiple p-values to account for sources of variance.
  ##
  for(r.p in 1:R.p) {
    futile.logger::flog.debug("Predator %d", r.p)

    ## Resample calibration factors
    #futile.logger::flog.debug("Resample calibration factors")
    if ( (is.vector(cal.mat)) || (nrow(cal.mat) == 1) || (ncol(cal.mat) == 1)) {
      #futile.logger::flog.debug("Only one calibration vector")
      data.star.cal <- matrix(cal.mat, ncol = 1)
    }
    else {
      cal.seq <- seq(1, ncol(cal.mat), 1)
      cal.sample <- sample(cal.seq, size = ncol(cal.mat), replace = T)
      data.star.cal <- cal.mat[, cal.sample]
    }
    ## Average resampled calibration coeffiecients
    data.star.cal.mean <- apply(data.star.cal, 1, mean)


    ## Resample prey and fat content
    #futile.logger::flog.debug("Resample prey and fat content")
    prey.seq <- seq(1, nrow(prey.mat), 1)
    prey.sample <- tapply(prey.seq, prey.mat[, 1], sample, replace = T)
    data.star.FC <- FC[unlist(prey.sample)]
    data.star.prey <- prey.mat[unlist(prey.sample),  ]
    data.star.prey.ext <- as.data.frame(data.star.prey)[c(dimnames(data.star.prey)[[2]][1], ext.fa)]
    data.star.prey.ext[, -1] <- data.star.prey.ext[, -1]/apply(data.star.prey.ext[, -1], 1, sum)

    ## Estimate diets of predators with FA signatures in predator.mat
    ## using the resampled calibration coefficients, prey
    ## signatures, and prey fat content
    p.mat.star <- as.matrix(p.QFASA(data.star.seals,
                                    MEANmeth(data.star.prey.ext),
                                    data.star.cal.mean,
                                    dist.meas,gamma=1,
                                    tapply(data.star.FC, prey.mat[, 1], mean, na.rm = T),
                                    ext.fa = ext.fa)[['Diet Estimates']])

    ## Generate R.ps pseudo-predators by sampling from the diets
    ## estimated above in p.mat.star. Then estimate the diet of
    ## these pseudo-predators, p.boot.mat, and use these estimates
    ## to bootstrap estimate nuisance parameters.
    #futile.logger::flog.debug("Estimate groups")
    R.ps.mod <- floor(R.ps/ns) * ns
    seq.vec <- rep(seq(1, ns, 1), floor(R.ps/ns))
    p.mat.initboot <- p.mat.star[seq.vec,  ]
    data.initboot <- comp.gen.pseudo.seals(data.star.prey,
                                           p.mat.initboot,
                                           data.star.cal,
                                           data.star.FC,
                                           0,
                                           nprey)
    seals.initboot <- data.initboot[[1]]
    cal.initboot <- data.initboot[[2]]
    FC.initboot <- data.initboot[[3]]
    p.boot.mat <- as.matrix(p.QFASA(seals.initboot,
                                    MEANmeth(data.star.prey.ext),
                                    cal.initboot,
                                    dist.meas,gamma=1,
                                    FC.initboot,
                                    ext.fa = ext.fa)[['Diet Estimates']])

    ## For each prey species, fit a beta distribution to the diets
    ## estimated on R.ps bootstrapped pseudo-predators generated above and
    ## estimate parameters theta and phi, the nuisance parameters.
    #futile.logger::flog.debug("Fit beta distributions")
    theta.beta <- rep(0, I)
    phi.beta <- rep(0, I)
    for (k in 1:I) {
      futile.logger::flog.debug("Prey species %d", k)
      if (sum(p.boot.mat[,k]!=0) >= 1) {

        est.gam <- gamlss::gamlss(p.boot.mat[, k] ~ 1,
                                  sigma.formula = ~ 1,
                                  nu.formula = ~ 1,
                                  control = (gamlss::gamlss.control(n.cyc=200,trace=F)),
                                  family = gamlss.dist::BEZI)

        phi.beta[k] <- exp(est.gam$sigma.coefficients)
        theta.beta[k] <- exp(est.gam$nu.coefficients)/(1+exp(est.gam$nu.coefficients))
      } else {
        theta.beta[k] = 1
      }
    }
    ## loop prey species

    ## Build nuisance parameter matrix: R.p x I x (phi, theta)
    if(r.p == 1) {
      par.list.prey <- list(rbind(theta.beta, phi.beta))
    }
    else {
      par.list.prey <- append(par.list.prey, list(rbind(theta.beta, phi.beta)))
    }
  }
  ## loop predators
  ##return(0)

  ## Confidence intervals
  ##
  ## Inputs:
  ## * beta.list: nuisance parameter matrix: R.p x I x (phi, theta).
  ## * R: number of bootstrap replicates to use in estimating p-values
  ## * alpha: confidence level for intervals.
  ## * p.mat: predator samples
  ###
  #futile.logger::flog.info("Estimate confidence intervals")
  beta.list[[1]] <- par.list.prey
  CI.L.1 <- rep(NA, I)
  CI.U.1 <- rep(NA, I)
  alpha1 <- alpha/I

  ## For each prey species
  for (k in 1:I) {
    futile.logger::flog.info("Computing confidence interval for prey species %d", k)
    CI <- opt.beta.lim(alpha1,
                       beta.list,
                       R,
                       p.mat,
                       k)
    #futile.logger::flog.info("Simultaneous %.4f confidence interval: [%.4f, %.4f]", alpha1, CI[[1]], CI[[2]])
    #futile.logger::flog.info("Individual %.4f confidence interval: [%.4f, %.4f]", alpha2, CI[[3]], CI[[4]])
    CI.L.1[k] <- CI[[1]]
    CI.U.1[k] <- CI[[2]]
  }

  out.mat <- rbind(CI.L.1, CI.U.1)
  colnames(out.mat) <- unique(prey.mat[,1])
  rownames(out.mat) <- c("Lower Limit","Upper Limit")

  est.bias <- bias.comp(p.mat,prey.mat,R.bias,dist.meas,ext.fa = ext.fa)

  out.mat[1,] <- out.mat[1,] - est.bias
  out.mat[2,] <- out.mat[2,] - est.bias

  out.mat[1,][out.mat[1,]<0] <- 0
  out.mat[2,][out.mat[2,]<0] <- 0

  out.mat[1,][out.mat[1,]>1] <- 1
  out.mat[2,][out.mat[2,]>1] <- 1

  # Return  simultaneous confidence intervals for each prey species
  return(out.mat)
}

#' Generate pseudo predators with ith predator having true diet
#' given by ith row of diet.null.mat.
#'
#' @param prey.mat matrix containing a representative FA signature
#'     from each prey group (usually the mean). The first column must
#'     index the prey group.
#' @param diet.null.mat null diet
#' @param cal.mat matrix of calibration coefficients where the \emph{i} th
#'     column is to be used with the \emph{i} th predator.
#' @param fat.cont vector of fat content
#' @param noise proportion of noise to add to the simulation.
#' @param nprey number of prey to sample.
#' @keywords internal
#'
comp.gen.pseudo.seals <-
  function(prey.mat, diet.null.mat, cal.mat, fat.cont, noise, nprey) {

    #futile.logger::flog.info("gen.pseudo.seals()")
    #futile.logger::flog.debug("Generating %d pseudo pradators.", nrow(diet.null.mat))
    ns <- nrow(diet.null.mat)
    I <- length(unique(prey.mat[, 1.]))
    nFA <- ncol(prey.mat) - 1.
    fat.sim <- rep(0., I)
    fat.mod <- matrix(rep(0., I * ns), byrow = T, ns, I)

    if ( !( (is.vector(cal.mat)) || (nrow(cal.mat) == 1.) || (ncol(cal.mat) == 1.)) )
    {
      ## Multiple calibration vectors provided: split calibration vectors into
      ## modeling and simulation set similarly to the prey database.
      ind.sim <- sample(seq(1., ncol(cal.mat), 1.), ns, replace = T)
      cal.sim <- as.matrix(cal.mat[, ind.sim])
      seals.pseudo <- matrix(rep(0., ns * nFA), byrow = T, ns, nFA)
      cal.mod <- matrix(rep(0., ns * nFA), byrow = T, nFA, ns)

      for(i in 1.:ns) {
        ## Split prey
        prey.out <- split_prey(prey.mat)
        prey.sim <- prey.out[[1]]
        prey.mod <- prey.out[[2]]

        ## Split fat content for each species
        for(k in 1.:I) {
          fat.split <- split_fatcont(fat.cont[prey.mat[, 1.] == unique(prey.mat[, 1.])[k]])
          fat.sim[k] <- fat.split[1]
          fat.mod[i, k] <- fat.split[2]
        }

        seals.pseudo[i,  ] <- pseudo.seal(prey.sim, diet.null.mat[i,  ], noise, nprey, cal.sim[, i], fat.sim, rep(F, I)) # is F supposed to be FALSE ???
        cal.mod[, i] <- apply(cal.mat[,  - ind.sim[i]], 1.,mean)
      }

      seals.pseudo <- as.data.frame(seals.pseudo)
      dimnames(seals.pseudo)[[2.]] <- dimnames(prey.mat[, -1.])[[2.]]
      cal.mod <- as.data.frame(cal.mod)
      dimnames(cal.mod)[[1.]] <- dimnames(prey.mat[, -1.])[[2.]]
      seals.pseudo <- as.matrix(seals.pseudo)
      cal.mod <- as.matrix(cal.mod)

    } else {

      ## Only a single calibration vector supplied so we cannot split.
      #futile.logger::flog.debug("WARNING IN gen.pseudo.seals:  ONE CALIBATION")
      seals.pseudo <- matrix(rep(0., ns * nFA), byrow = T, ns,nFA)

      for(i in 1.:ns) {
        ## Split prey
        prey.out <- split_prey(prey.mat)
        prey.sim <- prey.out[[1]]
        prey.mod <- prey.out[[2]]

        ## Split fat content for each species
        for(k in 1.:I) {
          fat.split <- split_fatcont(fat.cont[prey.mat[  , 1.] == unique(prey.mat[,1.])[k]])
          fat.sim[k] <- fat.split[1]
          fat.mod[i, k] <- fat.split[2]
        }

        seals.pseudo[i,  ] <- pseudo.seal(prey.sim, diet.null.mat[i,  ], noise, nprey, cal.mat, fat.sim, rep(F, I))
      }
      cal.mod <- cal.mat
      seals.pseudo <- as.data.frame(seals.pseudo)
      dimnames(seals.pseudo)[[2.]] <- dimnames(prey.mat[, -1.])[[2.]]
      seals.pseudo <- as.matrix(seals.pseudo)
    }

    return(list(seals.pseudo, cal.mod, fat.mod))
  }

#' Find simultaneous confidence intervals for diet
#' proportions of a single prey species i.e. solve f(pio) = PVAL(pio)
#' = alpha1.  Calls root finding function root.beta.
#'
#' @param alpha1 simultaneous confidence level
#' @param par.list a list of R.p lists of I beta distribution parameters phi
#'     and theta that define diet proportion estimates for each of the
#'     prey species. Effectively R.p beta distributions for each of the
#'     I prey species from which we bootstrap to calculate p-values.
#' @param R number of bootstrap replicates to use in p-value estimation.
#' @param p.mat predator diet estimates.
#' @param k prey species index 1...I
#' @return upper and lower limits simultaneous
#'     intervals respectively.
#' @keywords internal
#'
opt.beta.lim <- function(alpha1, par.list, R, p.mat, k) {

  #futile.logger::flog.debug("opt.beta.lim()")

  #cat("Finding lower limit")
  #futile.logger::flog.debug("Finding lower limit (%f)", 1-alpha1)
  CI.L = root.beta(0,mean(p.mat[,k]),alpha1,par.list,R,p.mat,k)
  #print(CI.L)

  #cat("Finding upper limit")
  #futile.logger::flog.debug("Finding U2 (%f)", 1-alpha1)
  CI.U = root.beta(mean(p.mat[,k]),1, alpha = alpha1,par.list,R,p.mat,k,low=FALSE)
  #print(CI.U)

  return(list(CI.L,CI.U))
}

#' Find root (i.e. solve f(pio) = PVAL(pio) = alpha1 using either uniroot or optimize.
#'
#' @param x1 lower starting value
#' @param x2 upper starting value
#' @param alpha simultaneous confidence level
#' @param par.list parameters needed to compute p-value by bootstrapping
#' @param R number of bootstraps
#' @param p.mat predator diet estimates
#' @param k prey species index 1...I
#' @param low TRUE if computing lower interval.
#' @return root which will be lower or upper CI limit
#' @keywords internal
root.beta <- function(x1, x2, alpha, par.list, R, p.mat, k,low=TRUE) {

  # cat("in root beta")


  f.lower = comp.beta.pval.k(par.list, R, x1, p.mat, k) - alpha
  f.upper = comp.beta.pval.k(par.list, R, x2, p.mat, k) - alpha


  if (f.lower*f.upper < 0 ) {

    r = stats::uniroot(function(x) { comp.beta.pval.k(par.list, R, x, p.mat, k) - alpha },
                       c(x1, x2), f.lower = f.lower, f.upper = f.upper,trace = 0)
    return(r$root)

  } else {

    if (low==TRUE && round(x2,3)==0) {

      #cat("Finding lower limit and x2 approx. equal to 0")

      return(0)

    } else if(low==FALSE && round(x1,3)==1) {

      #cat("Finding upper limit and x1 aprox. equal to 1")

      return(1)

    } else if (low==FALSE && round(x1,3)==0) {

      #cat("Finding upper limit x1 approx. equal to 0 so pass smaller interval to start.")

      r1 = stats::optimize(function(x) { abs(comp.beta.pval.k(par.list, R, x, p.mat, k) - alpha) },
                           c(x1, x1+0.05),tol=0.001)
      r2 = stats::optimize(function(x) { abs(comp.beta.pval.k(par.list, R, x, p.mat, k) - alpha) },
                           c(x1+0.05, x1+0.15),tol=0.001)
      r3 = stats::optimize(function(x) { abs(comp.beta.pval.k(par.list, R, x, p.mat, k) - alpha) },
                           c(x1+0.15, x1+0.35),tol=0.001)

      order.obj.ind <- order(r1$objective,r2$objective,r3$objective)[1]

      return(c(r1$minimum,r2$minimum,r3$minimum)[order.obj.ind])

    } else if (low==TRUE && round(x2,3)==1) {

      #cat("Finding lower limit and x2 aprox. equal to 1 so pass smaller interval to start.")

      r1 = stats::optimize(function(x) { abs(comp.beta.pval.k(par.list, R, x, p.mat, k) - alpha) },
                    c(x2-0.05, x2),tol=0.001)
      r2 = stats::optimize(function(x) { abs(comp.beta.pval.k(par.list, R, x, p.mat, k) - alpha) },
                   c(x2-0.15, x2-0.05),tol=0.001)
      r3 = stats::optimize(function(x) { abs(comp.beta.pval.k(par.list, R, x, p.mat, k) - alpha) },
                   c(x2-0.35, x2-0.15),tol=0.001)

      order.obj.ind <- order(r1$objective,r2$objective,r3$objective)[1]

      return(c(r1$minimum,r2$minimum,r3$minimum)[order.obj.ind])


    } else {

      cat("***Could not find limit***")

      return(NA)

    }

  }


}

#' Calculate p-value corresponding to a specified null value using bootstrapping
#'
#' @param par.list a list of R.p lists of I beta distribution parameters phi
#'     and theta that define diet proportion estimates for each of the
#'     prey species. Effectively R.p beta distributions for each of the
#'     I prey species which we bootstrap to calculate p-values.
#' @param R number of bootstrap replicates to use in p-value estimation.
#' @param diet.test.k null value
#' @param p.mat predator diet estimates.
#' @param k prey species index 1..I
#' @return p-value p-value
#' @keywords internal
comp.beta.pval.k <- function(par.list, R, diet.test.k, p.mat, k) {

  I <- ncol(p.mat)
  ns <- nrow(p.mat)
  T.obs <- abs(p.beta(p.mat[,k]) - diet.test.k)
  R.s <- length(par.list)
  R.p <- length(par.list[[1.]])
  pval.vec <- rep(0., R.s)

  for(r.s in 1.:R.s) {
    pval.vec.prey <- rep(0., R.p)

    for(r.p in 1.:R.p) {
      ## Thetas
      prob.zero.vec <- par.list[[r.s]][[r.p]][1.,  ]

      ## Phis
      phi.est <- par.list[[r.s]][[r.p]][2.,  ]

      if(round(diet.test.k, 6.) == 0.) {
        if(round(prob.zero.vec[k], 6.) == 1.) {
          pval.vec.prey[r.p] <- stats::runif(1.)
        } else {
          pval.vec.prey[r.p] <- 0.
        }
      }
      else if(round(prob.zero.vec[k], 6.) == 1.) {
        pval.vec.prey[r.p] <- 0.
      }
      else if(round(diet.test.k, 6.) >= (round(1. - prob.zero.vec[k], 6.))) {
        pval.vec.prey[r.p] <- 0.
      }
      else {
        ## Mean of ZIB is (1-theta).mu where mu is the mean of
        ## the underlying Beta distribution
        mu.null <- diet.test.k/(1-prob.zero.vec[k])


        T.star <- comp.Tstar.beta(p.mat[, k],
                                  prob.zero.vec[k],
                                  mu.null,
                                  phi.est[k],
                                  R)

        futile.logger::flog.debug("    T.star = %f", T.star)
        T.star <- abs(T.star - diet.test.k)

        ## p-value
        pval.vec.prey[r.p] <- sum(T.star >= T.obs)/R
      }
    }
    ## Average p-values over R.p ?
    pval.vec[r.s] <- mean(pval.vec.prey)
  }

  ## Average p-values over R.s ?
  pval <- mean(pval.vec)
  return(pval)
}

#' Generate bootstrap replicates of diet proportion estimates for
#' given prey species
#'
#' @param p.k list of predator diet proportions for prey species k
#' @param theta.boot distribution parameter
#' @param mu.null distribution parameter
#' @param phi.boot distribution parameter
#' @param R.boot number of bootstrap replicates to create
#' @return list of diet proportion replicates
#' @keywords internal
comp.Tstar.beta <- function(p.k, theta.boot, mu.null, phi.boot, R.boot) {
  futile.logger::flog.trace("Tstar.beta()")
  data.para <- boot::boot(data = p.k,
                          statistic = comp.p.beta,
                          R = R.boot,
                          sim = "parametric",
                          ran.gen = data.sim.beta,
                          mle = list(length(p.k), mu.null, phi.boot, theta.boot),
                          parallel = 'snow',
                          ncpus = getOption('qfasa.parallel.numcores', 1),
                          cl =  getOption('qfasa.parallel.cluster', NULL))

  return(data.para$t)
}

#' Bootstrap statistic function: in this case it is the mean
#' (meanBEZI) of the fitted (gamlss) ZIB distribution.
#'
#' @param data data
#' @return the expected value of the response for a fitted model
#' @keywords internal
comp.p.beta <- function(data){

  # FOUND ERRORS IN HOLLY'S CODE SO RE-WROTE

  p.vec <- data
  ns <- length(p.vec)

  no.zeros <- sum(round(p.vec,3)==0)

  if (no.zeros==ns)
  { return(0)
  }


  est.gam <- gamlss::gamlss(p.vec~1, sigma.formula=~1,nu.formula=~1,control=(gamlss::gamlss.control(n.cyc=1000,trace=F)),
                            family = gamlss.dist::BEZI)

  return(gamlss.dist::meanBEZI(est.gam)[1])
}


#' Bootstrap ran.gen function:
#' @keywords internal
#'
data.sim.beta <- function(data, mle) {
  if (mle[[4]]==0) { mle[[4]] <- 0.00001 }
  d <- gamlss.dist::rBEZI(mle[[1]], mle[[2]], mle[[3]], mle[[4]])
  d[d>0.999] <- 0.999
  return(d)
}

#' Calculate bias correction for confidence intervals.
#'
#' @param diet.est predator diet estimates
#' @param prey.mat matrix of prey FA signatures
#' @param R.bias number of replicates
#' @param dist.meas distance measure
#' @param ext.fa subset of FA's to use.
#' @return estimate of bias for each prey type
#' @keywords internal
bias.comp <- function(diet.est,prey.mat,R.bias,dist.meas, ext.fa, gamma = 1) {
  # Written February 3rd, 2021 to replace bias.all
  # Simulating and modeling with UNIT calibration and fat content vectors

  cal.vec <- rep(1,ncol(prey.mat)-1)
  fat.vec <- rep(1,ncol(diet.est))
  p.mat <- matrix(rep(NA,R.bias*ncol(diet.est)),nrow=R.bias)
  beta.est <- matrix(0, nrow(diet.est), ncol(diet.est))
  bias.mat <- matrix(0, nrow(diet.est), ncol(diet.est))

  for (i in 1:nrow(diet.est)) {
    #cat("i")
    #print(i)
    for (n in 1: R.bias) {
      pred.sig <- pseudo.pred(diet.est[i,],prey.mat,cal.vec,fat.vec)
      p.mat[n,] <- p.QFASA(pred.sig, MEANmeth(prey.mat),cal.vec,dist.meas,ext.fa=ext.fa)$`Diet Estimates`
    }

    bias.mat[i,] <- colMeans(p.mat) - diet.est[i,]

  }  # end i

  bias.est <- apply(bias.mat,2,stats::median)

  return(bias.est)

}

#' Splits prey database into a simulation set (1/3) and a modelling set (2/3).
#' Returns a list:
#'
#' 1. simulation prey database
#' 2. modelling prey database
#'
#' IF number of samples of a prey type <=5, then prey.mod AND prey.sim
#' are duplicated instead of split.
#'
#' @param prey.mat matrix of individual prey fatty acid signatures
#'     where the first column denotes the prey type
#'
split_prey <- function(prey.mat) {

  numprey <- tapply(prey.mat[, 1.], prey.mat[, 1.], length)

  I <- length(numprey)
  n.sim <- floor(numprey/3.)
  n.mod <- numprey - n.sim
  split.list <- tapply(seq(1, nrow(prey.mat), by=1), prey.mat[, 1.],sample)
  prey.mat.rand <- prey.mat[unlist(split.list),  ]

  if (numprey[1.] > 5.) {
    split.bool.sim <- c(rep(TRUE, n.sim[1.]), rep(FALSE, n.mod[1.]))
    split.bool.mod <- c(rep(FALSE, n.sim[1.]), rep(TRUE, n.mod[1.]))
  } else
  {
    split.bool.sim <- c(rep(TRUE, numprey[1.]))
    split.bool.mod <- c(rep(TRUE, numprey[1.]))
  }
  for (i in 2.:I) {
    if (numprey[i] > 5.) {
      split.bool.sim <- c(split.bool.sim, c(rep(TRUE, n.sim[i]), rep(FALSE, n.mod[i])))
      split.bool.mod <- c(split.bool.mod, c(rep(FALSE, n.sim[i]), rep(TRUE, n.mod[i])))
    }
    else      {
      split.bool.sim <- c(split.bool.sim, rep(TRUE, numprey[i]))
      split.bool.mod <- c(split.bool.mod, rep(TRUE, numprey[i]))
    }
  }

  prey.sim <- prey.mat.rand[split.bool.sim,  ]
  prey.mod <- prey.mat.rand[split.bool.mod,  ]

  split.prey.list <- list(prey.sim, prey.mod)
  return(split.prey.list)
}

split_fatcont <- function(fat.cont.k) {

  ## USED TO RANDOMLY SPLIT FAT CONTENT FOR SPECIES K IN TWO
  fat.cont.k <- fat.cont.k[!is.na(fat.cont.k)]
  fat.sample <- sample(fat.cont.k)
  fat.1 <- fat.sample[1.:round(length(fat.cont.k)/2.)]
  fat.1.mean <- mean(fat.1)
  fat.2 <- fat.sample[(round(length(fat.cont.k)/2.) + 1.):length(fat.cont.k)]
  fat.2.mean <- mean(fat.2)
  return(c(fat.1.mean, fat.2.mean))
}

##' Generate a single pseudo predator FA signature
##'
##' THIS IS THE NEW pseudo.seal FUNCTION THAT ALLOWS 1) FAT CONTENT
##' TO BE INCLUDED IN THE GENERATED SEALS AND 2) SOME SPECIES TO BE
##' TRULY ZERO (THAT IS,"ZERO SPECIES" DO NOT HAVE TO BE INCLUDED IN THE
##' "NOISE" )
##' NOTE:  IT IS ASSUMED THAT SUM(DIET) IS 1-NOISE
##'
##' @param prey.sim OUTPUT OF split.prey
##' @param diet DIET COMPOSITION VECTOR  (NOTE: THIS VECTOR SHOULD SUM TO 1-NOISE.  THE NOISE WILL BE ADDED TO THE diet VECTOR.)
##' @param noise AMOUNT OF NOISE
##' @param nprey nprey TOTAL NUMBER OF PREY TO BE SAMPLED
##' @param cal CALIBRATION FACTORS
##' @param fat.cont VECTOR OF FAT CONTENT OF LENGTH=I (# OF SPECIES)
##' @param specify.noise A BOOLEAN VECTOR WITH TRUES DENOTING SPECIES TO USE IN NOISE.
##' @keywords internal
##' @return seal.star SIMULATED SEAL FA SIGNATURE.
pseudo.seal <- function(prey.sim, diet, noise, nprey, cal, fat.cont, specify.noise) {

  ## MODIFYING DIET TO INCLUDE FAT CONTENT
  diet <- fat.cont * diet

  ## WANT sum(diet) + noise = 1
  ## MODIFYING NOISE TO INCLUDE FAT CONTENT

  if (noise != 0) {
    fat.cont.mean.noise <- mean(fat.cont[specify.noise])
    noise <- fat.cont.mean.noise * noise
  }

  diet.mod <- diet/(sum(diet) + noise)
  noise <- noise/(sum(diet) + noise)
  diet <- diet.mod
  numprey.diet <- round(nprey * diet)
  numprey.noise <- round(nprey * noise)

  ## FIRST SELECT FROM prey.sim, THE PREY INCLUDED IN THE DIET
  I <- length(unique(prey.sim[, 1.]))
  numprey.sim <- tapply(prey.sim[, 1.], prey.sim[, 1.], length)


  ## SAMPLE WITH REPLACEMENT numprey.diet FROM prey.sim.in
  prey.sim.in <- prey.sim[rep(numprey.diet != 0., numprey.sim),  ]
  numprey.diet.bool <- numprey.diet != 0.
  nonzero.ind <- seq(1., I, 1.)[numprey.diet.bool]

  ## HERE
  ind.list <- list(sample(seq(1., numprey.sim[nonzero.ind[1.]], 1.), as.numeric(numprey.diet[nonzero.ind[1.]]), replace = T))
  if (sum(numprey.diet != 0) > 1.) {
    for (i in 2.:sum(numprey.diet != 0.)) {
      ind.list <- append(ind.list, list(sample(seq(1., numprey.sim[nonzero.ind[i]], i),
                                               as.numeric(numprey.diet[nonzero.ind[i]]), replace = T)))
    }
  }

  numprey.sim.in <- numprey.sim[numprey.diet.bool]
  ind <- unlist(ind.list)
  ind <- ind + rep(cumsum(numprey.sim.in), numprey.diet[numprey.diet.bool]) -
    rep(numprey.sim.in, numprey.diet[numprey.diet.bool])

  ## SELECT FROM prey.sim, THE PREY INCLUDED IN NOISE
  prey.star <- prey.sim.in[ind,  ]
  if (noise != 0.) {
    prey.sim.out <- prey.sim[rep(specify.noise, numprey.sim),  ]
  }

  if ((noise != 0.) & (numprey.noise != 0.)) {
    sample.repl <- sample(seq(1., nrow(prey.sim.out), 1.), numprey.noise, replace
                          = T)
    if (length(sample.repl) == 1.) {
      noise.star <- prey.sim.out[sample.repl,  ]
      seal.star <- (apply(prey.star[, -1.], 2., sum) +  noise.star[1., -1.])/(nprey) * cal
    }
    else {
      noise.star <- prey.sim.out[sample.repl,  ]
      seal.star <- (apply(prey.star[, -1.], 2., sum) +
                      apply(noise.star[, -1.], 2., sum))/(nprey) * cal
    }
  }
  else {
    seal.star <- apply(prey.star[, -1.], 2., sum)/(nprey) * cal
  }

  seal.star <- seal.star/sum(seal.star)
  seal.star <- as.vector(seal.star)
  names(seal.star) <- dimnames(prey.sim[, -1.])[[2.]]
  return(seal.star)
}


