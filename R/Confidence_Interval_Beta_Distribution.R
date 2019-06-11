#'
#' Returns individual confidence intervals and simultaneous confidence intervals
#' based on the zero-inflated beta distribution (not bias corrected - see note below).
#' 
#' For details see:  
#'     Stewart, C. (2013)  Zero-Inflated Beta Distribution for Modeling the Proportions in 
#'     Quantitative Fatty Acid Signature Analysis.  
#'     Journal of Applied Statistics, 40(5), 985-992.
#' 
#' Note:
#' \itemize{
#'     \item These intervals are biased and should be corrected using the
#'           output from \code{\link{bias.all()}}.
#'     \item \code{CI.L.1} and \code{CI.U.1} contain the simultaneous
#'           confidence intervals.
#'     \item Slow because of bisection and lots of repetition.
#'     \item Need to replace \code{p.prey} with
#'           \code{\link{p.QFASA()}} eventually but just use
#'           p.prey() for now. Use example where we estimate a single
#'           predator diet to compare the estimates from each method.
#' }
#' 
#' @export
#' @param predator.mat matrix containing the fatty acid signatures of the predators. 
#' @param prey.mat prey database. A dataframe with first column a
#'     Species label and other columns fatty acid proportions. Fatty
#'     acid proportions are compositional. 
#' @param cal.mat matrix of calibration coefficients of
#'     predators. Each column corresponds to a different predator. At
#'     least one calibration coefficient vector must be supplied.
#' @param dist.meas distance measure to use for estimation: 1=KL,
#'     2=AIT or 3=CS
#' @param noise proportion of noise to include in the simulation.
#' @param nprey number of prey to sample from the the prey database
#'     when generating pseudo-predators for the nuisance parameter
#'     estimation. 
#' @param R.p number of beta diet distributions to generate for the
#'     nuisance parameters.
#' @param R.ps number of pseudo predators to generate when estimating
#'     nuisance parameters.
#' @param R number of bootstrap replicates to use when generating
#'     p-values for confidence interval estimation.
#' @param p.mat matrix of predator diet estimates for which we are
#'     trying to find confidence interavls. 
#' @param alpha confidence interval confidence level.
#' @param FC vector of prey fat content. Note that this vector is
#'     passed to the \code{\link{gen.pseudo.seals()}} which expects fat
#'     content values for individual prey samples while
#'     \code{\link{pseaudo.seal()}} and  \code{\link{p.QFASA()}}
#'     expect a species average.
#' @param ext.fa subset of fatty acids to be used to obtain QFASA diet
#'     estimates.
#' @return Individual confidence intervals and simultaneous confidence
#'     intervals based on the zero-inflated beta distribution. These
#'     intervals are biased and should be corrected using the output
#'     from \code{\link{bias.all()}}. \code{ci.l.1} and \code{ci.u.1}
#'     contain the simultaneous confidence intervals. 
#' @references Stewart, C. (2013) Zero-inflated beta distribution for
#'     modeling the proportions in quantitative fatty acid signature
#'     analysis. Journal of Applied Statistics, 40(5), 985-992. 
#'
#'
#' @examples
#'  ## Fatty Acids
#'  data(FAset)
#'  fa.set = as.vector(unlist(FAset))
#'  
#'  ## Predators
#'  data(predatorFAs)
#'  tombstone.info = predatorFAs[,1:4]
#'  predator.matrix = predatorFAs[,5:(ncol(predatorFAs))]
#'  npredators = nrow(predator.matrix)
#'
#'  ## Prey
#'  data(preyFAs)
#'  prey.sub=(preyFAs[,4:(ncol(preyFAs))])[fa.set]
#'  prey.sub=prey.sub/apply(prey.sub,1,sum) 
#'  group=as.vector(preyFAs$Species)
#'  prey.matrix=cbind(group,prey.sub)
#'  prey.matrix=MEANmeth(prey.matrix) 
#'
#' 
#' # Diet estimate
#' diet.est <- p.QFASA(predator.mat = predator.matrix,
#'                     prey.mat = prey.db.summarized,.matrix,
#'                     cal.mat = rep(1, nrow(FAset)),
#'                     dist.meas = 2,
#'                     ext.fa = colnames(prey.db.summarized))[['Diet Estimates']]
#' 
#' # Confidence intervals
#' ci = beta.meths.CI(predator.mat = predator.matrix,
#'                    prey.mat = prey.matrix,
#'                    cal.mat = rep(1, nrow(FAset)),
#'                    dist.meas = 2,
#'                    noise = 0,
#'                    nprey = 50,
#'                    R.p = 1,
#'                    R.ps = 10,
#'                    R = 10, 
#'                    p.mat = diet.est,
#'                    alpha = 0.05,
#'                    FC = rep(1, nrow(prey.db)),
#'                    ext.fa = FAset$FA)
#'
#' # Bias correction
#' bias <- bias.all(p.mat = diet.est,
#'                  prey.mat = prey.matrix,
#'                  cal.mat = as.matrix(rep(1, nrow(FAset))),
#'                  fat.cont = rep(1, nrow(prey.db)),
#'                  R.bias = 10,
#'                  noise = 0,
#'                  nprey = 50,
#'                  specify.noise = rep(FALSE, nspecies),
#'                  dist.meas = 2,
#'                  ext.fa = FAset$FA)
#' 
#' # SIMULTANEOUS CONFIDENCE INTERVALS:
#' # LOWER LIMIT
#' ci[[1]] - bias[3,]
#' 
#' # UPPER LIMIT 
#' ci[[2]] - bias[3,]
#'
beta.meths.CI <- function(predator.mat,
                          prey.mat,
                          cal.mat,
                          dist.meas,
                          noise,
                          nprey,
                          R.p,
                          R.ps,
                          R,
                          p.mat,
                          alpha,
                          FC,
                          ext.fa)
{
    futile.logger::flog.info("beta.meths.CI()")

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
        futile.logger::flog.debug("Resample calibration factors")
        if ( (is.vector(cal.mat)) || (nrow(cal.mat) == 1) || (ncol(cal.mat) == 1)) {
            futile.logger::flog.debug("Only one calibration vector")
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
        futile.logger::flog.debug("Resample prey and fat content")
        prey.seq <- seq(1, nrow(prey.mat), 1)
        prey.sample <- tapply(prey.seq, prey.mat[, 1], sample, replace = T)
        data.star.FC <- FC[unlist(prey.sample)]
        data.star.prey <- prey.mat[unlist(prey.sample),  ]
        data.star.prey.ext <- as.data.frame(data.star.prey)[c(dimnames(data.star.prey)[[2]][1], ext.fa)]
        data.star.prey.ext[, -1] <- data.star.prey.ext[, -1]/apply(data.star.prey.ext[, -1], 1, sum)
        
        ## Estimate diets of predators with FA signatures in predator.mat
        ## using the resampled calibration coefficients, prey
        ## signatures, and prey fat content
        p.mat.star <- as.matrix(p.prey(data.star.seals, 
                                       MEANmeth(data.star.prey.ext), 
                                       data.star.cal.mean, 
                                       dist.meas,
                                       tapply(data.star.FC, prey.mat[, 1], mean, na.rm = T),
                                       ext.fa = ext.fa))

        ## Generate R.ps pseudo-predators by sampling from the diets
        ## estimated above in p.mat.star. Then estimate the diet of
        ## these pseudo-predators, p.boot.mat, and use these estimates
        ## to bootstrap estimate nuisance parameters.
        futile.logger::flog.debug("Estimate groups")
        R.ps.mod <- floor(R.ps/ns) * ns
        seq.vec <- rep(seq(1, ns, 1), floor(R.ps/ns))
        p.mat.initboot <- p.mat.star[seq.vec,  ]
        data.initboot <- gen.pseudo.seals(data.star.prey,
                                          p.mat.initboot,
                                          data.star.cal,
                                          data.star.FC,
                                          0,
                                          nprey)
        seals.initboot <- data.initboot[[1]]
        cal.initboot <- data.initboot[[2]]
        FC.initboot <- data.initboot[[3]]
        p.boot.mat <- p.prey(seals.initboot,
                             MEANmeth(data.star.prey.ext),
                             cal.initboot,
                             dist.meas,
                             FC.initboot,
                             ext.fa = ext.fa)
        
        ## For each prey species, fit a beta distribution to the diets
        ## estimated on R.ps bootstrapped pseudo-predators generated above and
        ## estimate parameters theta and phi, the nuisance parameters.   
        futile.logger::flog.debug("Fit beta distributions")
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


    ## Confidence intervals
    ##
    ## Inputs:
    ## * beta.list: nuisance parameter matrix: R.p x I x (phi, theta).
    ## * R: number of bootstrap replicates to use in estimating p-values
    ## * alpha: confidence level for intervals.
    ## * p.mat: predator samples
    ###
    futile.logger::flog.info("Estimate confidence intervals")
    beta.list[[1]] <- par.list.prey
    CI.L.1 <- rep(NA, I)
    CI.U.1 <- rep(NA, I)
    CI.L.2 <- rep(NA, I)
    CI.U.2 <- rep(NA, I)
    alpha1 <- alpha/I
    alpha2 <- alpha
    futile.logger::flog.debug("Simultaneous confidence level: %.4f", alpha1)
    futile.logger::flog.debug("Individual confidence level: %.4f", alpha2)
   
    ## For each prey species
    for (k in 1:I) {
        futile.logger::flog.info("Prey species %d", k)
        CI <- bisect.beta.lim(alpha1,
                              alpha2,
                              beta.list,
                              R,
                              p.mat,
                              k)
        futile.logger::flog.info("Simultaneous %.4f confidence interval: [%.4f, %.4f]", alpha1, CI[[1]], CI[[2]])
        futile.logger::flog.info("Individual %.4f confidence interval: [%.4f, %.4f]", alpha2, CI[[3]], CI[[4]])
        CI.L.1[k] <- CI[[1]]
        CI.U.1[k] <- CI[[2]]
        CI.L.2[k] <- CI[[3]]
        CI.U.2[k] <- CI[[4]]

        ## Comparison with original method
        ##CI2 <- bisect.beta.lim.orig(alpha1,
        ##                            alpha2,
        ##                            beta.list,
        ##                            R,
        ##                            p.mat,
        ##                            k)
        ##futile.logger::flog.info("Orig Simultaneous %.4f confidence interval: [%.4f, %.4f]",
        ##           alpha1, CI2[[1]], CI2[[2]])
        ##futile.logger::flog.info("Orig Individual %.4f confidence interval: [%.4f, %.4f]",
        ##          alpha2, CI2[[3]], CI2[[4]])                  
    }

    # Return individual and simlutaneous confidence intervals for each prey species
    return(list(CI.L.1, CI.U.1, CI.L.2, CI.U.2))
}


MEANmeth <- function(prey.mat) {
    
    ## RETURNS THE MULTIVARIATE MEAN FA SIGNATURE FROM EACH PRETY TYPE
    
    ## INPUT:
    ## prey.mat --> MATRIX CONTAINING FA SIGNATURES OF THE PREY.  NOTE THAT
    ##                FIRST COLUMN INDEXES PREY TYPE.
    
    prey.means <- apply(prey.mat[, -1.], 2., tapply, prey.mat[, 1.], mean)
    return(prey.means)    
}


KLdist <- function(x, y)
{
    ## SIMILAR TO KL BUT REQUIRES TWO VECTORS AS INPUT
    ## RETURNS KULBACK-LIEBLER DISTANCE BETWEEN x AND y WHERE 
    ## INPUT:
    ## x --> VECTOR 
    ## y --> vector
    return(sum((x - y) * log(x/y)))
}


## Dont use this one...use QFASA::p.QFASA()
p.prey <- function(seal.mat, prey.mat, cal.mat, dist.meas, FC = rep(1., nrow(prey.mat)),
                   start.val = rep(1., nrow(prey.mat)), ext.fa ) {
  
  # COMPUTES THE DIET ESTIMATE FOR EACH SEAL IN seal.mat USING KL (1) OR
  # AIT DISTANCE MEASURE (2).
  
  # NOTE: seal.mat MUST BE MATRIX AND prey.mat CONTAINS A
  # REPRESENTATIVE FA SIGNATURE (SUCH AS THE MEAN) FROM
  # EACH SPECIES.
  
  # INPUT:
  # seal.mat   --> MATRIX CONTAINING FA SIGNATURES OF SEALS
  # prey.mat   --> MATRIX CONTAINING A REPRESENTATIVE FA SIGNATURE
  #                          FROM   EACH PREY TYPE.
  # *** ASSUMES THAT FIRST COLUMN CONTAINS NAME OF SPECIES ***
  # cal.mat     --> MATRIX OF CALIBRATION FACTORS WHERE THE iTH COLUMN IS TO
  #                      BE USED WITH THE itH SEAL.
  # dist.meas --> DISTANCE MEASURE TO USE: 1 = KL OR 2 = AITCHISON
  
  # OUTPUT:
  # p.mat  --> VECTOR THAT MINIMIZED THE FUNCTION (NORMALIZED)
  
  
  if ((is.vector(cal.mat)) || (nrow(cal.mat) == 1.) || (ncol(cal.mat ) == 1.)) {
    
    #cat("IN p.prey WARNING:  ONE SEAL")
    
    # ONLY ONE SEAL
    # REMOVE EXTENDED DIETARY
    
    seal.mat <- t(t(seal.mat)/as.vector(unlist(cal.mat)))
    seal.mat <- as.data.frame(seal.mat)[ext.fa]
    seal.mat <- seal.mat/apply(seal.mat, 1., sum)
    seal.mat <- as.matrix(seal.mat)
    
    
  }
  else {
    # REMOVE EXTENDED DIETARY
    
    seal.mat <- seal.mat/t(cal.mat)
    seal.mat <- as.data.frame(seal.mat)[ext.fa]
    seal.mat <- seal.mat/apply(seal.mat, 1., sum)
    seal.mat <- as.matrix(seal.mat)
  }
  
  # eps USED TO ADJUST DIET ESTIMATES THAT ARE 1.
  
  eps <- 0.0005
  I <- nrow(prey.mat)
  ns <- nrow(seal.mat)
  
  
  
  p.mat <- matrix(rep(0., nrow(prey.mat) * nrow(seal.mat)), byrow = T, 
                  nrow(seal.mat),
                  nrow(prey.mat))
  if (!(is.matrix(start.val))) {
    start.val <- matrix(rep(start.val, ns), byrow = T, ns, I)
  }
  
  if (dist.meas == 1.) {
    for (i in 1.:nrow(seal.mat)) {
      p.all <- nlminb(start = start.val[i,  ], objective = optquantile.obj2, 
                      gradient = optquantile.grad2,
                      seal = seal.mat[i,  ], prey.quantiles = 
                        prey.mat, lower =
                        rep(0.,nrow(prey.mat)), upper = rep(1., 
                                                            nrow(prey.mat)), control =
                        list(iter.max = 1000, eval.max = 1000.))
      
      p.mat[i,  ] <- p.all$par
      p.mat[i,  ] <- p.mat[i,  ]/sum(p.mat[i,  ])
      eq.one.bool <- p.mat[i,  ] == 1.
      if (any(eq.one.bool)) {
        p.mat[i,  ][eq.one.bool] <- 1. - eps
        p.mat[i,  ][!eq.one.bool] <- eps/I
      }
    } # END for
  }
  else {
    for (i in 1.:nrow(seal.mat)) {
      p.all <- nlminb(start = start.val[i,  ], objective = optcompdiff.obj,gradient = optcompdiff.grad, seal = 
                        seal.mat[i,  ],
                      prey.quantiles = prey.mat, lower = rep(0., nrow(prey.mat)),
                      upper = rep(1., nrow(prey.mat)), control = list(iter.max = 1000, eval.max = 1000.))
      p.mat[i,  ] <- p.all$par
      p.mat[i,  ] <- p.mat[i,  ]/sum(p.mat[i,  ])
      eq.one.bool <- p.mat[i,  ] == 1.
      if (any(eq.one.bool)) {
        p.mat[i,  ][eq.one.bool] <- 1. - eps
        p.mat[i,  ][!eq.one.bool] <- eps/I
      }
    } # END for
  }
  
  if (is.matrix(FC))
  { FC.mat <- FC
  } else { 
    FC.mat <- matrix(rep(FC, nrow(seal.mat)), byrow = T, nrow(seal.mat), I)
  } 
  
  p.mat <- p.mat/FC.mat
  p.mat <- p.mat/apply(p.mat, 1., sum)
  return(p.mat)
}
        
optquantile.obj2 <- function(alpha, seal, prey.quantiles){
  
  # USED IN optquantile (IN CALL TO nlminb) AS THE OBJECTIVE FUNCTION TO BE
  # MINIMIZED
  
  # INPUT:
  # alpha               --> VECTOR OVER WHICH MINIMIZATION TAKES PLACE.
  # seal                --> VECTOR OF  FATTY ACID COMPOSITIONS OF SEAL.
  # prey.quantiles      --> MATRIX OF FATTY ACID COMPOSITION OF PREY.
  #                         EACH ROW CONTAINS AN INDIVIDUAL PREY.
  #                         FROM A DIFFERENT SPECIES.
  
  no.zero <- sum(seal == 0.)
  seal[seal == 0.] <- 1e-05
  seal[seal > 0.] <- (1. - no.zero * 1e-05) * seal[seal > 0.]
  sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
  no.zero <- sum(sealhat == 0.)
  sealhat[sealhat == 0.] <- 1e-05
  sealhat[sealhat > 0.] <- (1. - no.zero * 1e-05) * sealhat[sealhat >0.]
  sealhat <- sealhat/sum(alpha)
  return(KLdist(seal, sealhat))
}

KLdist <- function(x, y)
{
  # SIMILAR TO KL BUT REQUIRES TWO VECTORS AS INPUT
  # RETURNS KULBACK-LIEBLER DISTANCE BETWEEN x AND y WHERE 
  # INPUT:
  # x --> VECTOR 
  # y --> vector
  return(sum((x - y) * log(x/y)))
}


optquantile.grad2 <- function(alpha, seal, prey.quantiles) {
  
  # USED IN optquantile (IN CALL TO nlminb) AS THE GRADIENT OF THE OBJECTIVE
  # FUNCTION TO BE MINIMIZED
  
  # INPUT:
  # alpha               --> VECTOR OVER WHICH MINIMIZATION TAKES PLACE.
  # seal                --> VECTOR OF  FATTY ACID COMPOSITIONS OF SEAL.
  # prey.quantiles      --> MATRIX OF FATTY ACID COMPOSITION OF PREY.  EACH ROW
  #                         CONTAINS AN INDIVIDUAL PREY FROM A DIFFERENT SPECIES.
  
  no.zero <- sum(seal == 0.)
  seal[seal == 0.] <- 1e-05
  seal[seal > 0.] <- (1. - no.zero * 1e-05) * seal[seal > 0.]
  sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
  no.zero <- sum(sealhat == 0.)
  sealhat[sealhat == 0.] <- 1e-05
  sealhat[sealhat > 0.] <- (1. - no.zero * 1e-05) * sealhat[sealhat >  0.]
  sumalpha <- sum(alpha)
  sealhat <- sealhat/sumalpha
  alphadum <- sealhat/(sumalpha)
  alphadum.mat <- matrix(rep(alphadum, nrow(prey.quantiles)), byrow = T,
                         nrow(prey.quantiles), length(alphadum))
  prey.quantilesdum <- prey.quantiles/sumalpha - alphadum.mat
  Term1 <- (prey.quantilesdum) %*% t(matrix(log(seal/sealhat), ncol = 
                                              length(seal)))
  Term2 <- (prey.quantilesdum) %*% t(matrix((seal - sealhat)/sealhat, ncol = 
                                              length(seal)))
  return( - Term1 - Term2)
}

optcompdiff.obj <- function(alpha, seal, prey.quantiles) {
  
  # USED IN optquantile (IN CALL TO nlminb) AS THE OBJECTIVE FUNCTION TO BE
  # MINIMIZED
  
  # INPUT:
  # alpha                --> VECTOR OVER WHICH MINIMIZATION TAKES PLACE
  # seal                 --> VECTOR OF  FATTY ACID COMPOSITIONS OF SEAL
  # prey.quantiles       --> MATRIX OF FATTY ACID COMPOSITION OF PREY.
  #                          EACH ROW CONTAINS AN INDIVIDUAL PREY
  #                          FROM A DIFFERENT SPECIES.
  
  no.zero <- sum(seal == 0.)
  seal[seal == 0.] <- 1e-05
  seal[seal > 0.] <- (1. - no.zero * 1e-05) * seal[seal > 0.]
  sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
  no.zero <- sum(sealhat == 0.)
  sealhat[sealhat == 0.] <- 1e-05
  sealhat[sealhat > 0.] <- (1. - no.zero * 1e-05) * sealhat[sealhat > 0.]
  sealhat <- sealhat/sum(alpha)
  return(compdiff(seal, sealhat))
}

optcompdiff.grad <- function(alpha, seal, prey.quantiles) {
  
  # USED IN optquantile (IN CALL TO nlminb) AS THE GRADIENT FUNCTION TO BE
  # MINIMIZED
  
  # INPUT:
  # alpha           --> VECTOR OVER WHICH MINIMIZATION TAKES PLACE
  # seal            --> VECTOR OF  FATTY ACID COMPOSITIONS OF SEAL
  # prey.quantiles  --> MATRIX OF FATTY ACID COMPOSITION OF PREY.EACH ROW CONTAINS AN INDIVIDUAL PREY
  #                     FROM A DIFFERENT SPECIES.
  
  nfa <- ncol(prey.quantiles)
  no.zero <- sum(seal == 0.)
  seal[seal == 0.] <- 1e-05
  seal[seal > 0.] <- (1. - no.zero * 1e-05) * seal[seal > 0.]
  sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
  no.zero <- sum(sealhat == 0.)
  sealhat[sealhat == 0.] <- 1e-05
  sealhat[sealhat > 0.] <- (1. - no.zero * 1e-05) * sealhat[sealhat > 0.]
  sealhat <- sealhat/sum(alpha)
  sumalpha <- sum(alpha)
  alphadum <- sealhat/(sumalpha)
  alphadum.mat <- matrix(rep(alphadum, nrow(prey.quantiles)), byrow = T,
                         nrow(prey.quantiles), length(alphadum))
  prey.quantilesdum <- t(t((prey.quantiles/sumalpha - alphadum.mat))/as.vector(sealhat))
  dist.dum <- log(seal/mean.geometric(seal)) - log(sealhat/mean.geometric(sealhat))
  dist.dum <- matrix(dist.dum, byrow = T, 1., length(seal))
  alphadum.mat2 <- matrix(rep(apply(prey.quantilesdum, 1., sum), nfa), byrow = F,
                          nrow(prey.quantiles), nfa)
  Term1 <- prey.quantilesdum %*% t(dist.dum)
  Term2 <- -1./nfa * alphadum.mat2 %*% t(dist.dum)
  return(-2. * (Term1 + Term2))
}

compdiff <- function(x, bigX)
{
  # COMPUTES THE DIFFERENCE BETWEEN TO VECTORS OF COMPOSITIONAL DATA
  # AS DESCRIBED IN AITCHISON (1992) "MEASURES OF COMPOSITIONAL DIFFERENCE"
  
  return(sum((log(x/mean.geometric(x)) - log(bigX/mean.geometric(bigX)))^
                    2.))
}



#' Generate pseudo predators with ith predator having true diet
#' given by ith row of diet.null.mat.
#'
#' @export
#' @param prey.mat matrix containing a representative FA signature
#'     from each prey group (usually the mean). The first column must
#'     index the prey group.
#' @param diet.null.mat null diet
#' @param cal.mat matrix of calibration coefficients where the \emph{i} th
#'     column is to be used with the \emph{i} th predator.
#' @param fat.cont vector of fat content
#' @param noise proportion of noise to add to the simulation.
#' @param nprey number of prey to sample.
#'
gen.pseudo.seals <- function(prey.mat, diet.null.mat, cal.mat, fat.cont, noise, nprey) {

    futile.logger::flog.info("gen.pseudo.seals()")
    futile.logger::flog.debug("Generating %d pseudo pradators.", nrow(diet.null.mat))
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
            prey.out <- split.prey(prey.mat)
            prey.sim <- prey.out[[1]]
            prey.mod <- prey.out[[2]]

            ## Split fat content for each species
            for(k in 1.:I) {
                fat.split <- split.fatcont(fat.cont[prey.mat[, 1.] == unique(prey.mat[, 1.])[k]])
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
        futile.logger::flog.debug("WARNING IN gen.pseudo.seals:  ONE CALIBATION")
        seals.pseudo <- matrix(rep(0., ns * nFA), byrow = T, ns,nFA)
        
        for(i in 1.:ns) {
            ## Split prey 
            prey.out <- split.prey(prey.mat)
            prey.sim <- prey.out[[1]]
            prey.mod <- prey.out[[2]]

            ## Split fat content for each species
            for(k in 1.:I) {
                fat.split <- split.fatcont(fat.cont[prey.mat[  , 1.] == unique(prey.mat[,1.])[k]])
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
split.prey <- function(prey.mat) {
    
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

split.fatcont <- function(fat.cont.k) {
    
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
##'
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


##---------------------------------
## Old bisection method
##---------------------------------
bisect.beta.find <- function(x1, x2, f1, f2, alpha, par.list, R, p.mat, k)
{
  futile.logger::flog.debug("bisect.beta.find(x1=%f, x2=%f, f(x1)=%f, f(x2)=%f, alpha=%f)", x1, x2, f1, f2, alpha)
    
  # FINDS x1 AND x2 SUCH THAT f(x1)*f(x2)<0
  # NOTE: f1 IS PVALUE AT X1 (WITHOUT ALPHA SUBTRACTED)
  
  int <- x2 - x1
  if((f1 - alpha) == 0.) {
    return(list(x1, x2, f1, f2))
  }
  
  if((f2 - alpha) == 0.) {
    return(list(x1, x2, f1, f2))
  }
  count <- 0.
  while((f1 - alpha) * (f2 - alpha) > 0.) {
    if(round(x2, 6.) == 1.) {
      futile.logger::flog.debug("NEED SMALLER INCREMENT")
      if(count > 3.) {
        futile.logger::flog.debug("FOR ONE SPECIES, NO PVALUES > alpha")
        x1 <- NA
        x2 <- NA
        f1 <- NA
        f2 <- NA
        return(list(x1, x2, f1, f2))
      }
      count <- count + 1.
      int <- int/10.
      x1 <- 0.
      x2 <- x1 + int
      f1 <- beta.pval.k(par.list, R, x1, p.mat, k)
      f2 <- beta.pval.k(par.list, R, x2, p.mat, k)
    }
    else {
      x1 <- x1 + int
      x2 <- x2 + int
      # (x1 <- x2)
      f1 <- f2
      f2 <- beta.pval.k(par.list, R, x2, p.mat, k)
    }
  }
  return(list(x1, x2, f1, f2))
}


bisect.beta <- function(x1, x2, f1, f2, alpha, par.list, R, p.mat, k)
{
  # GIVEN x1 AND x2 SUCH THAT (f(x1)-alpha)*(f(x2)-alpha) < 0, FINDS m such
  # THAT f(m)-alpha = 0.
  # NOTE: f1 and f2 WERE COMPUTED BEFORE alpha WAS SUBTRACTED
  
  if(x2 <= 0.1) {
    tol <- 1e-05
  }
  else {
    tol <- 0.001
  }
  if(abs(x2 - x1) < tol) {
    return(x1)
  }
  if(abs(f1 - alpha) <= tol) {
    return(x1)
  }
  if(abs(f2 - alpha) <= tol) {
    return(x2)
  }
  m <- (x1 + x2)/2.
  fm <- beta.pval.k(par.list, R, m, p.mat, k)
  while(abs(fm - alpha) > tol) {
    if((fm - alpha) * (f1 - alpha) < 0.) {
      x2 <- m
    }
    else {
      x1 <- m
      f1 <- fm
    }
    if(abs(x1 - x2) < tol) {
      m <- (x1 + x2)/2.
      return(m)
    }
    m <- (x1 + x2)/2.
    fm <- beta.pval.k(par.list, R, m, p.mat, k)
  }
  return(m)
}


bisect.beta.lim.orig <- function(alpha1, alpha2, par.list, R, p.mat, k) { 

  futile.logger::flog.debug("bisect.beta.lim.rig()")
    
  # NOVEMBER 15TH, 2011
    
  # SOLVES f(pio)=PVAL(PIO)=alpha1 AND f(PIO)=PVAL(PIO)=alpha2 USING
  # BISECTION.
  # ASSUME alpha1 < alpha2
  # NOTE:  TRIED TO MINIMIZE NUMBER OF TIMES HAVE TO COMPUTE A PVALUE
  # SINCE VERY SLOW.
  # CALLS bisect.parmeth.2find TO FIND x1 AND x2 SUCH THAT f(x1)*f(x2)<0.
  # CALLS bisect.parmeth.2 TO DO BISECTION.
  # COMPUTES 2 UPPER AND 2 LOWER LIMITS
  
  x1 <- 0.
  
  # IF CHANGE x1 AND x2 MUST BE CAREFUL SINCE WE ASSUME 
  # THAT f1.start IS FOR 0 AND f2.start IS FOR 0.10
  
  x2 <- 0.1
  int <- x2 - x1
  f1 <- beta.pval.k(par.list, R, x1, p.mat, k)
  f2 <- beta.pval.k(par.list, R, x2, p.mat, k)
  f1.start <- f1
  
  # FIND x1 and x2 FOR FIRST INTERVAL
  
  f2.start <- f2
  futile.logger::flog.debug("Finding 1")
  
  # FIND x1 AND x2 FOR SECOND INTERVAL
  
  if(f1.start >= alpha1) {
    x1.1.L <- 0.
    x2.1.L <- 0.
    f1.1.L <- f1.start
    f2.1.L <- f1.start
  }
  else {
    val.start <- bisect.beta.find(x1, x2, f1, f2, alpha1,
                                  par.list, R, p.mat, k)
    x1.1.L <- val.start[[1]]
    x2.1.L <- val.start[[2]]
    f1.1.L <- val.start[[1]]
    f2.1.L <- val.start[[2]]
    if(is.na(x1.1.L)) {
      futile.logger::flog.debug("IN BISECT.BETA.LIM - INT TOO SMALL")
      CI.L.1 <- NA
      CI.U.1 <- NA
      CI.L.2 <- NA
      CI.U.2 <- NA
      return(list(CI.L.1, CI.U.1, CI.L.2, CI.U.2))
    }
    int <- x2.1.L - x1.1.L
  }
  futile.logger::flog.debug("Finding 2")
  
  # FIND x1 AND x2 FOR THIRD INTERVAL
  
  if(f1.start >= alpha2) {
    x1.2.L <- 0.
    x2.2.L <- 0.
    f1.2.L <- f1.start
    f2.2.L <- f1.start
  }
  else {
    if((x1.1.L == 0.) && (x2.1.L == 0.)) {
      x1 <- 0.
      x2 <- 0.1
    }
    else {
      x1 <- x1.1.L
      x2 <- x2.1.L
      f1 <- f1.1.L
      f2 <- f2.1.L
    }
    val.start <- bisect.beta.find(x1, x2, f1, f2, alpha2,
                                  par.list, R, p.mat, k)
    x1.2.L <- val.start[[1]]
    x2.2.L <- val.start[[2]]
    f1.2.L <- val.start[[1]]
    f2.2.L <- val.start[[2]]
    if(is.na(x1.2.L)) {
      futile.logger::flog.debug("IN BISECT.BETA.LIM - INT TOO SMALL")
      CI.L.1 <- NA
      CI.U.1 <- NA
      CI.L.2 <- NA
      CI.U.2 <- NA
      return(list(CI.L.1, CI.U.1, CI.L.2, CI.U.2))
    }
    int <- x2.2.L - x1.2.L
  }
  
  futile.logger::flog.debug("Finding 3")
  
  if((x1.2.L == 0.) && (x2.2.L == 0.)) {
    x1 <- 0.
    x2 <- 0.1
    f1 <- f1.start
    if(f1 == alpha2) {
      # MAKE POSITIVE
      f1 <- f1 + 0.0001
    }
    f2 <- f2.start
  }
  else {
    if(f2.2.L == alpha2) {
      # MAKE POSITIVE
      f2.2.L <- f2.2.L + 0.0001
    }
    x1 <- x2.2.L
    x2 <- x2.2.L + int
    f1 <- f2.2.L
    f2 <- beta.pval.k(par.list, R, x2, p.mat, k)
  }
  val.start <- bisect.beta.find(x1, x2, f1, f2, alpha2, par.list,
                                R, p.mat, k)
  x1.2.U <- val.start[[1]]
  x2.2.U <- val.start[[2]]
  f1.2.U <- val.start[[1]]
  f2.2.U <- val.start[[2]]
  if(is.na(x1.2.U)) {
    futile.logger::flog.debug("IN BISECT.BETA.LIM - INT TOO SMALL")
    CI.L.1 <- NA
    CI.U.1 <- NA
    CI.L.2 <- NA
    CI.U.2 <- NA
    return(list(CI.L.1, CI.U.1, CI.L.2, CI.U.2))
  }
  # FIND x1 AND x2 FOR FOURTH INTERVAL
  int <- x2.2.U - x1.2.U
  
  futile.logger::flog.debug("Finding 4")
  
  x1 <- x1.2.U
  x2 <- x2.2.U
  f1 <- f1.2.U
  f2 <- f2.2.U
  val.start <- bisect.beta.find(x1, x2, f1, f2, alpha1, par.list,
                                R, p.mat, k)
  x1.1.U <- val.start[[1]]
  x2.1.U <- val.start[[2]]
  f1.1.U <- val.start[[1]]
  f2.1.U <- val.start[[2]]
  if(is.na(x1.1.U)) {
    futile.logger::flog.debug("IN BISECT.BETA.LIM - INT TOO SMALL")
    CI.L.1 <- NA
    CI.U.1 <- NA
    CI.L.2 <- NA
    CI.U.2 <- NA
    return(list(CI.L.1, CI.U.1, CI.L.2, CI.U.2))
  }
  futile.logger::flog.debug("x1 and x2 to use in lowest limit")
  futile.logger::flog.debug(c(x1.1.L, x2.1.L))
  CI.L.1 <- bisect.beta(x1.1.L, x2.1.L, f1.1.L, f2.1.L, alpha1,
                        par.list, R, p.mat, k)
  futile.logger::flog.debug("x1 and x2 to use in second limit")
  futile.logger::flog.debug(c(x1.2.L, x2.2.L))
  CI.L.2 <- bisect.beta(x1.2.L, x2.2.L, f1.2.L, f2.2.L, alpha2,
                        par.list, R, p.mat, k)
  futile.logger::flog.debug("x1 and x2 to use in third limit")
  futile.logger::flog.debug(c(x1.2.U, x2.2.U))
  CI.U.2 <- bisect.beta(x1.2.U, x2.2.U, f1.2.U, f2.2.U, alpha2,
                        par.list, R, p.mat, k)
  futile.logger::flog.debug("x1 and x2 to use in fourth limit")
  futile.logger::flog.debug(c(x1.1.U, x2.1.U))
  CI.U.1 <- bisect.beta(x1.1.U, x2.1.U, f1.1.U, f2.1.U, alpha1,
                        par.list, R, p.mat, k)
  return(list(CI.L.1, CI.U.1, CI.L.2, CI.U.2))
}


#' Find simultaneous and individual confidence intervals for diet
#' proportions of a single prey species i.e. solve f(pio) = PVAL(pio)
#' = alpha1 and f(pio) = PVAL(pio) = alpha2 using bisection. 
#' 
#' Assumptions: alpha1 < alpha2
#' 
#' Note:  Tried to minimize number of times have to compute a pvalue
#' since very slow.
#' 
#' @param alpha1 simultaneous confidence level
#' @param alpha2 individual confidence level
#' @param par.list a list of R.p lists of I beta distribution parameters phi
#'     and theta that define diet proportion estimates for each of the
#'     prey species. Effectively R.p beta distibutions for each of the
#'     I prey species from which we bootstrap to calculate p-values.
#' @param R number of bootstrap replicates to use in p-value estimation.
#' @param p.mat predator diet estimates.
#' @param k prey species index 1..I
#' @return upper and lower limits for individual and simultaneous
#'     intervals respectively. 
#' 
bisect.beta.lim <- function(alpha1, alpha2, par.list, R, p.mat, k) { 
    
    futile.logger::flog.debug("bisect.beta.lim()")
    

    ##----------------------------------------------------------
    ## L1: return lower limit if there is an issue finding roots
    ##----------------------------------------------------------
    futile.logger::flog.debug("Finding L1 (%f)", alpha1)
    CI.L.1 = tryCatch({ uniroot.beta(x1 = 0,
                                     x2 = colMeans(p.mat)[k],
                                     alpha = alpha1,
                                     par.list,
                                     R,
                                     p.mat,
                                     k) },
                      error = function(err) {
                          futile.logger::flog.warn("Issue finding L1 roots. Using lower limit: %s", err)
                          return(0) })
    
    ##----------------------------------------------------------
    ## L2: return lower limit if there is an issue finding roots
    ##----------------------------------------------------------
    futile.logger::flog.debug("Finding L2 (%f)", alpha2)
    CI.L.2 = tryCatch({ uniroot.beta(0,
                                     colMeans(p.mat)[k],
                                     alpha2,
                                     par.list,
                                     R,
                                     p.mat,
                                     k) },
                      error = function(err) {
                          futile.logger::flog.warn("Issue finding L2 roots. Using lower limit: %s", err)
                          return(0) })

    ##----------------------------------------------------------
    ## U2: return upper limit if there is an issue finding roots   
    ##----------------------------------------------------------
    futile.logger::flog.debug("Finding U2 (%f)", alpha2)
    CI.U.2 = tryCatch({ uniroot.beta(x1 = colMeans(p.mat)[k], x2 = 1, alpha = alpha2, par.list, R, p.mat, k) },
                      error = function(err) {
                          futile.logger::flog.warn("Issue finding roots. Using upper limit: %s", err)
                          return(1) })

    ##----------------------------------------------------------
    ## U1: return upper limit if there is an issue finding roots
    ##----------------------------------------------------------
    futile.logger::flog.debug("Finding U1 (%f)", alpha1)
    CI.U.1 = tryCatch({ uniroot.beta(colMeans(p.mat)[k], 1, alpha1, par.list, R, p.mat, k) },
                      error = function(err) {
                          futile.logger::flog.warn("Issue finding roots. Using upper limit: %s", err)
                          return(1) })

    return(list(CI.L.1, CI.U.1, CI.L.2, CI.U.2))
}



#' Calculate p-value of a given prey type diet proportion under the
#' predator diets and estimated diet distribution provided.
#'
#' @param par.list a list of R.p lists of I beta distribution parameters phi
#'     and theta that define diet proportion estimates for each of the
#'     prey species. Effectively R.p beta distibutions for each of the
#'     I prey species which we bootstrap to calculate p-values. 
#' @param R number of bootstrap replicates to use in p-value estimation.
#' @param diet.test.k diet proportion for which we want the p-value
#' @param p.mat predator diet estimates.
#' @param k prey species index 1..I 
#' @return p-value
#' 
beta.pval.k <- function(par.list, R, diet.test.k, p.mat, k) { 
    
    ## NOVEMBER 15TH, 2011
    I <- ncol(p.mat)
    ns <- nrow(p.mat)
    T.obs <- abs(p.beta(p.mat[,k]) - diet.test.k)
    R.s <- length(par.list)
    R.p <- length(par.list[[1.]])
    pval.vec <- rep(0., R.s)

    futile.logger::flog.debug("beta.pval.k(T.obs=%f, diet.test.k=%f, R=%d, R.s=%d, R.p=%d", T.obs, diet.test.k, R, R.s, R.p)
    
    for(r.s in 1.:R.s) {
        pval.vec.prey <- rep(0., R.p)
        
        for(r.p in 1.:R.p) {
            ## Thetas
            prob.zero.vec <- par.list[[r.s]][[r.p]][1.,  ]

            ## Phis 
            phi.est <- par.list[[r.s]][[r.p]][2.,  ]
            
            if(round(diet.test.k, 6.) == 0.) {
                if(round(prob.zero.vec[k], 6.) == 1.) {
                    pval.vec.prey[r.p] <- runif(1.)
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
                T.star <- Tstar.beta(p.mat[, k], 
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
        ## Averate p-values over R.p ?
        pval.vec[r.s] <- mean(pval.vec.prey)
    }

    ## Average p-values over R.s ?
    pval <- mean(pval.vec)
    futile.logger::flog.debug("    beta.pval.k(%.12f) = %.12f", diet.test.k, pval)
    return(pval)
}


#' Generate bootstrap replicates of diet proportion estimates for
#' given prey species
#'
#' @param p.k list of predator diet proportions for prey species k
#' @param theta.boot 
#' @param mu.null 
#' @param phi.boot 
#' @param R.boot number of bootstrap replicates to create
#' @return list of diet proportion replicates

Tstar.beta <- function(p.k, theta.boot, mu.null, phi.boot, R.boot) {
    futile.logger::flog.trace("Tstar.beta()")
    data.para <- boot::boot(data = p.k,
                            statistic = p.beta,
                            R = R.boot, 
                            sim = "parametric",
                            ran.gen = data.sim.beta,
                            mle = list(length(p.k), mu.null, phi.boot, theta.boot))
                            
    return(data.para$t)
}


#' Bootstrap ran.gen function: 
data.sim.beta <- function(data, mle) {
    if (mle[[4]]==0) { mle[[4]] <- 0.00001 }
    d <- gamlss.dist::rBEZI(mle[[1]], mle[[2]], mle[[3]], mle[[4]])
    d[d>0.999] <- 0.999    
    return(d)
}


#' Bootstrap statistic function: in this case it is the mean
#' (meanBEZI) of the fitted (gamlss) ZIB distribution.
#'
#' @param data data
#' @return the expected value of the response for a fitted model
#' 
p.beta <- function(data) {
    
    ## FOUND ERRORS IN HOLLY'S CODE SO RE-WROTE
    p.vec <- data
    ns <- length(p.vec)
    
    no.zeros <- sum(round(p.vec,4)==0)
    
    if ( (no.zeros==ns) || (no.zeros==(ns-1)) )
    { 
        return(0)
    }

    futile.logger::flog.trace('gamlss::gamlss()')
    est.gam <- gamlss::gamlss(p.vec~1,
                      sigma.formula=~1, 
                      nu.formula=~1,
                      control=(gamlss::gamlss.control(n.cyc=1000,trace=F)),
                      family = gamlss.dist::BEZI)
        
    return(gamlss.dist::meanBEZI(est.gam)[1])
}


#' Use uniroot() to find roots and compare with bisection outcome
uniroot.beta <- function(x1, x2, alpha, par.list, R, p.mat, k)
{

    futile.logger::flog.debug("uniroot.beta(x1=%.12f, x2=%.12f, alpha=%.6f)", x1, x2, alpha)
    if (x2 <= 0.1) { tol <- 1e-05 } else { tol <- 0.001 }
    r = uniroot(function(x) { beta.pval.k(par.list, R, x, p.mat, k) - alpha }, c(x1, x2), tol = tol, trace=1)
    futile.logger::flog.debug("Uniroot found root f(%f) = %f in %d iterations", r$root, r$f.root, r$iter)
    return(r$root)
}


#' Calculate bias correction for confidence intervals from \code{\link{beta.meths.CI()}}.
#'
#' @export
#' @param p.mat matrix containing the fatty acid signatures of the predators.
#' @param prey.mat matrix containing a representative fatty acid signature
#' @param cal.mat matrix of calibration factors where the \emph{i} th
#'     column is to be used with the \emph{i} th predator. If modelling is
#'     to be done without calibration coefficients, simply pass a
#'     vector or matrix of ones. 
#' @param fat.cont prey fat content
#' @param R.bias botstrap replicates
#' @param noise noise
#' @param nprey number of prey
#' @param specify.noise noise
#' @param dist.meas distance measure
#' @param ext.fa subset of FA's to use.
#' @return Row 1 is Lambda CI, row 2 is Lambda skew, and row 3 is Beta CI
#' 
bias.all <- function(p.mat,
                     prey.mat,
                     cal.mat,
                     fat.cont, 
                     R.bias,
                     noise,
                     nprey,
                     specify.noise,
                     dist.meas,
                     ext.fa = ext.common.fa.list) {
        
    ## Returned : Row 1 is Lamda CI, Row 2 is Lamda Skew CI, Row 3 is Beta CI
    lambda.est <- matrix(0, nrow(p.mat), ncol(p.mat))
    lambda.skew <- matrix(0, nrow(p.mat), ncol(p.mat))
    beta.est <- matrix(0, nrow(p.mat), ncol(p.mat))
    
    bias1 <- matrix(0, nrow(p.mat), ncol(p.mat))
    bias2 <- matrix(0, nrow(p.mat), ncol(p.mat))
    bias3 <- matrix(0, nrow(p.mat), ncol(p.mat))

    for (i in 1:nrow(p.mat)) {
        
        diet.true <- as.matrix(p.mat[i,])
        I <- length(unique(prey.mat[,1]))
        
        if ( !((nrow(cal.mat) == 1.) || (ncol(cal.mat) == 1.)) ) {
            
            ind.sim <- sample(seq(1., ncol(cal.mat), 1.), R.bias, replace = T)
            cal.sim <- as.matrix(cal.mat[, ind.sim])
            
            cal.mod <- matrix(rep(0., nrow(cal.mat) * R.bias), byrow = T, nrow(cal.mat), R.bias)
        }
        
        p.mat.seal <- matrix(rep(0., R.bias * I), byrow = T, R.bias, I)
        
        ## GENERATE ns PSEUDO SEALS WITH CALIBRATION FACTORS IN cal.sim
        ## AND CALCULATE ns DIET ESTIMATES USIN cal.mod.  INCORPORATE FAT
        ## CONTENT
        seals.pseudo <- matrix(rep(0., R.bias * (ncol(prey.mat) - 1.)),byrow = T, R.bias,
        (ncol(prey.mat) - 1.))
        fat.sim <- rep(0,I)
        fat.mod <- rep(0,I)
        
        for (n in 1.:R.bias) {
            
            prey.split <- split.prey(prey.mat)
            prey.sim <- prey.split[[1]]
            prey.mod <- prey.split[[2]]
            
            prey.mod.ext <- as.data.frame(prey.mod)[c(dimnames(prey.mod)[[2]][1],ext.fa)]
            
            prey.mod.ext[, -1.] <- prey.mod.ext[, -1.]/apply(prey.mod.ext[, -1.], 1., sum)
            for (k in 1.:I) {
                
                fat.split <- split.fatcont(fat.cont[prey.mat[  , 1.] == unique(prey.mat[, 1.])[k]])
                fat.sim[k] <- fat.split[1]
                fat.mod[k] <- fat.split[2]
            } # END k
            
            if ((nrow(cal.mat) == 1.) || (ncol(cal.mat) == 1.)) {
                
                
                seals.pseudo[n,  ] <- pseudo.seal(prey.sim, diet.true, noise, nprey,
                                                  cal.mat, fat.sim, specify.noise)
                seal.to.use <- matrix(seals.pseudo[n,  ], nrow = 1.)
                seal.to.use <- as.data.frame(seal.to.use)
                dimnames(seal.to.use)[[2.]] <- dimnames(prey.mat[,-1.])[[2.]]
                p.mat.seal[n,  ] <- p.prey(seal.to.use,MEANmeth(prey.mod.ext),cal.mat,
                                           dist.meas, fat.mod,ext.fa=ext.fa)
                
            } else {
                seals.pseudo[n,  ] <- pseudo.seal(prey.sim, diet.true, noise, nprey,
                                                  cal.sim[,n], fat.sim, specify.noise)
                seal.to.use <- matrix(seals.pseudo[n,  ], nrow = 1.)
                seal.to.use <- as.data.frame(seal.to.use)
                dimnames(seal.to.use)[[2.]] <- dimnames(prey.mat[,-1.])[[2.]]
                cal.mod[, n] <- apply(cal.mat[,  - ind.sim[n]], 1.,mean)
                p.mat.seal[n,  ] <- p.prey(seal.to.use,MEANmeth(prey.mod.ext),cal.mod[,n],
                                           dist.meas, fat.mod,ext.fa=ext.fa)
            }
        } # END n

        ## Lamda skew
        for (j in 1:I) {
            p.vec <- p.mat.seal[,j]
            non.zeros <- p.mat.seal[,j][p.mat.seal[,j] != 0]
            
            #' Lamda
            trans.data <-  log(non.zeros/(1-non.zeros))
            theta <- (R.bias-length(non.zeros))/R.bias
            mu <- mean(trans.data)
            if (theta==1)
            { lambda.est[i,j]<-0
            } else {
                
                lambda.est[i,j] <- (1-theta)*exp(mu)/(1+exp(mu))
            }
            
            #' Lambda Skew
            #' NEED TO FIX THE SKEW FUNCTIONS.
            #'    lambda.skew[i,j] <- p.skew(p.vec)
            lambda.skew[i,j] <- 0
            
            #' Beta
            beta.est[i,j] <- p.beta(p.vec)
            
        } # END j

        bias1[i,] <- lambda.est[i,] - diet.true
        bias2[i,] <- lambda.skew[i,] - diet.true
        bias3[i,] <- beta.est[i,] - diet.true
        
    } # END i

    
    bias1 <- apply(bias1,2,median)
    bias2 <- apply(bias2,2,median)
    bias3 <- apply(bias3,2,median)
    
    return(rbind(bias1, bias2, bias3))
}


#' This is a hidden QFASA function so duplicate here...
mean.geometric <- function(x) {
    D <- length(x)
    return(prod(x)^(1./D))
}
