#' QFASA: A package for Quantitative Fatty Acid Signature Analysis
#'
#' Accurate estimates of the diets of predators are required
#' in many areas of ecology, but for many species current methods are
#' imprecise, limited to the last meal, and often biased. The diversity
#' of fatty acids and their patterns in organisms, coupled with the
#' narrow limitations on their biosynthesis, properties of digestion in
#' monogastric animals, and the prevalence of large storage reservoirs of
#' lipid in many predators, led us to propose the use of quantitative
#' fatty acid signature analysis (QFASA) to study predator diets.
#'
#' @aliases QFASA-package
#' @name QFASA
#'
NULL

#' Returns QFASA diet estimates corresponding to a sample of predators.
#'
#' Computes the diet estimate for each predator in \emph{predator.mat} using
#' either the Kullback-Leibler Distance (KL), the Aitchison Distance
#' (AIT) or the Chi-Square Distance (CS).
#'
#' @export
#' @param predator.mat matrix containing the FA signatures of the predators.
#' @param prey.mat matrix containing a representative FA signature
#'     from each prey group (usually the mean). The first column must
#'     index the prey group. Note can use function \emph{MEANmeth} to calculate the means.
#' @param cal.mat matrix of calibration factors where the \emph{i} th
#'     column is to be used with the \emph{i} th predator. If modelling is to be done without
#'     calibration coefficients, simply pass a vector or matrix of ones.
#' @param dist.meas distance measure to use for estimation: 1=KL,
#'     2=AIT or 3=CS
#' @param gamma parameter required for calculations using CS distance
#'     (passed to CS.obj). Currently being set to 1.
#' @param FC vector of fat content of length equal to the number of prey groups
#'           or species.
#' @param start.val initial vector of parameters to be optimized
#' @param ext.fa subset of fatty acids to be used to obtain QFASA diet estimates.
#'
#' @details Before carrying out an analysis using QFASA, rows of prey database must be normalized to sum to one.
#'          See Example for code that extracts a subset of FAs and then normalizes the prey database signatures.
#'
#' @return A list with components:
#' \item{Diet Estimates}{A matrix of the diet estimates for each predator where each row corresponds to a predator and the columns to prey species. The estimates are expressed as proportions summing to one.}
#' \item{Additional Measures}{For each predator for which a diet estimate was obtained:}
#' \item{ModFAS}{the value of the modelled fatty acid. These are expressed as proportions summing to one.}
#' \item{DistCont}{The contribution of each fatty acid to the final minimized distance.}
#' \item{PropDistCont}{The contribution of each fatty acid to the final minimized distance as a proportion of the total.}
#' \item{MinDist}{The final minimized distance.}
#'
#' @references Iverson, Sara J., Field, Chris, Bowen, W. Don and
#' Blanchard, Wade (2004) Quantitative Fatty Acid Signature Analysis: A New
#' Method of Estimating Predator Diets. Ecological Monographs, 74(2), 211-235
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
#'  ## Fat Content
#'
#'  FC = preyFAs[,c(2,3)]
#'  FC = as.vector(tapply(FC$lipid,FC$Species,mean,na.rm=TRUE))
#'
#'  ## Calibration Coefficients
#'  data(CC)
#'  cal.vec = CC[,2]
#'  cal.mat = replicate(npredators, cal.vec)
#'
#'  ## Run QFASA
#'  Q = p.QFASA(predator.matrix,
#'              prey.matrix,
#'              cal.mat,
#'              dist.meas = 1,
#'              gamma=1,
#'              FC,
#'              start.val = rep(1,nrow(prey.matrix)),
#'              fa.set)
#'
#' ## Diet Estimates
#' DietEst = Q$'Diet Estimates'
#'
p.QFASA <- function(predator.mat,
                    prey.mat,
                    cal.mat,
                    dist.meas,
                    gamma = 1,
                    FC = rep(1, nrow(prey.mat)),
                    start.val = rep(0.99999, nrow(prey.mat)),
                    ext.fa) {

    seal.mat <- predator.mat

    # Cast inputs to prevent broadcasting
    seal.mat = as.matrix(seal.mat)
    prey.mat = as.matrix(prey.mat)
    cal.mat = as.matrix(cal.mat)

    # CALIBRATING SEAL FA SIGNATURES AND THEN EXTRACTING EXTENDED DIETARY FAS
    if ((is.vector(cal.mat)) || (nrow(cal.mat) == 1.) || (ncol(cal.mat ) == 1.)) {

        ## IF ONLY ONE SEAL
        seal.mat <- t(t(seal.mat)/as.vector(unlist(cal.mat)))
        seal.mat <- as.data.frame(seal.mat)[ext.fa]
        seal.mat <- seal.mat/apply(seal.mat, 1., sum)
        seal.mat <- as.matrix(seal.mat)
    }
    else {
        seal.mat <- seal.mat/t(cal.mat)
        seal.mat <- as.data.frame(seal.mat)[ext.fa]
        seal.mat <- seal.mat/apply(seal.mat, 1., sum)
        seal.mat <- as.matrix(seal.mat)
    }

    I <- nrow(prey.mat)
    ns <- nrow(seal.mat)
    p.mat <- matrix(rep(0, nrow(prey.mat) * nrow(seal.mat)),
                    byrow = T,
                    nrow(seal.mat),
                    nrow(prey.mat))
    more.list <- vector("list",ns)

    # USING SAME STARTING VECTOR FOR EACH OF THE ns SEALS
    if (!(is.matrix(start.val))) {
        start.val <- matrix(rep(start.val, ns), byrow = T, ns, I)
    }

    if (dist.meas == 1) { # KL Distance
        for (i in 1.:nrow(seal.mat)) {
            p.all <- Rsolnp::solnp(pars = start.val[i,  ],
                                   fun  = KL.obj,
                                   predator = seal.mat[i,  ],
                                   prey.quantiles = prey.mat,
                                   eqfun = QFASA.const.eqn,
                                   eqB = 1,
                                   LB = rep(0, nrow(prey.mat)),
                                   UB = rep(0.999999, nrow(prey.mat)),
                                   control=list(trace=0))

            if (p.all$convergence!=0) {
                cat("WARNING: DID NOT CONVERGE")
            }

            p.mat[i,  ] <- p.all$par
            more.list[[i]] <- KL.more(p.mat[i,],seal.mat[i,],prey.mat)
        }

    }
    else if (dist.meas==2) { # AIT distance
        for (i in 1.:nrow(seal.mat)) {
            p.all <- Rsolnp::solnp(pars = start.val[i,  ],
                                   fun  = AIT.obj,
                                   predator = seal.mat[i,  ],
                                   prey.quantiles = prey.mat,
                                   eqfun = QFASA.const.eqn,
                                   eqB=1,
                                   LB = rep(0, nrow(prey.mat)),
                                   UB = rep(0.999999, nrow(prey.mat)),
                                   control=list(trace=0))

            if (p.all$convergence!=0) {
                cat("WARNING: DID NOT CONVERGE")
            }

            p.mat[i,  ] <- p.all$par
            more.list[[i]] <- AIT.more(p.mat[i,],seal.mat[i,],prey.mat)
        }
    }
    else { # CS distance
        for (i in 1.:nrow(seal.mat)) {
            p.all <- Rsolnp::solnp(pars = start.val[i,  ],
                                   fun  = CS.obj,
                                   predator = seal.mat[i,  ],
                                   prey.quantiles = prey.mat,
                                   gamma = gamma,
                                   eqfun= QFASA.const.eqn,
                                   eqB = 1,
                                   LB = rep(0, nrow(prey.mat)),
                                   UB = rep(0.999999, nrow(prey.mat)),
                                   control=list(trace=0))

            if (p.all$convergence!=0) {
                cat("WARNING: DID NOT CONVERGE")
            }

            p.mat[i,  ] <- p.all$par
            more.list[[i]] <- CS.more(p.mat[i,],seal.mat[i,],prey.mat,gamma)

        }

    }

    if (is.matrix(FC))
    { FC.mat <- FC
    } else {
        FC.mat <- matrix(rep(FC, nrow(seal.mat)), byrow = T, nrow(seal.mat), I)

    }

    p.mat <- p.mat/FC.mat
    p.mat <- p.mat/apply(p.mat, 1, sum)
    colnames(p.mat) <- as.vector(rownames(prey.mat)) #add column names
    out.list <- list(p.mat,more.list)
    names(out.list) <- c("Diet Estimates","Additional Measures")
    return(out.list)
}


#' Returns the distance between two compositional vectors using Aitchison's distance measure.
#'
#' @export
#' @param x.1 compositional vector
#' @param x.2 compositional vector
#' @references Aitchison, J., (1992) On criteria for measures of compositional difference. Mathematical Geology, 24(4), pp.365-379.
#' @references Connie Stewart (2017) An approach to measure distance between compositional diet estimates containing essential zeros, Journal of Applied Statistics, 44:7, 1137-1152, DOI: 10.1080/02664763.2016.1193846
#'
AIT.dist <- function(x.1, x.2) {
    return(sqrt(sum((log(x.1/mean_geometric(x.1)) - log(x.2/mean_geometric(x.2)))^2.)))
}


#' Used to provide additional information on various model components
#' evaluated at the optimal solution, i.e., using the QFASA diet estimates and
#' Aitchison distance measure.
#'
#' @export
#' @param alpha compositional QFASA diet estimate.
#' @param predator fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#'
AIT.more <- function(alpha, predator, prey.quantiles) {

    seal <- predator
    no.zero <- sum(seal == 0.)
    seal[seal == 0.] <- 1e-05
    seal[seal > 0.] <- (1. - no.zero * 1e-05) * seal[seal > 0.]

    sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
    no.zero <- sum(sealhat == 0.)
    sealhat[sealhat == 0.] <- 1e-05
    sealhat[sealhat > 0.] <- (1. - no.zero * 1e-05) * sealhat[sealhat > 0.]

    AIT.sq.vec <-
        ( log(seal/mean_geometric(seal)) - log(sealhat/mean_geometric(sealhat)) )^2

    dist <- (sum(AIT.sq.vec))^(1/2)

    out.list <- list(sealhat,AIT.sq.vec,AIT.sq.vec/sum(AIT.sq.vec),dist)
    names(out.list) <- c("ModFAS","DistCont", "PropDistCont","MinDist")

    return(out.list)
}


#' Used in \code{solnp()} as the objective function to be minimized when
#' Aitchison distance measure is chosen.
#'
#' @export
#' @param alpha vector over which minimization takes place.
#' @param predator fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#'
AIT.obj <- function(alpha, predator, prey.quantiles) {

    seal <- predator
    no.zero <- sum(seal == 0.)
    seal[seal == 0.] <- 1e-05
    seal[seal > 0.] <- (1. - no.zero * 1e-05) * seal[seal > 0.]

    sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
    no.zero <- sum(sealhat == 0.)
    sealhat[sealhat == 0.] <- 1e-05
    sealhat[sealhat > 0.] <- (1. - no.zero * 1e-05) * sealhat[sealhat > 0.]

    return(AIT.dist(seal, sealhat))
}


#' Returns the distance between two compositional vectors using the chi-square distance.
#'
#' @export
#' @param x.1 compositional vector
#' @param x.2 compositional vector
#' @param gamma power transform taken to be 1.
#' @references Stewart, C., Iverson, S. and Field, C. (2014) Testing for a change in
#' diet using fatty acid signatures.  Environmental and Ecological Statistics 21, pp. 775-792.
#' @references Connie Stewart (2017) An approach to measure distance between compositional diet estimates containing essential zeros, Journal of Applied Statistics, 44:7, 1137-1152, DOI: 10.1080/02664763.2016.1193846
chisq.dist <- function(x.1, x.2, gamma) {

    nfa <- length(x.1)

    y.1 <- x.1^(gamma)
    y.1 <- y.1/sum(y.1)

    y.2 <- x.2^(gamma)
    y.2 <- y.2/sum(y.2)

    d.sq <- (y.1-y.2)^2
    c.vec <- y.1+y.2

    if ( any(d.sq!=0) ) {

        d.sq[d.sq!=0] <- d.sq[d.sq!=0]/c.vec[d.sq!=0]
    }

    CS.dist <- 1/gamma*sqrt(2*nfa)*sqrt(sum(d.sq))
    return(CS.dist)
}


#' Used to provide additional information on various model components
#' evaluated at the optimal solution, i.e., using the QFASA diet estimates and
#' chi-square distance measure.
#'
#' @export
#' @param alpha compositional QFASA diet estimate.
#' @param predator fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#' @param gamma power transform exponent (see \code{chisq.dist()}).
#'
CS.more <- function(alpha, predator, prey.quantiles, gamma) {

    seal <- predator
    sealhat <- t(as.matrix(alpha)) %*% prey.quantiles

    nfa <- length(seal)
    y.1 <- seal^(gamma)
    y.1 <- y.1/sum(y.1)
    y.2 <- sealhat^(gamma)
    y.2 <- y.2/sum(y.2)

    d.sq <- (y.1-y.2)^2
    c.vec <- y.1+y.2

    if ( any(d.sq!=0) ) {

        d.sq[d.sq!=0] <- d.sq[d.sq!=0]/c.vec[d.sq!=0]
    }

    CS.vec.sq <- d.sq
    dist <- 1/gamma*sqrt(2*nfa)*sqrt(sum(d.sq))

    out.list <- list(sealhat,CS.vec.sq,CS.vec.sq/sum(CS.vec.sq),dist)
    names(out.list) <- c("ModFAS","DistCont", "PropDistCont","MinDist")
    return(out.list)
}


#' Used in \code{solnp()} as the objective function to be minimized when
#' chi-square distance measure is chosen. Unlike \code{AIT.obj()} and \code{KL.obj()},
#' does not require modifying zeros.
#'
#' @export
#' @param alpha vector over which minimization takes place.
#' @param predator fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#' @param gamma power transform exponent (see \code{chisq.dist()}).
CS.obj <- function(alpha, predator, prey.quantiles, gamma){

    seal <- predator
    sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
    return(chisq.dist(seal, sealhat, gamma))
}


#' Returns the distance between two compositional vectors using Kullback--Leibler
#' distance measure.
#'
#' @export
#' @param x.1 compositional vector
#' @param x.2 compositional vector
#' @references S.J. Iverson, C. Field, W.D. Bowen, and W. Blanchard (2004)  Quantitative
#' fatty acid signature analysis: A new method of estimating predator diets, Ecological
#' Monographs 72, pp. 211-235.
#'
KL.dist <- function(x.1, x.2) {
    return(sum((x.1 - x.2) * log(x.1/x.2)))
}


#' Used to provide additional information on various model components
#' evaluated at the optimal solution, i.e., using the QFASA diet estimates and
#' Kullback-Leibler distance measure.
#'
#' @export
#' @param alpha compositional QFASA diet estimate.
#' @param predator fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#'
KL.more <- function(alpha, predator, prey.quantiles) {

    no.zero <- sum(predator == 0.)
    predator[predator == 0.] <- 1e-05
    predator[predator > 0.] <- (1. - no.zero * 1e-05) * predator[predator > 0.]

    predatorhat <- t(as.matrix(alpha)) %*% prey.quantiles
    no.zero <- sum(predatorhat == 0.)
    predatorhat[predatorhat == 0.] <- 1e-05
    predatorhat[predatorhat > 0.] <- (1. - no.zero * 1e-05) * predatorhat[predatorhat >0.]

    KL.vec <- (predator - predatorhat) * log(predator/predatorhat)

    dist <- sum(KL.vec)

    out.list <- list(predatorhat,KL.vec,KL.vec/sum(KL.vec),dist)
    names(out.list) <- c("ModFAS","DistCont", "PropDistCont","MinDist")

    return(out.list)
}


#' Used in \code{solnp()} as the objective function to be minimized when
#' Kullback--Leibler distance measure is chosen.
#'
#' @export
#' @param alpha vector over which minimization takes place.
#' @param predator fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#'
KL.obj <- function(alpha, predator, prey.quantiles) {

    predator <- predator
    no.zero <- sum(predator == 0.)
    predator[predator == 0.] <- 1e-05
    predator[predator > 0.] <- (1. - no.zero * 1e-05) * predator[predator > 0.]
    predatorhat <- t(as.matrix(alpha)) %*% prey.quantiles
    no.zero <- sum(predatorhat == 0.)
    predatorhat[predatorhat == 0.] <- 1e-05
    predatorhat[predatorhat > 0.] <- (1. - no.zero * 1e-05) * predatorhat[predatorhat >0.]
    return(KL.dist(predator, predatorhat))
}


#' Returns the geometric mean of a compositional vector
#'
#' @param x compositional vector
#'
mean_geometric <- function(x) {
    D <- length(x)
    return(prod(x)^(1./D))
}


#' Returns the multivariate mean FA signature of each prey group
#' entered into the QFASA model. Result can be passed to prey.mat in
#' \code{p.QFASA()}.
#'
#' @export
#' @param prey.mat matrix containing the FA signatures of the
#'     prey. The first column indexes the prey group.
#'
MEANmeth <- function(prey.mat) {
    prey.mat[,-1] <- prey.mat[,-1]/apply(prey.mat[,-1],1,sum)
    prey.means <- apply(prey.mat[, -1], 2, tapply, prey.mat[, 1], mean)
    return(prey.means)
}


#' Returns \code{sum(alpha)} and used in \code{solnp()}.
#'
#' @param alpha vector over which minimization takes place.
#' @param predator fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#' @param gamma power transform exponent (see chisq.dist).
#'
QFASA.const.eqn <- function(alpha, predator, prey.quantiles, gamma) {
    return(sum(alpha))
}


#' Test for a difference between two independent samples of compositional data.
#' Zeros of any type are allowed.
#' @export
#' @param compdata.1 sample of compositional data.
#' @param compdata.2 sample of compositional data.
#' @param R number of bootstrap samples, default is 500.
#'
#' @return  p-value obtained through a multivariate permutation test with test statistic based
#' on chi-square distances.
#'
#' @examples
#'
#' ## Prey
#' data(preyFAs)
#'
#' ## Capelin FA sig
#' capelin.sig=preyFAs[preyFAs$Species=="capelin",4:(ncol(preyFAs))]
#' capelin.sig=capelin.sig/apply(capelin.sig,1,sum)
#'
#' ## Sandlance FA sig
#' sandlance.sig=preyFAs[preyFAs$Species=="sandlance",4:(ncol(preyFAs))]
#' sandlance.sig=sandlance.sig/apply(sandlance.sig,1,sum)
#'
#' # Note: uncomment examples to run. CRAN tests fail because execution time > 5 seconds
#' # testfordiff.ind.pval(as.matrix(capelin.sig),as.matrix(sandlance.sig))
#'
#'
#' @references Stewart, C., Iverson, S. and Field, C. (2014) Testing for a change in
#' diet using fatty acid signatures.  Environmental and Ecological Statistics 21, pp. 775-792.
#'
testfordiff.ind.pval <- function(compdata.1, compdata.2, R=500) {

  compdata.1 <- as.matrix(compdata.1)
  ns1 <- nrow(compdata.1)
  compdata.2 <- as.matrix(compdata.2)
  boot.out <- testfordiff.ind.boot(rbind(compdata.1, compdata.2), ns1, R)
  T.orig <- boot.out$t0
  T.vec <- boot.out$t

  pval.vec <- T.vec >= T.orig

  pval <- mean(pval.vec)

  return(list(T.orig, T.vec, pval))

}


#' Called by \code{testfordiff.ind.pval()}.
#' @export
#' @param data sample of compositional data
#' @param ns1 sample size of compdata.1
#' @param R number of bootstrap samples.  default is 500.
#'

testfordiff.ind.boot <- function(data, ns1, R) {

  data.boot <- boot::boot(data = data, statistic = testfordiff.ind.boot.fun,
                    ns1 = ns1, sim = "permutation", R = R)
  return(data.boot)
}


#' Called by \code{testfordiff.ind.boot()}.
#' @export
#' @param data sample of compositional data
#' @param i row index
#' @param ns1 sample size of compdata.1
#' @param change.zero tolerance
#'
testfordiff.ind.boot.fun <- function(data, i, ns1, change.zero = 1e-05) {

  d <- data[i,  ]
  ns2 <- nrow(data) - ns1
  Y.1 <- d[1.:ns1,  ]
  Y.2 <- d[(ns1 + 1.):nrow(data),  ]

  alpha <- 1
  Y.1.t <- Y.1^(1/alpha)
  Y.1.t <- Y.1.t/apply(Y.1.t,1,sum)

  Y.2.t <-  Y.2^(1/alpha)
  Y.2.t <- Y.2.t/apply(Y.2.t,1,sum)

  d.mat <- alpha * create.d.mat(Y.1.t,Y.2.t)*sqrt(ncol(Y.1))


  T.chisq <- sum(d.mat)
  return(T.chisq)
}

#' Called by \code{testfordiff.ind.boot.fun()} to create a matrix of distances.
#' @export
#' @param Y.1 vector
#' @param Y.2 vector
#'
create.d.mat <- function(Y.1,Y.2) {

  ns1 <- nrow(Y.1)
  ns2 <- nrow(Y.2)
  nFA <- ncol(Y.1)

  ind.vec <-
    as.vector(unlist(tapply(seq(1,ns1,1),seq(1,ns1,1),rep,ns2)))

  Y.1.rep <- Y.1[ind.vec, ]

  Y.2.rep <- rep(t(Y.2),ns1)
  Y.2.rep <- matrix(Y.2.rep,ncol=ncol(Y.2),byrow=T)

  Y.1.split <- split(Y.1.rep,seq(1,nrow(Y.1.rep),1))
  Y.2.split <- split(Y.2.rep,seq(1,nrow(Y.2.rep),1))

  d.mat <- matrix(mapply(chisq.CA,Y.1.split,Y.2.split),byrow=T,ns1,ns2)

  return(d.mat)
}


#' Called by \code{create.d.mat()} to compute the chi-square distance.
#' @export
#' @param x1 vector
#' @param x2 vector
#'
chisq.CA <- function(x1,x2) {

  d.sq <- (x1-x2)^2
  c.vec <- x1+x2

  if ( any(d.sq!=0)) {

    d.sq[d.sq!=0] <- d.sq[d.sq!=0]/c.vec[d.sq!=0]

  }

  d.sq <- 4*sum(d.sq)

  d <- sqrt(d.sq/2)

  return(d)

}


