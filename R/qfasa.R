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
#' @docType package
#' @name QFASA
#'
NULL


#' Computes the diet estimate for each predator in seal.mat using
#' either the Kullback-Leibler Distance (KL), the Aitchison Distance
#' (AIT) or the Chi-Square Distance (CS).
#' 
#' @export
#' @param seal.mat matrix containing the FA signatures of the predators.
#' @param prey.mat matrix containing a representative FA signature
#'     from each prey group (usually the mean). The first column must
#'     index the prey group.
#' @param cal.mat matrix of calibration factors where the \emph{i} th
#'     column is to be used with the \emph{i} th seal
#' @param dist.meas distance measure to use for estimation: 1=KL,
#'     2=AIT or 3=CS
#' @param gamma parameter required for calculations using CS distance
#'     (passed to CS.obj). Currently being set to 1.
#' @param FC vector of fat content
#' @param start.val initial vector of parameters to be optimized
#' @param ext.fa subset of FA's to be used to obtain QFASA diet estimates.
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
#'  FC = preyFAs[,c(2,3)] 
#'  FC = as.vector(tapply(FC$lipid,FC$Species,mean,na.rm=TRUE))
#'
#'  ## Calibration Coefficients
#'  data(CC)
#'  cal.vec = CC[,2]
#'  cal.mat = replicate(npredators, cal.vec)
#'
#'  # Run QFASA
#'  Q = p.QFASA(predator.matrix,
#'              prey.matrix,
#'              cal.mat,
#'              dist.meas=1, 
#'              gamma=1,
#'              FC,
#'              start.val = rep(1,nrow(prey.matrix)),
#'              fa.set)
#'  
p.QFASA <- function(seal.mat,
                    prey.mat,
                    cal.mat,
                    dist.meas,
                    gamma,
                    FC = rep(1., nrow(prey.mat)),
                    start.val = rep(0.99999, nrow(prey.mat)),
                    ext.fa) {
    
    # CALIBRATING SEAL FA SIGNATURES AND THEN EXTRACTING EXTENDED DIETARY FAS
    if ((is.vector(cal.mat)) || (nrow(cal.mat) == 1.) || (ncol(cal.mat ) == 1.)) {
        
        # IF ONLY ONE SEAL
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
                                   seal = seal.mat[i,  ], 
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
                                   seal = seal.mat[i,  ],
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
                                   seal = seal.mat[i,  ],
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
#' @references Stewart, C. (2016) An approach to measure distance between compositional diet estimates containing essential
#'     zeros. Journal of Applied Statistics, 10.1080/02664763.2016.119384.
#' 
AIT.dist <- function(x.1, x.2) {
    return(sqrt(sum((log(x.1/mean.geometric(x.1)) - log(x.2/mean.geometric(x.2)))^2.)))
}


#' Used to provide additional information on various model components
#' evaluated at the optimal solution i.e. using the QFASA diet estimates and
#' Aitchison distance measure. 
#'
#' @export
#' @param alpha compositional QFASA diet estimate.
#' @param seal fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#' 
AIT.more <- function(alpha, seal, prey.quantiles) {
    
    no.zero <- sum(seal == 0.)
    seal[seal == 0.] <- 1e-05
    seal[seal > 0.] <- (1. - no.zero * 1e-05) * seal[seal > 0.]
    
    sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
    no.zero <- sum(sealhat == 0.)
    sealhat[sealhat == 0.] <- 1e-05
    sealhat[sealhat > 0.] <- (1. - no.zero * 1e-05) * sealhat[sealhat > 0.]
    
    AIT.sq.vec <- 
        ( log(seal/mean.geometric(seal)) - log(sealhat/mean.geometric(sealhat)) )^2
    
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
#' @param seal fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#' 
AIT.obj <- function(alpha, seal, prey.quantiles) {
    
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
#' @references Stewart, C. (2016) An approach to measure distance between compositional
#' diet estimates containing essential zeros.  Journal of Applied Statistics, 10.1080/02664763.2016.119384.
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
#' evaluated at the optimal solution i.e. using the QFASA diet estimates and
#' chi-square distance measure. 
#'
#' @export
#' @param alpha compositional QFASA diet estimate.
#' @param seal fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#' @param gamma power transform exponent (see chisq.dist).
#'
CS.more <- function(alpha, seal, prey.quantiles, gamma) {
    
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
#' chi-square distance measure is chosen. Unlike AIT.obj and KL.obj, does
#' not require modifying zeros.
#'
#' @export
#' @param alpha vector over which minimization takes place.
#' @param seal fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#' @param gamma power transform exponent (see chisq.dist).
CS.obj <- function(alpha, seal, prey.quantiles, gamma){

    sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
    return(chisq.dist(seal, sealhat, gamma))
}


#' Returns the distance between two compositional vectors using Kullback-Leibler
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
#' evaluated at the optimal solution i.e. using the QFASA diet estimates and
#' Kullback-Leibler distance measure. 
#'
#' @export
#' @param alpha compositional QFASA diet estimate.
#' @param seal fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#'
KL.more <- function(alpha, seal, prey.quantiles) {
    
    no.zero <- sum(seal == 0.)
    seal[seal == 0.] <- 1e-05
    seal[seal > 0.] <- (1. - no.zero * 1e-05) * seal[seal > 0.]
    
    sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
    no.zero <- sum(sealhat == 0.)
    sealhat[sealhat == 0.] <- 1e-05
    sealhat[sealhat > 0.] <- (1. - no.zero * 1e-05) * sealhat[sealhat >0.]
    
    KL.vec <- (seal - sealhat) * log(seal/sealhat)
    
    dist <- sum(KL.vec)
    
    out.list <- list(sealhat,KL.vec,KL.vec/sum(KL.vec),dist)
    names(out.list) <- c("ModFAS","DistCont", "PropDistCont","MinDist")
    
    return(out.list)
}


#' Used in \code{solnp()} as the objective function to be minimized when
#' Kullback-Leibler distance measure is chosen.
#'
#' @export
#' @param alpha vector over which minimization takes place.
#' @param seal fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#' 
KL.obj <- function(alpha, seal, prey.quantiles) {
    
    no.zero <- sum(seal == 0.)
    seal[seal == 0.] <- 1e-05
    seal[seal > 0.] <- (1. - no.zero * 1e-05) * seal[seal > 0.]
    sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
    no.zero <- sum(sealhat == 0.)
    sealhat[sealhat == 0.] <- 1e-05
    sealhat[sealhat > 0.] <- (1. - no.zero * 1e-05) * sealhat[sealhat >0.]
    return(KL.dist(seal, sealhat))
}


#' Returns the geometric mean of a compositional vector
#'
#' @param x compositional vector
#' 
mean.geometric <- function(x) {
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
    prey.means <- apply(prey.mat[, -1], 2, tapply, prey.mat[, 1], mean)
    return(prey.means)
}


#' Returns \code{sum(alpha)} and used in \code{solnp}.
#'
#' @param alpha vector over which minimization takes place.
#' @param seal fatty acid signature of predator.
#' @param prey.quantiles matrix of fatty acid signatures of
#'     prey. Each row contains an individual prey signature from a different
#'     species.
#' @param gamma power transform exponent (see chisq.dist).
#' 
QFASA.const.eqn <- function(alpha, seal, prey.quantiles, gamma) {
    return(sum(alpha))
}


#' Multiplot
#'
#' @export
#' @param ... TODO
#' @param plotlist TODO
#' @param file TODO
#' @param cols TODO
#' @param layout TODO
#' 
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {

    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols,
                         nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout),
                                                                     ncol(layout))))
        
         # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                                  layout.pos.col = matchidx$col))
        }
    }
}
