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


#' Returns the multivariate mean FA signature of each prey group
#' entered into the model. Result can be passed to prey.mat in
#' \code{qfasa()}.
#'
#' @export
#' @param prey.mat matrix containing the FA signatures of the
#'     prey. The first column indexes the prey group.
#' 
MEANmeth <- function(prey.mat) {
  prey.means <- apply(prey.mat[, -1], 2, tapply, prey.mat[, 1], mean)
  return(prey.means)
}


#' Returns the geometric mean.
#'
#' @export
#' @param x vector
#' 
mean.geometric <- function(x) {
  D <- length(x)
  return(prod(x)^(1./D))
}


#' Computes the diet estimate for each predator in seal.mat using
#' either the Kullback-Leibler Distance (KL), the Aitchison Distance
#' (AIT) or the Chi-Squared Distance (CS).
#' 
#' @export
#' @param seal.mat matrix containing the FA signatures of the predator
#' @param prey.mat matrix containing a representative FA signature
#'     from each prey group (usually the mean) -> assumes that the
#'     first column contains the name of the prey group
#' @param cal.mat matrix of calibration factors where the *i* th
#'     column is to be used with the *i* th seal
#' @param dist.measure distance measure to use for estimation: 1=KL,
#'     2=AIT or 3=CS
#' @param gamma parameter required for calculations using CS distance
#'     (passed to CS.obj). Currently being set to =1.
#' 
p.QFASA <- function(seal.mat, prey.mat, cal.mat, dist.meas, gamma, FC = rep(1., nrow(prey.mat)),
                    start.val = rep(0.99999, nrow(prey.mat)), ext.fa = ext.common.fa.list) {
  
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
  p.mat <- matrix(rep(0, nrow(prey.mat) * nrow(seal.mat)), byrow = T, nrow(seal.mat),nrow(prey.mat))
  
  more.list <- vector("list",ns)
    
  # USING SAME STARTING VECTOR FOR EACH OF THE ns SEALS
  
  if (!(is.matrix(start.val))) {
    start.val <- matrix(rep(start.val, ns), byrow = T, ns, I)
  }
  
  if (dist.meas == 1) { # KL Distance
    
    for (i in 1.:nrow(seal.mat)) {
      p.all <- Rsolnp::solnp(pars = start.val[i,  ], fun  = KL.obj, 
                     seal = seal.mat[i,  ], 
                     prey.quantiles = prey.mat, eqfun=QFASA.const.eqn, eqB=1,
                     LB = rep(0, nrow(prey.mat)),
                     UB = rep(0.999999, nrow(prey.mat)),control=list(trace=0))
      if (p.all$convergence!=0){
        
        cat("WARNING: DID NOT CONVERGE")
      }
          
      p.mat[i,  ] <- p.all$par
      
      more.list[[i]] <- KL.more(p.mat[i,],seal.mat[i,],prey.mat)
          
    }
    
  }
  else if (dist.meas==2) { # AIT distance
    for (i in 1.:nrow(seal.mat)) {
      p.all <- solnp(pars = start.val[i,  ], fun  = AIT.obj,
                     seal = seal.mat[i,  ],
                     prey.quantiles = prey.mat, eqfun=QFASA.const.eqn, eqB=1,
                     LB = rep(0, nrow(prey.mat)),
                     UB = rep(0.999999, nrow(prey.mat)),control=list(trace=0))
      if (p.all$convergence!=0){ 
        
        cat("WARNING: DID NOT CONVERGE")
      }
            
      p.mat[i,  ] <- p.all$par
      
      more.list[[i]] <- AIT.more(p.mat[i,],seal.mat[i,],prey.mat)
            
    } 
    
  }
  else { # CS distance
    for (i in 1.:nrow(seal.mat)) {
      
      p.all <- solnp(pars = start.val[i,  ], fun  = CS.obj,
                     seal = seal.mat[i,  ],
                     prey.quantiles = prey.mat, gamma = gamma, eqfun=QFASA.const.eqn, eqB=1,
                     LB = rep(0, nrow(prey.mat)),
                     UB = rep(0.999999, nrow(prey.mat)),control=list(trace=0))
      if (p.all$convergence!=0){
        
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


#' The objective function to be minimized by code{solnp()}. Similar to
#' optquantile.obj2 but does not normalize alpha.
#'
#' @export
#' @param alpha vector over which minimization takes place.
#' @param seal vector of fatty acid compositions of seal.
#' @param prey.quantiles matrix of fatty acid composition of
#'     prey. Each row contains an individual prey from a different species.
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


#' Similar to KL but requires two vectors as input.
#'
#' @export
#' @param x vector
#' @param y vector
#' @return Kulback-Liebler distance between x and y.
#' 
KL.dist <- function(x, y) {
  return(sum((x - y) * log(x/y)))
}


#' Used to provide additional information on model components when
#' alpha corresponds to the QFASA diet estimates (i.e. estimates
#' that minimized the ait distance). Used in \code{solnp()} as the
#' objective function to be minimized.
#' 
#' @param alpha vector over which minimization takes place.
#' @param seal vector of  fatty acid compositions of seal.
#' @param prey.quantiles matrix of fatty acid composition of
#'     prey. Each row contains an individual prey from a different species.
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


#' Similar to \code{optcompdiff.obj()} but does not normalize alpha.
#' Used in \code{solnp()} as the objective function to be minimized.
#'
#' @export
#' @param alpha vector over which minimization takes place.
#' @param seal vector of  fatty acid compositions of seal.
#' @param prey.quantiles matrix of fatty acid composition of
#'     prey. Each row contains an individual prey from a different
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


# Computes the difference between to vectors of compositional data
# as described in aitchison (1992) "measures of compositional difference".
#' Note that this function is different to \code{compdiff()} in that it takes square root.
#' \code{compdiff()} is used in calculating diet estimates.
#'
#' @export
#' @param x
#' param bigX
#' 
AIT.dist <- function(x, bigX) {
  return(sqrt(sum((log(x/mean.geometric(x)) - log(bigX/mean.geometric(bigX)))^2.)))
}


#' Used to provide additional information on model components when
#' alpha corresponds to the QFASA diet estimates (i.e. estimates
#' that minimized the ait distance). Used in \code{solnp()} as the
#' objective function to be minimized.
#'
#' @export
#' @param alpha vector over which minimization takes place
#' @param seal vector of  fatty acid compositions of seal
#' @param prey.quantiles matrix of fatty acid composition of
#'     prey. Each row contains an individual prey from a different
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


#' Similar to AIT.obj and KL.obj but does not require changing zeros.
#' Used in \code{solnp()} as the objective function to be minimized.
#'
#' @export
#' @param alpha vector over which minimization takes place.
#' @param seal vector of  fatty acid compositions of seal.
#' @param prey.quantiles matrix of fatty acid composition of
#'     prey. Each row contains an individual prey from a different
#'     species.
#' @param gamma parameter passed to chisq.dist
#'
CS.obj <- function(alpha, seal, prey.quantiles, gamma){

  sealhat <- t(as.matrix(alpha)) %*% prey.quantiles
  return(chisq.dist(seal,sealhat,gamma))
}


#' Calculates CS distance as described in testing for diff in diet paper.
#' This is different to chisq.ca in directory run.oct11.r.1 because
#' x1 and x2 are not assumed to be transformed.
#'
#' @export
#' @param x.1 vector
#' @param x.2 vector
#' @param alpha xxx
#' 
chisq.dist <- function(x.1,x.2,alpha) {
  
  nfa <- length(x.1)
  
  y.1 <- x.1^(alpha)
  y.1 <- y.1/sum(y.1)
  
  y.2 <- x.2^(alpha)
  y.2 <- y.2/sum(y.2)
  
  d.sq <- (y.1-y.2)^2
  c.vec <- y.1+y.2
  
  if ( any(d.sq!=0) ) {
    
    d.sq[d.sq!=0] <- d.sq[d.sq!=0]/c.vec[d.sq!=0]
  }
  
  CS.dist <- 1/alpha*sqrt(2*nfa)*sqrt(sum(d.sq))
  return(CS.dist)
}


#' Used to provide additional information on model components when
#' alpha corresponds to the QFASA diet estimates (i.e. estimates
#' that minimized the cs distance). Used in \code{solnp()} as the
#' objective function to be minimized.
#'
#' @export
#' @param alpha vector over which minimization takes place.
#' @param seal vector of  fatty acid compositions of seal.
#' @param prey.quantiles matrix of fatty acid composition of
#'     prey. Each row contains an individual prey from a different
#'     species.
#' @param gamma parameter passed to chisq.dist
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


#' QFASA.const.eqn
#'
#' @param alpha xxx
#' @param seal xxx
#' @param prey.quantiles xxx
#' @param gamma xxx
#' 
QFASA.const.eqn <- function(alpha, seal=seal.mat[i,], prey.quantiles=prey.mat, gamma=gamma) {
  
  return(sum(alpha))
}


#' Multiplot
#'
#' @export
#' @param ... xxx
#' @param plotList xxx
#' @param file xxx
#' @param cols xxx
#' @param layout xxx
#' 
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
