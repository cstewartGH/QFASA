#' Computes within species variance-covariance matrices on transformed scaled,
#' along wtih a pooled estimate.
#'
#' @param prey.mat matrix containing transformed FA signatures of the prey.
#'                 Note that the first column indexes prey type.
#'
#' @return Returns the variance-covariance matrix of each prey type as well as
#'         a pooled estimate of the variance-covariance matrix
#' @import dplyr
#' @export
#'
POOLVARmeth <- function(prey.mat) {

  groups <- unique(prey.mat[,1])
 colnames(prey.mat)[1] <- "group"
 n <- tapply(prey.mat[,1], prey.mat[,1], length)
 prey.var <- vector("list", length(groups))
 Spool <- matrix(0, nrow=(ncol(prey.mat)-1), ncol=(ncol(prey.mat)-1))

   for(i in 1:length(groups)){
    prey.sub <- prey.mat %>% filter(.data$group==groups[i])
    prey.var[[i]] <- stats::var(prey.sub[,-1])
    Spool <- Spool + (n[[i]]-1)*prey.var[[i]]
   }
 Spool <- Spool/(sum(n) - length(n))

  return(list(prey.var=prey.var, Spool=Spool))
}

#' Generate a pseudo predator parametrically from multivariate normal
#' distributions.
#'
#' @param mu.mat matrix where each row represents the mean transformed FA
#'               signature of each prey type
#' @param sigma.pool pooled variance-covariance matrix of the transformed
#'                   fatty acid signatures of prey types
#' @param diet vector of proportions of prey species in diet (true diet)
#' @details Similar to \emph{pseudo.pred} but instead generates the
#' pseudo-predators parametrically by assuming ilr transformed FA signatures
#' have a multivariate normal distribution.
#'
#' @return A simulated predator FA signature.  See \emph{pseudo.pred} for an
#' example illustrating how to generate a sample of pseudo predators.
#' @export
#'
pseudo.pred.norm <- function(mu.mat, sigma.pool, diet){

  J <- length(diet)

  x.mat <- matrix(0, nrow=J, ncol=ncol(mu.mat))
  for (j in 1:J){
    x.mat[j,] <- MASS::mvrnorm(1, mu=mu.mat[j,], Sigma=sigma.pool)
  }

  x.mat.o <- compositions::ilrInv(x.mat)
  x.mat.o <- as.data.frame(x.mat.o)
  x.mat.o <- as.matrix(x.mat.o)
  Y <- diet %*% x.mat.o
  return(Y)

}


#' Modifies the zeros for a sample of FA signatures using the multiplicative
#' replacement method (and same delta for every zero).
#'
#' @keywords internal
#'
mod.zeros.FA.sig.mat <- function(y.mat,delta) {

   y.mat <- t(apply(y.mat,1,mod.zeros.FA.sig,delta=delta))

  return(y.mat)
}

#' Modifies the zeros for a single FA signature using the multiplicative
#' replacement method (and same delta for every zero).
#'
#' @keywords internal
#'
mod.zeros.FA.sig <- function(y,delta) {

  # MODIFIES THE ZEROS IN A SINGLE FA SIGNATURES USING THE MULTIPLICATIVE REPLACEMENT
  # STRATEGY

  no.zero <- sum(y==0)
  y[y>0] <- (1 - no.zero*delta)*y[y>0]
  y[y == 0] <- delta

  return(y)
}

#' Returns MUFASA diet estimates corresponding to a sample of predators.
#'
#' Computes the diet estimate for each predator in \emph{pred.mat} using MLE
#' method.
#'
#' @export
#' @param pred.mat matrix containing the FA signatures of the predators.
#' @param prey.mat matrix containing a representative FA signature
#'     from each prey group (usually the mean). The first column must
#'     index the prey group or species.
#' @param cal.mat matrix of calibration factors where the \emph{i} th
#'     column is to be used with the \emph{i} th predator. If modelling is to be
#'     done without calibration coefficients, simply pass a vector or matrix of
#'     ones.
#' @param FC vector of fat content of length equal to the number of prey groups
#'           or species.
#' @param ext.fa subset of fatty acids to be used to obtain QFASA diet estimates.
#' @import dplyr compositions
#' @return A list with components:
#' \item{Diet_Estimates}{A matrix of the diet estimates for each predator where
#' each row corresponds to a predator and the columns to prey species.
#' The estimates are expressed as proportions summing to one.}
#' \item{nll}{Negative log likelihood values.  As per \emph{solnp} documentation,
#'            \emph{nll} is "Vector of function values during optimization with
#'            last one the value at the optimal".}
#' \item{Var_Epsilon}{Optimized values of error variance.  See reference.}
#'
#'
#' @references Steeves, Holly (2020) Maximum likelihood approach to diet
#' estimation and inference based on fatty acid signatures. PhD thesis available
#' at https://dalspace.library.dal.ca/handle/10222/80034.
#'@examples
#'
#'  ##  This example takes some time to run.
#'  ## Please uncomment code below to run.
#'
#'#library(dplyr)
#'#library(compositions)
#'  ## Fatty Acids
#'#data(FAset)
#'#fa.set = as.vector(unlist(FAset))
#'
#' ## Predators
#'#data(predatorFAs)
#'#tombstone.info = predatorFAs[,1:4]
#'#predator.matrix = predatorFAs[,5:(ncol(predatorFAs))]
#'#npredators = nrow(predator.matrix)
#'
#' ## Prey
#' ## Extracting a small number of species to speed up calculations for the example.
#'#data(preyFAs)
#'#prey.matrix = preyFAs[,-c(1,3)]
#'#spec.red <-c("capelin", "herring", "mackerel", "pilchard", "sandlance")
#'#spec.red <- sort(spec.red)
#'#prey.red <- prey.matrix %>% filter(Species %in% spec.red)
#'
#'## Fat content
#'#FC = preyFAs[,c(2,3)]
#'#FC = FC %>%  arrange(Species)
#'#FC.vec = tapply(FC$lipid,FC$Species,mean,na.rm=TRUE)
#'#FC.red <- FC.vec[spec.red]
#'
#'## Calibration Coefficients
#'#data(CC)
#'#cal.vec = CC[,2]
#'#cal.m = replicate(npredators, cal.vec)
#'#rownames(cal.m) <- CC$FA
#'
#'#M <- p.MUFASA(predator.matrix, prey.red, cal.m, FC.red, fa.set)
#'
#'## Diet EStimates
#'
#'#M$Diet_Estimates
#'

p.MUFASA <- function(pred.mat,
                        prey.mat,
                        cal.mat,
                        FC,
                        ext.fa) {

  npred = nrow(pred.mat)

  colnames(prey.mat)[1] <- "Species"

  # Removing extended dietary FAs
  prey.mat= prey.mat %>% arrange(.data$Species)
  prey.sub=prey.mat[,ext.fa]
  prey.sub=prey.sub/apply(prey.sub,1,sum)
  group=as.vector(prey.mat$Species)
  prey.mat=cbind(group,prey.sub)

  prey.mean=MEANmeth(prey.mat)

  pred.mat = pred.mat[,ext.fa]
  pred.mat=pred.mat/apply(pred.mat,1,sum)

  cal.mat <- cal.mat[ext.fa,]

  ## Transforming Prey
  if(any(prey.sub==0)) {
    prey.mat.t <- mod.zeros.FA.sig.mat(prey.sub,0.00005) #analytic precision is two decimals after the percentage.
    prey.mat.t <- compositions::ilr(prey.mat.t)
  } else{
    prey.mat.t <- compositions::ilr(prey.sub)
  }
  prey.mat.t <- data.frame(group=group,mat=prey.mat.t)

  prey.mean.t=MEANmeth(prey.mat.t)
  rownames(prey.mean.t) <- prey.mean.t[,1]
  prey.var.t=POOLVARmeth(prey.mat.t)$Spool

  # Running original QFASA
  Q = p.QFASA(pred.mat,
              prey.mean,
              cal.mat,
              dist.meas = 1,
              gamma=1,
              FC,
              start.val = rep(1,nrow(prey.mean)),
              ext.fa)

  # Preds after calibration coefficients

  if ((is.vector(cal.mat)) || (nrow(cal.mat) == 1.) || (ncol(cal.mat ) == 1.)) {

    ## IF ONLY ONE SEAL
    pred.mat <- t(t(pred.mat)/as.vector(unlist(cal.mat)))
    pred.mat <- pred.mat/apply(pred.mat, 1., sum)
    pred.mat <- as.matrix(pred.mat)
  }
  else {
    pred.mat <- pred.mat/t(cal.mat)
    pred.mat <- pred.mat/apply(pred.mat, 1., sum)
    pred.mat <- as.matrix(pred.mat)
  }

  pred <- mod.zeros.FA.sig.mat(pred.mat,0.00005)
  pred.mat.t <- compositions::ilr(pred)

  ## Estimating start values for the error
  ers <- matrix(NA, nrow = 50*npred, ncol = (ncol(prey.mat)-1))

  # Creating 50 pseudo-predators based on  QFASA diet estimates and prey info
  for(j in 1:50){
    yest <- matrix(NA, nrow = npred, ncol = (ncol(prey.mat)-1))
    for(i in 1:nrow(yest)){
      yest[i,] <- pseudo.pred.norm(prey.mean.t, prey.var.t, Q$`Diet Estimates`[i,])
    }

    pred.mat <- compositions::acomp(pred.mat)

    yest <- compositions::acomp(yest) #without error, untransformed
    lb <- (j-1)*npred + 1
    ub <- j*npred
    ers[lb:ub,] = -yest + pred.mat
  }

  ers <- compositions::acomp(ers)
  V <-compositions::ilrBase(D=length(ext.fa))
  G <- V%*%t(V)
  sep.start <- diag(-1/2*t(V)%*%G%*%compositions::variation(ers)%*%G%*%V)

  # Getting the quantiles of the diagonal
  quan <- stats::quantile(sep.start)
  quan.start <- numeric(4)
  groupind <- numeric(length(sep.start))
  for(j in 1:4){
    quan.start[j] <- mean(sep.start[sep.start>=quan[j] & sep.start<=quan[j+1]])
    groupind[sep.start>=quan[j] & sep.start<=quan[j+1]] <- j
  }

  spec <- unique(prey.mat$group)
  spec <- spec[order(spec)]
  I <- length(spec)

  # getting the number sampled from each species
  n <- tapply(prey.mat$group, prey.mat$group, length)

  parameters <- list(alpha = Q$`Diet Estimates`[,-I],
                     z  = array(rep(prey.mean.t), c(nrow(prey.mean.t), ncol(prey.mean.t), npred)),
                     sepsilon = quan.start
  )

  #Data to send to tmb
  data <- list(y = pred.mat.t,
               n=n,
               varz=prey.var.t,
               mu = prey.mean.t,
               V=V,
               sind=groupind
  )

  objnt <- TMB::MakeADFun(data,parameters,random="z",DLL="QFASA")
  npars <- length(objnt$par)


  # constraining alphas to be between 0 and 1
  # Full parameter list
  lb <- rep(0,npars)
  ub <- c(rep(1,(npars-length(quan.start))), rep(Inf,length(quan.start)))


  # Full parameter set
  al.sum <- function(pars){
    npars <- length(pars)
    alpha <- matrix(pars[1:(npars-length(quan.start))], ncol=I-1)
    return(apply(alpha,1,sum))
  }

  # Full parameter set
  optnt <- Rsolnp::solnp(pars=objnt$par, fun=objnt$fn, ineqfun=al.sum, ineqLB = rep(0,npred), ineqUB = rep(1,npred), LB=lb, UB=ub, control=list(delta=0.0001, tol=0.00001))

  L <- optnt$values

  alpha <- matrix(optnt$pars[1:((I-1)*(npred))],nrow=npred)
  seps <- optnt$pars[((I-1)*(npred)+1):length(optnt$pars)]
  alpha <- cbind(alpha,1-apply(alpha,1,sum))

  if (is.matrix(FC)) {
    FC.mat <- FC
  }
  else {
    FC.mat <- matrix(rep(FC, npred), byrow = T, nrow=npred, ncol=I)
  }

  alpha <- alpha/FC.mat
  alpha <- alpha/rowSums(alpha)
  colnames(alpha) <- colnames(Q$`Diet Estimates`)

  return(list(Diet_Estimates = alpha, nll = L, Var_Epsilon = seps))
}








