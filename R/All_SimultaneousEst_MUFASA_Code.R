#' Simultaneous maximum unified fatty acid signature analysis
#'
#' Returns SMUFASA calibration coefficient estimates and an average diet among
#' a sample of predators.
#'
#' @export
#' @param pred.mat matrix containing the FA signatures of the predators.
#' @param prey.mat matrix containing FA signatures from each prey group.
#'                 The first column must index the prey group.
#'                 \emph{prey.mat} is the prey database.
#' @param FC vector of fat content of length equal to the number of prey groups
#'           or species.
#' @param ext.fa subset of fatty acids to be used to obtain estimates.
#' @details Calibration coefficients (CCs) are not supplied but are
#' instead estimated. While one overall diet is computed, the CCs can
#' be used in p.QFASA or p.MUFASA to estimate individual diet
#' estimates.
#' @import dplyr compositions
#' @return A list with components:
#' \item{Cal_Estimates}{A vector of estimated calibration coefficients common
#' to all predators. The \emph{k} th value corresponds to the \emph{k} th
#' fatty acid. The estimates sum to the number of fatty acids.}
#' \item{Diet_Estimate}{A vector of estimates of the average diet among the
#' predators. The estimates are expressed as proportions summing to one.}
#' \item{Var_Epsilon}{Optimized values of error variance.}
#' \item{nll}{Negative log likelihood values.  As per \emph{solnp} documentation,
#'            \emph{nll} is "Vector of function values during optimization with
#'            last one the value at the optimal".}
#'
#'
#'
#'
#'
#'@examples
#' ## This example takes some time to run.
#' ## Please uncomment code below to run.
#'

#'#library(dplyr)
#'#library(compositions)
#'## Fatty Acids
#'#data(FAset)
#'#fa.set = as.vector(unlist(FAset))
#'
#'## Predators
#'#data(predatorFAs)
#'#tombstone.info = predatorFAs[,1:4]
#'#predator.matrix = predatorFAs[,5:(ncol(predatorFAs))]
#'#npredators = nrow(predator.matrix)
#'
#'## Prey
#'## Extracting a small number of species to speed up calculations for the example.
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
#'#out <- p.SMUFASA(predator.matrix, prey.red, FC.red, fa.set)
#'
#'#out$Cal_Estimates
#'

p.SMUFASA <- function(pred.mat, prey.mat, FC, ext.fa)
{
  npred = nrow(pred.mat)
  colnames(prey.mat)[1] <- "Species"
  prey.mat = prey.mat %>% arrange(.data$Species)
  prey.sub = prey.mat[, ext.fa]
  prey.sub = prey.sub/apply(prey.sub, 1, sum)
  group = as.vector(prey.mat$Species)
  prey.mat = cbind(group, prey.sub)
  prey.mean = MEANmeth(prey.mat)
  pred.mat = pred.mat[, ext.fa]
  pred.mat = pred.mat/apply(pred.mat, 1, sum)
  if (any(prey.sub == 0)) {
    prey.mat.t <- mod.zeros.FA.sig.mat(prey.sub, 5e-05)
    prey.mat.t <- compositions::ilr(prey.mat.t)
  }
  else {
    prey.mat.t <- compositions::ilr(prey.sub)
  }
  prey.mat.t <- data.frame(group = group, mat = prey.mat.t)
  prey.mean.t = MEANmeth(prey.mat.t)
  rownames(prey.mean.t) <- prey.mean.t[, 1]
  prey.var.t = POOLVARmeth(prey.mat.t)$Spool
  Q = p.QFASA(pred.mat, prey.mean, rep(1,length(ext.fa)), dist.meas = 2,
              gamma = 1, FC, start.val = rep(1, nrow(prey.mean)), ext.fa)
  pred <- mod.zeros.FA.sig.mat(pred.mat, 5e-05)

  pred.mat.t <- compositions::ilr(pred)

  ers <- matrix(NA, nrow = 100 * npred, ncol = (ncol(prey.mat) - 1))
  for (j in 1:100) {
    yest <- matrix(NA, nrow = npred, ncol = (ncol(prey.mat) - 1))
    for (i in 1:nrow(yest)) {
      yest[i, ] <- pseudo.pred.norm(prey.mean.t, prey.var.t,
                                    Q$`Diet Estimates`[i, ])
    }

    pred.mat <- compositions::acomp(pred.mat)

    yest <- compositions::acomp(yest)

    lb <- (j - 1) * npred + 1
    ub <- j * npred
    ers[lb:ub, ] = -yest + pred.mat
  }

  ers <- compositions::acomp(ers)

  V <- compositions::ilrBase(D = length(ext.fa))

  G <- V %*% t(V)
  sep.start <- diag(-1/2 * t(V) %*% G %*% compositions::variation(ers) %*%
                      G %*% V)


  quan <- stats::quantile(sep.start)
  quan.start <- numeric(4)
  groupind <- numeric(length(sep.start))
    for (j in 1:4) {
    quan.start[j] <- mean(sep.start[sep.start >= quan[j] &
                                      sep.start <= quan[j + 1]])
    groupind[sep.start >= quan[j] & sep.start <= quan[j +  1]] <- j
  }
  spec <- unique(prey.mat$group)
  spec <- spec[order(spec)]
  I <- length(spec)
  n <- tapply(prey.mat$group, prey.mat$group, length)
  K <- length(ext.fa)-1


  my.start <- apply(Q$`Diet Estimates`[, -I],2,mean)
  #parameters <- list(alpha = my.start,
  #                   cal= rep(1,K),
  #                   z  = array(rep(prey.mean.t),
  #                              c(nrow(prey.mean.t), ncol(prey.mean.t),
  #                                npred)),
  #                   sepsilon = quan.start)

  parameters <- list(alpha = my.start,
                     cal= rep(1/(K+1),K),
                     z  = array(rep(prey.mean.t),
                                c(nrow(prey.mean.t), ncol(prey.mean.t),
                                  npred)),
                     sepsilon = quan.start)


  data <- list(y = pred.mat.t, n=n, varz=prey.var.t, mu = prey.mean.t,
               V=V, sind=groupind)

  objnt <- TMB::MakeADFun(data,parameters,random="z",DLL="CommonDiet")



  LB = c(rep(0,(I-1)), rep(0,K),rep(0,length(quan.start)))
  UB = c(rep(1,(I-1)), rep(1,K),rep(Inf,length(quan.start)))


  al.cal.sum <- function(pars){
    npars <- length(pars)
    al <- pars[1:(npars-length(quan.start)-K)]
    cal <- pars[I:(npars-length(quan.start))]
    return(c(sum(al),sum(cal)))
  }

  optnt <- Rsolnp::solnp(pars=objnt$par, fun=objnt$fn, ineqfun=al.cal.sum,
                         ineqLB = c(0,0), ineqUB = c(1,1), LB=LB, UB = UB,
                         control=list(delta=1e-04, tol=1e-04))

  npars <- length(objnt$par)
  L <- optnt$values
  alpha <- optnt$pars[1:(I-1)]
  alpha  <-  c(alpha,(1 - sum(alpha)))
  cal.est <- optnt$pars[I:(npars-length(quan.start))]
  #cal.est <- c(cal.est,(K+1 - sum(cal.est)))

  cal.est <- c(cal.est,(1 - sum(cal.est)))*(K+1)


  seps <- optnt$pars[(npars-length(quan.start)+1):length(optnt$pars)]
  alpha <- alpha/FC
  alpha <- alpha/sum(alpha)
  names(cal.est)  <- ext.fa
  names(alpha) <- colnames(Q$`Diet Estimates`)
  return(list(Cal_Estimates = cal.est, Diet_Estimates = alpha, nll = L, Var_Epsilon = seps))
}




