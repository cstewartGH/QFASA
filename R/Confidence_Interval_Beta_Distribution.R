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
                          ext.fa) {
    
    seal.mat <- predator.mat
}



#' Splits prey database into a simulation set (1/3) and a modelling
#' set (2/3). If number of samples of a prey type is less than or equal to 5,
#' then \code{prey.mod} and \code{prey.sim} are duplicated instead of split.
#'
#' @export
#' @param prey.mat matrix of individual prey fatty acid signatures
#'     where the first column denotes the prey type
#' @return list of simulation prey database and modelling prey database.
#'
split.prey <- function(prey.mat) {}


#' Calculate bias correction for confidence intervals from \code{\link{beta.meths.CI()}}.
#'
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
                     ext.fa = ext.common.fa.list) {}


