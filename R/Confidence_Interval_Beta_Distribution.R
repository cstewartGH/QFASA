#'
#' Returns indvidual confidence intervals and simultaneous confidence intervals
#' based (not bias corrected - see note below) on the zero-inflated beta distribution.
#' 
#' For details see:  
#'     Stewart, C. (2013)  Zero-Inflated Beta Distribution for Modeling the Proportions in 
#'     Quantitative Fatty Acid Signature Analysis.  
#'     Journal of Applied Statistics, 40(5), 985-992.
#' 
#' Note:
#' \itemize{
#'     \item These intervals are biased and should be corrected using the
#'           output from \code{\link{bias.all}}.
#'     \item \code{CI.L.1} and \code{CI.U.1} contain the simultaneous
#'           confidence intervals.
#'     \item Slow because of bisection and lots of repetition.
#'     \item Need to replace \code{p.prey} with
#'           \code{\link{p.QFASA}} eventually but just use
#'           p.prey for now. Use example where we estimate a single
#'           seal diet to compare the estimates from each method.
#' }
#' 
#' @export
#' @param seal.mat matrix containing the FA signatures of the predators.
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
#'     passed to the \code{\link{gen.pseudo.seals}} which expects fat
#'     content values for individual prey samples while
#'     \code{\link{pseaudo.seal}} and  \code{\link{p.QFASA}}
#'     expect a species average.
#' @param ext.fa subset of FA's to be used to obtain QFASA diet
#'     estimates.
#' @return indvidual confidence intervals and simultaneous confidence
#'     intervals based on the zero-inflated beta distribution. These
#'     intervals are biased and should be corrected using the output
#'     from \code{\link{bias.all}}. \code{ci.l.1} and \code{ci.u.1}
#'     contain the simultaneous confidence intervals. 
#' @references Stewart, C. (2013) Zero-inflated beta distribution for
#'     modeling the proportions in quantitative fatty acid signature
#'     analysis. Journal of Applied Statistics, 40(5), 985-992. 
#' 
beta.meths.CI <- function(seal.mat,
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
                          ext.fa) {}


#' Generate pseudo predators
#'
#' @export
#' @param prey.mat matrix containing a representative FA signature
#'     from each prey group (usually the mean). The first column must
#'     index the prey group.
#' @param diet.null.mat diet null matrix
#' @param cal.mat matrix of calibration coefficients where the \emph{i} th
#'     column is to be used with the \emph{i} th predator.
#' @param fat.cont vector of fat content
#' @param noise proportion of noise to add to the simulation.
#' @param nprey number of prey to sample.
#'
gen.pseudo.seals <- function(prey.mat,
                             diet.null.mat,
                             cal.mat,
                             fat.cont,
                             noise,
                             nprey) {}
    

#' Splits prey database into a simulation set (1/3) and a modelling
#' set (2/3). If number of samples of a prey type <=5, then prey.mod
#' AND prey.sim are duplicated instead of split.
#'
#' @export
#' @param prey.mat matrix of individual prey fatty acid signatures
#'     where the first column denotes the prey type
#' @return list of simulation prey database and modelling prey database.
#'
split.prey <- function(prey.mat) {}


#' Generates pseudo-seals. Allows:
#'
#' 1) Fat content to be included in the generated seals and
#' 2) Some species to be truly zero (that is,"zero species" do not
#' have to be included in the "noise")
#' 
#' Note:  it is assumed that sum(diet) is 1-noise
#'
#' @param prey.sim OUTPUT OF split.prey
#' @param diet DIET COMPOSITION VECTOR  (NOTE: THIS VECTOR SHOULD SUM
#'     TO 1-NOISE.  THE NOISE WILL BE ADDED TO THE diet VECTOR.) 
#' @param noise AMOUNT OF NOISE
#' @param nprey nprey TOTAL NUMBER OF PREY TO BE SAMPLED
#' @param cal CALIBRATION FACTORS
#' @param fat.cont VECTOR OF FAT CONTENT OF LENGTH=I (# OF SPECIES)
#' @param specify.noise A BOOLEAN VECTOR WITH TRUES DENOTING SPECIES
#'     TO USE IN NOISE. 
#' @return seal.star SIMULATED SEAL FA SIGNATURE.
pseudo.seal <- function(prey.sim,
                        diet,
                        noise,
                        nprey,
                        cal,
                        fat.cont,
                        specify.noise) {}



#' TODO
#'
#' @param p.mat matrix containing the FA signatures of the predators.
#' @param prey.mat matrix containing a representative FA signature
#' @param cal.mat matrix of calibration factors where the \emph{i} th
#'     column is to be used with the \emph{i} th seal. If modelling is
#'     to be done without calibration coefficients, simply pass a
#'     vector or matrix of ones. 
#' @param fat.cont prey fat content
#' @param R.bias botstrap replicates
#' @param noise noise
#' @param nprey number of prey
#' @param specify.noise noise
#' @param dist.meas distance measure
#' @param ext.fa subset of FA's to use.
#' @return row 1 is Lambda CI, row 2 is Lambda skew, and row 3 is Beta CI
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



