#' Repeatability in Diet Estimates
#'
#' Computes a measure of repeatability for a sample of predators with repeated
#' diet estimate measurements.
#'
#' @param data data frame of diet estimates. First column must denote the predator
#'       and second column the second factor (eg. year or season).
#' @param prey.database prey data base that was used to compute the QFASA diet
#'       estimates in \code{data}. Will be used to generate pseudo predators.
#' @param fatcont.mat data frame or matrix of length equal to the number of prey
#'      FA signatures in prey data base.  First column is name of species and
#'      second column is lipid.
#' @param dist.meas distance measure to use in \code{p.QFASA}.
#' @param ext.fa subset of FAs to use.
#' @param B number of pseudo predators samples to generate for bias calculation.
#'       Default is set to 50 because is slow to run.
#' @param R number of bootstrap samples (i.e. \code{R} samples for each generated sample
#'       of pseudo predators). Default is set to 100 because it is slow to run.
#' @param CI indicates if a confidence interval for rho is to be calculated.
#'      Default is FALSE since this is time consuming to obtain.
#' @param alpha a (1-alpha/2)X100 percent confidence interval is calculated for rho
#'      if \code{CI=TRUE}.
#' @param gamma.QFASA if \code{dist.meas=3}, gamma is required.  Default is 1.
#' @param gamma.rho value of gamma to be used to compute CS distance
#'       in repeatablity functions.  Default is 1.
#' @return Bias corrected measure of repeatability, estimate of the bias and
#'      (if \code{CI=TRUE}) a confidence interval for the true repeatability.
#' @references "Repeatability for Compositional Diet Estimates with Zeros".
#' Contact Connie Stewart (cstewart@unb.ca).
#' @export
#'
#' @examples
#'
#'  ##  These examples take some time to run.
#'  ##  Please uncomment code below to run them.
#'
#' # data(preyFAs)
#' # data(FAset)
#'
#' ## Balanced Diet Data
#'
#' #my.preybase <- preyFAs[,-c(1,3)]
#' #my.preybase[,-1] <- my.preybase[,-1]/rowSums(my.preybase[,-1])
#'
#' #set.seed(10)
#'
#' #comp.rep(data = bal.diet.data,prey.database=my.preybase,
#' #fatcont.mat = as.data.frame(preyFAs[,c(2,3)]),dist.meas=2,
#' #ext.fa = as.vector(unlist(FAset)))
#'
#' ## Unbalanced Diet Data
#'
#' # my.preybase <- preyFAs[,-c(1,3)]
#' # my.preybase[,-1] <- my.preybase[,-1]/rowSums(my.preybase[,-1])
#'
#' # set.seed(10)
#'
#' # comp.rep(unbal.diet.data,my.preybase,as.data.frame(preyFAs[,c(2,3)]),2,
#' # as.vector(unlist(FAset)))

comp.rep <- function(data, prey.database, fatcont.mat, dist.meas, ext.fa, B=50,
                     R=100, CI=FALSE, alpha=0.05, gamma.QFASA = 1,
                     gamma.rho = 1) {

  # BALANCED CASE

  # DETERMINE IF DESIGN IS BALANCED

  k.vec <- as.vector(table(factor(data[,2])))

  if ( any(k.vec!=k.vec[1]) ) {

    # UNBALANCED DESIGN

    cat("Unbalanced Case:",fill=TRUE)

    out <- unbal.rep.CI(data,B,R,CI,alpha,prey.database,fatcont.mat,dist.meas,ext.fa)

  } else {

    # BALANCED DESIGN

    cat("Balanced Case:", fill=TRUE)

    out <- bal.rep.CI(data,B,R,CI,alpha,prey.database,fatcont.mat,dist.meas,ext.fa)
  }

  return(out)
}

#' Measure of Repeatability for Diet Estimates (Balanced Case)
#'
#' Computes a measure of repeatability for a sample of predators with repeated
#' diet estimate measurements (Balanced Case)
#'
#' @param data data frame of diet estimates. First column must denote the predator
#'        and second column the second factor (eg. year or season).  Must contain
#'        the same number of repeated measurements for each predator.
#' @param B number of pseudo predators samples to generate for bias calculation.
#'       Default is set to 50 because is slow to run.
#' @param R number of bootstrap samples (i.e. \code{R} samples for each generated sample
#'       of pseudo predators). Default is set to 100 because it is slow to run.
#' @param CI indicates if a confidence interval for rho is to be calculated.
#'      Default is FALSE since this is time consuming to obtain.
#' @param alpha a (1-alpha/2)X100 percent confidence interval is calculated for rho
#'      if \code{CI=TRUE}.
#' @param prey.database prey data base that was used to compute the QFASA diet
#'       estimates in \code{data}. Will be used to generate pseudo predators.
#' @param fatcont.mat data frame or matrix of length equal to the number of prey
#'      FA signatures in prey data base.  First column is name of species and
#'      second column is lipid.
#' @param dist.meas distance measure to use in \code{p.QFASA}.
#' @param ext.fa subset of FAs to use.
#' @param gamma.QFASA if \code{dist.meas=3}, gamma is required.  Default is 1.
#' @param gamma.rho value of gamma to be used to compute CS distance
#'       in repeatablity functions.  Default is 1.
#' @return Bias corrected measure of repeatability, estimate of the bias and
#'      (if \code{CI=TRUE}) a confidence interval for the true repeatability.
#' @keywords internal

bal.rep.CI <- function(data, B, R, CI=TRUE, alpha, prey.database, fatcont.mat,
                       dist.meas, ext.fa, gamma.QFASA = 1, gamma.rho = 1) {

  # Similar to two.factor.rep.CI but allows only the bias to be calculated if CI=FALSE.
  # Also, only returning T intervals since they performed well in simulations
  # and were similar to BCa intervals.

  seal.ids <- unique(data[,1])
  ns <- length(seal.ids)

  rho.out <- two.factor.rep(data,gamma.rho)
  rho.t0 <- rho.out[[2]]
  k <- rho.out[[3]]

  #print(rho.t0)

  rho.mat.star<- matrix(rep(NA,B*R),nrow = B,ncol = R, byrow=T)
  rho.mat.star.se <- matrix(rep(NA,B*R),nrow = B,ncol = R, byrow=T)

  # TO ESTIMATE THE BIAS
  rho.b <- rep(NA,B)

  for (b in 1:B) {

    p.mat.b <- pseudo.pred.table(data, prey.database, fatcont.mat, dist.meas,
                                 ext.fa, gamma.QFASA)

    rho.b[b] <- two.factor.rep(p.mat.b,gamma.rho)[[2]]

    if (CI==TRUE) {
      for (r in 1:R) {
        ind <- sample(seal.ids,replace=TRUE)
        rho.star.out <- rho.boot.fun(ind, p.mat.b, gamma.rho)
        rho.mat.star[b,r] <- rho.star.out[[1]]
        rho.mat.star.se[b,r] <- bootstrap::jackknife(1:ns, rho.jack.fun,
                                                     data = rho.star.out[[2]],
                                                     gamma.rho = gamma.rho)$jack.se
      }
    } # End if

  } # End b loop


  if (CI==TRUE) {
    t.star.1 <- as.vector( rho.mat.star - matrix(rep(rho.b,R),nrow=B) )
    t.star.2 <- as.vector( (rho.mat.star -
                              matrix(rep(rho.b,R),nrow=B))/rho.mat.star.se )

    # REMOVE OUTLIERS

    #cat("outliers")
    t.out <- graphics::boxplot(t.star.2,plot=FALSE)$out
    t.star.2 <- t.star.2[!t.star.2 %in% t.out ]

    q.t.star.1 <- stats::quantile(t.star.1,c(alpha/2,1-alpha/2), type=1)
    q.t.star.2 <- stats::quantile(t.star.2,c(alpha/2,1-alpha/2), type=1)

    #cat("BASIC BOOTSTRAP CI")

    #NEW
    #CI.L.boot <- 2*mean(rho.b) - q.t.star.1[2]
    #CI.U.boot <- 2*mean(rho.b) - q.t.star.1[1]
    #print(c(CI.L.boot,CI.U.boot))

    cat("BOOTSTRAP T CI")

    #NEW
    CI.L.boot.t <- mean(rho.b) - q.t.star.2[2]*stats::sd(rho.mat.star)
    CI.U.boot.t <- mean(rho.b) - q.t.star.2[1]*stats::sd(rho.mat.star)

    print(c(CI.L.boot.t,CI.U.boot.t))

    #cat("BCa CI")
    # FROM STATISTICAL COMPUTING WITH R BY MARIA L. RIZZO - 2008

    # NEW

    #if ( sum( as.vector( rho.mat.star < matrix(rep(rho.b,R),nrow=B) ) )/(R*B)==1 ) {
    #
    #    z0 <- qnorm(0.99999999)
    #
    #  } else {
    #
    #    z0 <- qnorm( sum( as.vector( rho.mat.star < matrix(rep(rho.b,R),nrow=B) ))/(R*B))
    #
    #  }

    #jack.val <- bootstrap::jackknife(seal.ids, rho.jack.fun, data = data, gamma.rho = gamma.rho)$jack.values

    #L <- mean(jack.val) - jack.val
    #a <- sum(L^3)/( 6* sum(L^2)^1.5 )

    #alpha.vec <- c(alpha/2,1-alpha/2)
    #zalpha <- qnorm(alpha.vec)

    #adj.alpha <- pnorm(z0 +(z0 + zalpha)/(1-a*(z0+zalpha)))

    #q.BCa <- quantile(as.vector(rho.mat.star),adj.alpha,type=6)
    #CI.L.BCa <- q.BCa[1]
    #CI.U.BCa <- q.BCa[2]

    #print(c(CI.L.BCa,CI.U.BCa))

  }  # End if

  # Bias

  bias <- mean(rho.b) - rho.t0

  # Adjusted rho.hat

  rho.hat.adj <- rho.t0-bias

  if (CI==TRUE) {

    CI.L.adj <- CI.L.boot.t - 2*bias
    CI.U.adj <- CI.U.boot.t - 2*bias
    out.list <- list(rho.hat.adj,c(CI.L.adj,CI.U.adj),bias)
    names(out.list) <- c("Bias-corrected estimate of repeatability",
                         "Bias-corrected bootstrap T interval","Estimated Bias")
    return(out.list)

  }  else {

    out.list <- c(rho.hat.adj,bias)
    names(out.list) <- c("Bias-corrected estimate of repeatability","Estimated bias")
    return(out.list)
  }

}

  #' Measure of Repeatability for Diet Estimates (Unbalanced/Missing Value Case)
  #'
  #' Computes a measure of repeatability for a sample of predators with repeated
  #' diet estimate measurements (Unbalanced/Missing Value Case)
  #'
  #' @param data data frame of diet estimates. First column must denote the predator
  #'        and second column the second factor (eg. year or season).  Need not contain
  #'        the same number of repeated measurements for each predator.
  #' @param B number of pseudo predators samples to generate for bias calculation.
  #'       Default is set to 50 because is slow to run.
  #' @param R number of bootstrap samples (i.e. \code{R} samples for each generated sample
  #'       of pseudo predators). Default is set to 100 because it is slow to run.
  #' @param CI indicates if a confidence interval for rho is to be calculated.
  #'      Default is FALSE since this is time consuming to obtain.
  #' @param alpha a (1-alpha/2)X100 percent confidence interval is calculated for rho
  #'      if \code{CI=TRUE}.
  #' @param prey.database prey data base that was used to compute the QFASA diet
  #'       estimates in \code{data}. Will be used to generate pseudo predators.
  #' @param fatcont.mat data frame or matrix of length equal to the number of prey
  #'      FA signatures in prey data base.  First column is name of species and
  #'      second column is lipid.
  #' @param dist.meas distance measure to use in \code{p.QFASA}.
  #' @param ext.fa subset of FAs to use.
  #' @param gamma.QFASA if \code{dist.meas=3}, gamma is required.  Default is 1.
  #' @param gamma.rho value of gamma to be used to compute CS distance
  #'       in repeatablity functions.  Default is 1.
  #' @return Bias corrected measure of repeatability, estimate of the bias and
  #'      (if \code{CI=TRUE}) a confidence interval for the true repeatability.
  #' @keywords internal

  unbal.rep.CI <- function(data, B, R, CI=TRUE, alpha, prey.database, fatcont.mat,
                         dist.meas, ext.fa, gamma.QFASA = 1, gamma.rho = 1) {

    # Similar to unbal.two.factor.rep.CI but allows only the bias to be calculated
    # if CI=FALSE. Also, only returning T intervals since they performed well
    # in simulations and were similar to BCa intervals.


    seal.ids <- unique(data[,1])
    ns <- length(seal.ids)

    rho.out <- unbal.two.factor.rep(data,1,gamma.rho)
    rho.t0 <- rho.out[[3]]

    rho.mat.star<- matrix(rep(NA,B*R),nrow = B,ncol = R, byrow=T)
    rho.mat.star.se <- matrix(rep(NA,B*R),nrow = B,ncol = R, byrow=T)

    # TO ESTIMATE THE BIAS
    rho.b <- rep(NA,B)

    for (b in 1:B) {

      p.mat.b <- pseudo.pred.table(data, prey.database, fatcont.mat, dist.meas, ext.fa, gamma.QFASA)
      rho.b[b] <- unbal.two.factor.rep(p.mat.b,1,gamma.rho)[[3]]
      #cat("rho.b")
      #print(rho.b[b])

      if (CI==TRUE) {

        for (r in 1:R) {

          ind <- sample(seal.ids,replace=TRUE)
          rho.star.out <- unbal.rho.boot.fun(ind, p.mat.b, gamma.rho)
          rho.mat.star[b,r] <- rho.star.out[[1]]
          rho.mat.star.se[b,r] <- bootstrap::jackknife(1:ns, unbal.rho.jack.fun, data = rho.star.out[[2]], gamma.rho = gamma.rho)$jack.se
        }
      } # End if
    } # End b loop

    # NEW

    if (CI==TRUE) {

      t.star.1 <- as.vector( rho.mat.star - matrix(rep(rho.b,R),nrow=B) )
      t.star.2 <- as.vector( (rho.mat.star - matrix(rep(rho.b,R),nrow=B))/rho.mat.star.se )

      # REMOVE OUTLIERS

      t.out <- graphics::boxplot(t.star.2,plot=FALSE)$out
      t.star.2 <- t.star.2[!t.star.2 %in% t.out ]

      q.t.star.1 <- stats::quantile(t.star.1,c(alpha/2,1-alpha/2), type=1)
      q.t.star.2 <- stats::quantile(t.star.2,c(alpha/2,1-alpha/2), type=1)


      #cat("BASIC BOOTSTRAP CI")

      ## NEW

      #CI.L.boot <- 2*mean(rho.b) - q.t.star.1[2]
      #CI.U.boot <- 2*mean(rho.b) - q.t.star.1[1]

      #print(c(CI.L.boot,CI.U.boot))

      cat("BOOTSTRAP T CI")

      #NEW

      CI.L.boot.t <- mean(rho.b) - q.t.star.2[2]*stats::sd(rho.mat.star)
      CI.U.boot.t <- mean(rho.b) - q.t.star.2[1]*stats::sd(rho.mat.star)


      print(c(CI.L.boot.t,CI.U.boot.t))

      #cat("BCa CI")
      ## FROM STATISTICAL COMPUTING WITH R BY MARIA L. RIZZO - 2008

      # NEW

      #if ( sum( as.vector( rho.mat.star < matrix(rep(rho.b,R),nrow=B) ) )/(R*B)==1 ) {

      #  z0 <- qnorm(0.99999999)
      #
      #} else {
      #
      #  z0 <- qnorm( sum( as.vector( rho.mat.star < matrix(rep(rho.b,R),nrow=B) ))/(R*B))
      #
      #}

      #jack.val <- bootstrap::jackknife(seal.ids, unbal.rho.jack.fun, data = data, gamma.rho = gamma.rho)$jack.values
      #
      #L <- mean(jack.val) - jack.val
      #a <- sum(L^3)/( 6* sum(L^2)^1.5 )

      #alpha.vec <- c(alpha/2,1-alpha/2)
      #zalpha <- qnorm(alpha.vec)

      #adj.alpha <- pnorm(z0 +(z0 + zalpha)/(1-a*(z0+zalpha)))
      #print(adj.alpha)

      #q.BCa <- quantile(as.vector(rho.mat.star),adj.alpha,type=6)
      #CI.L.BCa <- q.BCa[1]
      #CI.U.BCa <- q.BCa[2]

      #print(c(CI.L.BCa,CI.U.BCa))

    }



    # Bias

    bias <- mean(rho.b) - rho.t0

    # Adjusted rho.hat

    rho.hat.adj <- rho.t0-bias

    if (CI==TRUE) {

      CI.L.adj <- CI.L.boot.t - 2*bias
      CI.U.adj <- CI.U.boot.t - 2*bias
      out.list <- list(rho.hat.adj,c(CI.L.adj,CI.U.adj),bias)
      names(out.list) <- c("Bias-corrected estimate of repeatability", "Bias-corrected bootstrap T interval","Estimated Bias")
      return(out.list)

    }  else {

      out.list <- c(rho.hat.adj,bias)
      names(out.list) <- c("Bias-corrected estimate of repeatability","Estimated bias")
      return(out.list)
    }

  }

  #' Measure of (unadjusted) repeatability
  #'
  #' @param dataset data frame of diet estimates. Columns 1 and 2 of
  #'        \code{dataset} must be the first and second factors in your two factor
  #'        ANOVA model where column 1 specifies the predator and column 2 the
  #'        year, season, etc...
  #' @param gamma needed to compute CS distance.  Default is 1.
  #'
  #' @return unadjusted consistent repeatability, unadjusted absolute repeatability
  #'         and the number of levels of second factor (eg. year, seasons, etc..)
  #' @details We use absolute repeatability (or intraclass correlation coefficient).
  #'          The difference in the two types is discussed in:
  #'          McGraw and Wong (1996) Forming Inferences about some intraclass
  #'          correlation coefficients.  Psychological Methods. 1(1):30-46.
  #' @keywords internal

  two.factor.rep<- function(dataset, gamma=1){

    # CALCULATES CONSISTENT AND ABSOLUTE REPEATABILITY FROM A MATRIX OF PREDATOR
    # DIETS. THIS FUNCTION ASSUMES A BALANCED DESIGN.

    # COLUMNS 1 AND 2 OF DATASET MUST BE THE FIRST AND SECOND FACTORS IN YOUR TWO FACTOR
    # ANOVA MODEL WHERE COLUMN 1 SPECIFIES THE PREDATOR AND COLUMN 2 THE YEAR, SEASON, ETC...

    # THE FUNCTION WILL RETURN TWO VALUES IN A LIST; R.CONSISTENT FIRST AND R.ABSOLUTE SECOND.


    # CHECK THAT DESIGN IS BALANCED

    k.vec <- as.vector(table(factor(dataset[,2])))


    if ( any(k.vec!=k.vec[1]) ) {

      cat("ERROR: Design is not balanced and/or missing values")

      return()
    }


    # COMPUTE k

    k <- length(unique(dataset[,2]))



    # COMPUTE n

    n <- nrow(dataset)/k


    # CHECK THAT n IS AN INTEGER USING MODULUS FUNCTION

    if ( n%%1!=0 ) {
      cat("ERROR:  n is not a multiple of k")
      return()

    }

    # CALCULATING CS DISTANCES BETWEEN EACH PREDATOR

    d<-cs_distance.mat(dataset)

    model <- vegan::adonis(d~as.factor(dataset[,1]) + as.factor(dataset[,2]), data=dataset, permutations = 0)


    MSr <- model$aov.tab[1,3]
    MSe <- model$aov.tab[3,3]
    MSc <- model$aov.tab[2,3]


    rep.consistent <- (MSr - MSe)/(MSr + (k-1)*MSe)
    rep.absolute <- (MSr - MSe)/(MSr + (k-1)*MSe + (k/n)*(MSc-MSe))

    return(list(rep.consistent, rep.absolute, k))
  }

  #' Measure of (unadjusted) repeatability with missing values
  #'
  #' @param dataset data frame of diet estimates. Columns 1 and 2 of
  #'        \code{dataset} must be the first and second factors in your two factor
  #'        ANOVA model where column 1 specifies the predator and column 2 the
  #'        year, season, etc...  Missing values are allowed.
  #' @param gamma needed to compute CS distance.  Default is 1.
  #' @param k.est.meth default is 1 (harmonic mean).  See details below.
  #' @return unadjusted consistent repeatability, unadjusted absolute repeatability
  #'         and a modified measure of repeatability.
  #' @details We use the modified measure of repeatability that uses an adjusted
  #'          degrees of freedom and takes into account the estimate of \code{k}.
  #'          When \code{k.est.meth=1}, the harmonic mean is used (see
  #'          p. 212 from Biometry Fourth Edition by Sokal and Rohlf) and when
  #'          \code{k.est.meth=2},estimate is n_0 from Lessels and Boag (1987).
  #'
  #' @keywords internal

  unbal.two.factor.rep <- function(dataset,k.est.meth = 1, gamma=1){

    # CALCULATES CONSISTENT AND ABSOLUTE REPEATABILITY FROM A MATRIX OF PREDATOR
    # DIETS. THIS FUNCTION ALLOWS FOR MISSING VALUES.

    # COLUMNS 1 AND 2 OF DATASET MUST BE THE FIRST AND SECOND FACTORS IN YOUR TWO FACTOR
    # ANOVA MODEL WHERE COLUMN 1 SPECIFIES THE PREDATOR AND COLUMN 2 THE YEAR, SEASON, ETC...

    # k.est.meth = 1 if harmonic mean is to be used (see p. 212 from Biometry Fourth Edition by Sokal and Rohlf).
    # Otherwise, estimate is n_0 from Lessels and Boag (1987)

    # THE FUNCTION WILL RETURN TWO VALUES IN A LIST; R.CONSISTENT FIRST AND R.ABSOLUTE SECOND.


    n <- length(unique(dataset[,1]))

    d<-cs_distance.mat(dataset)

    model <- vegan::adonis(d~as.factor(dataset[,1]) + as.factor(dataset[,2]), data=dataset, permutations = 0)


    MSr <- model$aov.tab[1,3]
    MSe <- model$aov.tab[3,3]
    MSc <- model$aov.tab[2,3]
    SSc <- model$aov.tab[2,2]
    SSe <- model$aov.tab[3,2]

    counts <- data.frame(table(dataset[,1]))

    if (k.est.meth ==1) {

      k.est <- n/sum(1/counts[,2])

    } else {

      k.est <- (1/(nrow(counts)-1))*(sum(counts[,2])-sum(counts[,2]^2)/sum(counts[,2]))

    }


    df.c <- k.est-1
    df.e <- model$aov.tab[4,1]-model$aov.tab[1,1]-df.c


    unbal.r.consistent <- ((MSr - MSe)/k.est)/( (MSr-MSe)/k.est + MSe )
    unbal.r.absolute <- ( (MSr - MSe)/k.est )/( (MSr - MSe)/k.est + (MSc - MSe)/n + MSe )
    unbal.r.absolute.mod <- ( (MSr - (SSe/df.e))/k.est )/( (MSr - (SSe/df.e))/k.est + ( (SSc/df.c) - (SSe/df.e) )/n + (SSe/df.e) )

    return(list(unbal.r.consistent, unbal.r.absolute,unbal.r.absolute.mod))
  }

  #' CALCULATE CS DISTANCE BETWEEN EACH TWO PREDATORS IN dataset
  #'
  #' @keywords internal
  #'
  cs_distance.mat<-function(dataset,alpha=1) {

    # WRITTEN BY HONGCHANG BAO (JUNE 2018)

    # CALCULATE CS DISTANCE BETWEEN EACH TWO PREDATORS IN dataset

    #COLUMNS 1 AND 2 OF DATASET MUST BE THE FIRST AND SECOND FACTORS

    #store i,j into rec (i,j from the loop)
    res<-merge(data.frame(first=1:nrow(dataset)),
               data.frame(second=1:nrow(dataset)))

    dataset<-dataset[,-c(1,2)]

    dataset <- dataset/apply(dataset,1,sum)


    d.sq<-(dataset[res[,2],]-dataset[res[,1],])^2
    #d.sq
    c.vec <-dataset[res[,2],]+dataset[res[,1],]
    #c.vec

    if ( any(d.sq!=0) ) {
      d.sq[d.sq!=0] <- d.sq[d.sq!=0]/c.vec[d.sq!=0]
    }

    f.m<-function(d.sq.d,alpha=1)
    {
      nfa <- length(d.sq.d)
      CS.dist <- 1/alpha*sqrt(2*nfa)*sqrt(sum(d.sq.d))
      return(CS.dist)
    }

    r.m.y1 <- apply(d.sq,1,f.m)
    r.m.y1<-matrix(r.m.y1,nrow=nrow(dataset),byrow=TRUE)
    return(r.m.y1)
  }



  #'FUNCTION TO GENERATE ns PSEUDO SEALS WHERE ns IS THE NUMBER OF ROWS IN
  #'\code{diet.mat.seal} AND RETURNS THE CORRESPONDING QFASA DIET ESTIMATES.
  #'ASSUMES EACH ROW IS A DIET VECTOR AND EACH COLUMN CORRESPONDS TO 1 PREY
  #'SPECIES IN THE DIET
  #'
  #' @keywords internal
  pseudo.pred.rep <- function(diet.mat.seal, prey.database, fatcont.mat, dist.meas, ext.fa, gamma=1) {


    # SEPTEMBER 13TH, 2017
    # CHANGING TO JUST PASS MATRIX OF CALIBRATION FACTORS (QFASA DIET ESTIMATES APPEARED TO BE TOO VARIABLE AND INACCURATE)


    # WRITTEN BY MAGGIE BROWN (SUMMER 2015) AND MODIFIED BY MYSELF
    # AS OF AUGUST, 2017, CALLS A DIFFERENT pseudo.pred.  (OLD ONE IS CALLED pseudo.pred.traditional.)
    # USING MY pseudo.pred (SO NEEDED TO ADD specify.noise AND n.prey)

    # FUNCTION TO GENERATE ns PSEUDO SEALS WHERE ns IS THE NUMBER OF ROWS IN diet.mat.seal
    # AND RETURNS THE CORRESPONDING QFASA DIET ESTIMATES.
    # ASSUMES EACH ROW IS A DIET VECTOR AND EACH COLUMN CORRESPONDS TO 1 PREY SPECIES IN THE DIET

    # CALLIBRATIION COEFFICIENTS DIFFERENT FOR EACH INDIVIDIUAL SEAL AND INDIVIDUAL YEAR

    # INPUT:
    # diet.mat.seal: matrix where each row corresponds to one diet
    # cal.mat: matrix of calibration coefficients where each row is 1 FAs
    # n.noise: percent noise to be used in adding error (ie. 0.1)
    # fatcont.mat: matrix where first column is prey species and second column is corresponding fat content value
    # dist.meas: distance measure used in p.QFASA: KL (1), AIT (2), CS (3)
    # gamma: value for gamma used in CS distance measure in p.QFASA
    # ext.fa: extended dietary fatty acids
    # specify.noise: see pseudo.pred
    # tot.prey:  see pseudo.pred.table


    # determine species in prey database

    species <- unique(prey.database[,1])
    nfa <- ncol(prey.database[,-1])


    # initiate empty matrix where each row i will contain the generated FAs of diet i in diet.mat.seal

    diet.rep.mat <- matrix(rep(NA, nrow(diet.mat.seal)*ncol(diet.mat.seal)), nrow=nrow(diet.mat.seal), byrow=T)


    pseudo.pred.m.fa <- matrix(rep(NA, nrow(diet.mat.seal)*ncol(prey.database[,-1])), nrow=nrow(diet.mat.seal), byrow=T)

    for (m in 1:nrow(diet.mat.seal)) {

      #CALLIBRATION COEFFICIENTS DIFFERENT FOR EACH INDIVIDUAL SEAL AND INDIVIDUAL YEAR
      #randomly select 1 herring FAs and 1 seal FAs



      # NEW WAY

      cal.sim <- jitter(rep(1,nfa))
      cal.mod <- jitter(rep(1,nfa))




      #FAT CONTENT
      #intiate vectors to contain fat content values to generate seal i (fat.sim) and estimate diet of seal i (fat.mod)

      fat.sim <- rep(NA, length(species))
      fat.mod <- rep(NA, length(species))
      split.species.fat.content <- split(fatcont.mat, fatcont.mat[,1])

      for (j in 1:length(species)) {

        #select all the fat content values for species j

        fat.species.j <- data.frame(split.species.fat.content[j])
        colnames(fat.species.j) <- colnames(fatcont.mat)
        fatcont.species.j <- fat.species.j[,2]

        #remove any missing values present in fat content vector of species j

        fatcont.species.j <- fatcont.species.j[!is.na(fatcont.species.j)]

        #take a random sample of fat content values and split in half

        fatcont.species.sample <- sample(fatcont.species.j)
        fatcont.species.sample.1 <- fatcont.species.sample[1.:round(length(fatcont.species.j)/2.)]
        fatcont.species.sample.2 <- fatcont.species.sample[round(length(fatcont.species.j)/2.):length(fatcont.species.j)]
        fatcont.species.sample.1.mean <- mean(fatcont.species.sample.1)
        fatcont.species.sample.2.mean <- mean(fatcont.species.sample.2)

        #use half for generating pseudo seal FAs and half for estimating pseudo seal FAs
        fat.sim[j] <- fatcont.species.sample.1.mean
        fat.mod[j] <- fatcont.species.sample.2.mean
      }


      #for each row in diet.mat.seal, extract the diet of row i, generate FAs of seal i and estimate diet of seal i using p.QFASA
      #assign the first row as a diet vector

      #cat("diet.vec")
      diet.vec <- as.vector(unlist(diet.mat.seal[m,]))

      #generate FAs of pseudo seal m

      #pseudo.pred.m.fa<- pseudo.pred.traditional(prey.database, diet.vec, n.noise, tot.prey, cal.sim, fat.sim, specify.noise)

      #print(round(diet.vec,3))

      pseudo.pred.m.fa <- pseudo.pred(diet.vec,prey.database,cal.sim,fat.sim)


      #pseudo.pred.m.fa <- matrix(pseudo.pred.m.fa,nrow = 1)
      #pseudo.pred.m.fa <- as.data.frame(pseudo.pred.m.fa)
      #dimnames(pseudo.pred.m.fa)[[2]] <- dimnames(prey.database[,-1])[[2]]

      prey.base.ext.fattyacids<-prey.database[ext.fa]
      prey.base.ext.normalize <- prey.base.ext.fattyacids/(apply(prey.base.ext.fattyacids,1,sum))
      prey.base.ext <- cbind(prey.database[,1], prey.base.ext.normalize)
      prey.mat.ext <- MEANmeth(prey.base.ext)

      diet.rep.mat[m,]<- (p.QFASA(pseudo.pred.m.fa, prey.mat.ext, cal.mod, dist.meas, gamma,
                                  FC=fat.mod, ext.fa=ext.fa))[[1]]
      #cat("gen diet")
      #print(round(diet.rep.mat[m,],3))

    }

    return(diet.rep.mat)
  }


  #'FUNCTION TO GENERATE ns SEALS FOR EACH OF ny YEARS
  #' @keywords internal

  pseudo.pred.table <- function(diet.rep.mat,prey.database, fatcont.mat, dist.meas, ext.fa, gamma=1) {


    #function(diet.rep.mat, prey.database, n.noise, cal.herring.mat,
    #                              cal.seal.mat, fatcont.mat, dist.meas, gamma, ext.fa, tot.prey,
    #                              specify.noise) {


    # THIS IS MAGGIE BROWN'S (SUMMER 2015) var.seals.sim.err.noise FUNCTION RENAMED.
    # NOTE:  HER FUNCTION GENERATED n.sim "TABLES". THIS FUNCTION ONLY GENERATES ONE.


    #FUNCTION TO GENERATE ns SEALS FOR EACH OF ny YEARS

    #FUNCTION DIFFERENT FROM var.seals.sim.R/var.seals.sim.2.R AS THIS FUNCTION ADDS ERROR TO
    #SEALS DIETS BY ADDING NOISE

    # INPUT:
    # diet.mat: matrix where each row corresponds to one diet, assumes first column is seal ID and second column is year ID
    # NOTE:  SEAL ID MUST BE OF FORM 1,1,1,1....,2,2,2,2... AND YEAR MUST BE 1,2,3,...,1,2,3...
    # prey.database: prey database used to generate/estimate seal FAs/diet, assumes first column is species ID
    # n.noise: percent noise to be used in adding error (ie. 0.1)
    # cal.herring.mat: matrix containing herring FAs (each row is 1 FAs)
    # cal.seal.mat: matrix containing seal FAs (each row is 1 FAs)
    # fatcont.mat: matrix of fat content values where first column is prey species and second column is corresponding fat content value
    # dist.meas: distance measure used in p.QFASA: KL (1), AIT (2), CS (3)
    # gamma: value for gamma used in CS distance measure in p.QFASA
    # ext.fa: extended dietary fatty acids
    # tot.prey:  Total number of prey to use in creating pseudo seal (use to be called nprey)
    # specify.noise:  see pseudo.pred.rep


    #determine number of seals and prey in diet matrix

    n.seals <- length(unique(diet.rep.mat[,1]))
    n.prey <- ncol(diet.rep.mat[,c(-1,-2)])

    #split diet matrix (diet.rep.mat) according to seal ID

    split.seals <- split(as.data.frame(diet.rep.mat), diet.rep.mat[,1])


    #initiate matrix to contain sample
    sim.seal.diet.mat <- matrix(rep(NA), nrow=0, ncol=n.prey)

    for (p in 1:n.seals) {


      #cat("Seal")

      #extract seal ID p
      seal.p <- data.frame(split.seals[p])

      colnames(seal.p) <- colnames(diet.rep.mat)

      seal.p <- seal.p[,c(-1,-2)]


      # generate FAs and estimate diet using p.QFASA



      sim.seal.p.diet.mat<- pseudo.pred.rep(seal.p, prey.database, fatcont.mat, dist.meas, ext.fa, gamma)


      sim.seal.diet.mat <- rbind(sim.seal.diet.mat, sim.seal.p.diet.mat)


    }

    #build sample  and store in list

    # COMMENTED CODE WAS ADDED BY ME IN CASE SEAL ID AND YEAR ID NOT
    # IN PROPER FORMAT BUT DON'T THINK IT WILL WORK WHEN USING MISSING VALUES.

    #n.years <- length(unique(diet.rep.mat[,2]))
    #seal.col <- rep(1,n.years)
    #for (i in 2:n.seals) {
    #seal.col <- c(seal.col,rep(i,n.years))
    #}

    #year.col <- rep(seq(1,n.years,1),n.seals)


    sim.seal.diet.sample.mat <- cbind(diet.rep.mat[,c(1,2)], sim.seal.diet.mat)

    #sim.seal.diet.sample.mat <- cbind(seal.col,year.col, sim.seal.diet.mat)


    colnames(sim.seal.diet.sample.mat) <- colnames(diet.rep.mat)

    return(sim.seal.diet.sample.mat)

  }

  #'RETURNS RHO AND THE BOOTSTRAP SAMPLE DETERMINED BY ind.
  #' @keywords internal


  rho.boot.fun <- function(ind,data,gamma.rho) {

    # RETURNS RHO AND THE BOOTSTRAP SAMPLE DETERMINED BY ind.

    # THIS IS THE FUNCTION THAT WILL BE CALCULATED
    # ind IS THE INDICES NEEDED FOR THE BOOTSTRAPPING
    # dat seal.ids <- unique(data[,1])


    ns <- length(ind)
    k <- length(unique(data[,2]))

    p.mat.star<- data[data[,1]==ind[1],]

    for (j in 2:ns) {

      p.mat.star <- rbind(p.mat.star,data[data[,1]==ind[j],])

    }

    # MAKING SURE SAMPLE TREATED AS A RANDOM SAMPLE OF SIZE ns

    seq.seal<- rep(1,k)
    for (i in 2:ns)
    { seq.seal <- c(seq.seal,rep(i,k))
    }

    p.mat.star[,1] <- seq.seal

    rho.hat <- two.factor.rep(p.mat.star,gamma.rho)[[2]]


    return(list(rho.hat, p.mat.star))


  }

  #'FUNCTION TO BE PASSED TO bootstrap::jackknife.
  #' @keywords internal
  #'
  rho.jack.fun <- function(ind,data,gamma.rho) {

    # FUNCTION TO BE PASSED TO bootstrap::jackknife.

    # COMPUTES rho USING OBSERVATIONS CORRESPONDING TO INDICES IN ind.

    # RETURNS rho.


    k <- length(unique(data[,2]))

    p.mat.star<- data[data[,1]==ind[1],]

    for (j in 2:length(ind)) {

      p.mat.star <- rbind(p.mat.star,data[data[,1]==ind[j],])

    }

    # MAKING SURE SAMPLE TREATED AS A RANDOM SAMPLE OF SIZE ns


    rho.hat <- two.factor.rep(p.mat.star,gamma.rho)[[2]]


    return(rho.hat)


  }

  #' RETURNS RHO AND THE BOOTSTRAP SAMPLE DETERMINED BY ind.
  #' @keywords internal
  #'

  unbal.rho.boot.fun <- function(ind,data,gamma.rho) {

    # RETURNS RHO AND THE BOOTSTRAP SAMPLE DETERMINED BY ind.

    # THIS IS THE FUNCTION THAT WILL BE CALCULATED
    # ind IS THE INDICES NEEDED FOR THE BOOTSTRAPPING
    # dat seal.ids <- unique(data[,1])

    ns <- length(ind)
    k.vals <- rep(NA,ns)
    p.mat.star<- data[data[,1]==ind[1],]

    k.vals[1] <- nrow(p.mat.star)


    for (j in 2:ns) {

      p.mat.star.curr <- data[data[,1]==ind[j],]
      k.vals[j] <- nrow(p.mat.star.curr)

      p.mat.star <- rbind(p.mat.star,p.mat.star.curr)

    }

    # MAKING SURE SAMPLE TREATED AS A RANDOM SAMPLE OF SIZE ns


    p.mat.star[,1] <- rep(seq(1,ns,1),k.vals)


    rho.hat <- unbal.two.factor.rep(p.mat.star,gamma.rho)[[3]]


    return(list(rho.hat, p.mat.star))


  }

  #'FUNCTION TO BE PASSED TO bootstrap::jackknife.
  #' @keywords internal
  #'
  unbal.rho.jack.fun <- function(ind,data,gamma.rho) {

  # FUNCTION TO BE PASSED TO bootstrap::jackknife.

  # COMPUTES rho USING OBSERVATIONS CORRESPONDING TO INDICES IN ind.

  # RETURNS rho.

  p.mat.star<- data[data[,1]==ind[1],]

  for (j in 2:length(ind)) {

    p.mat.star <- rbind(p.mat.star,data[data[,1]==ind[j],])

  }


  rho.hat <- unbal.two.factor.rep(p.mat.star,gamma.rho)[[3]]


  return(rho.hat)


  }






