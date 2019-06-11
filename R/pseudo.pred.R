#' Generate a pseudo predator by sampling with replacement from prey
#' database.
#'
#' Note: if preysize=1, then one prey is selecting from each species.
#'       otherwise, a sample of size n_k (number of species k) is
#'       sampled with replacement.
#'
#' @export 
#' @param diet compositional vector of proportions that sums to one.
#'     Length is equal to the number of prey species.  
#' @param preybase prey database with first column providing the species name.
#' @param cal.vec vector of calibration coefficients.
#' @param fat.vec vector of fat content whose length is the same as the number of species.
#' @param preysize number of prey to sample from prey database.
#' @return a simulated predator FA signature
#'
#' @examples
#' data(preyFAs)
#' p.mat <- matrix(rep(NA,100*11),nrow=100)
#' for (i in 1: 100) {
#'     my.seal <- pseudo.pred(rep(1/11,11),
#'                            preyFAs[,-c(1,3)],
#'                            rep(1,ncol(preyFAs[,-c(1,3)])-1),
#'                            rep(1,11))
#'     p.mat[i,] <- p.QFASA(my.seal,
#'                          MEANmeth(preyFAs[,-c(1,3)]),
#'                          rep(1,length(my.seal)),
#'                          2,
#'                          ext.fa=colnames(preyFAs[,-c(1:3)]))$`Diet Estimates`
#' }
#' 
#' # Average diet estimate 
#' round(apply(p.mat,2,mean),3)
#' 
pseudo.pred <- function(diet, preybase, cal.vec, fat.vec, preysize=2) {
    
  preybase[,-1] <- preybase[,-1]/apply(preybase[,-1],1,sum)
  sort.preytype <- order(preybase[,1])
  preybase <- preybase[sort.preytype,]
  
  
  diet <- fat.vec * diet
  diet <- diet/sum(diet)
  
  
  I <- length(unique(preybase[,1]))
  nk.vec <- tapply(preybase[,1],preybase[,1],length)
  pred <- matrix(rep(NA,(ncol(preybase)-1)*I),nrow=ncol(preybase)-1,ncol=I,byrow=T)
  diet <- unlist(diet)
  
  cum <- cumsum(nk.vec)
  
  preybase.mat <- as.matrix(preybase[,-1])
  
  if (preysize==1) {
    
    # SELECTING ONE PREY SIGNATURE
    
    for (k in 1:I) {
      
      if (k==1) {
        
        
        ind <- sample(seq(1,cum[k],1),size=1)
        pred[,k] <- diet[k] * preybase.mat[ind, ]
        
      } else {
        
        
        # SELECTING ONE PREY SIGNATURE
        ind <- sample(seq((cum[k-1]+1),cum[k],1),size=1)
        pred[,k] <- diet[k] * preybase.mat[ind,]
        
      }
      
    } # end loop
    
  } else {
    
    # SAMPLING PREY DATA BASE
    
    for (k in 1:I) {
      
      if (k==1) {
        
        inds <- sample(seq(1,cum[k],1),size=nk.vec[k],replace=TRUE)
        pred[,k] <- diet[k] * apply(preybase.mat[inds, ],2,mean)
        
      } else {
        
        inds <- sample(seq((cum[k-1]+1),cum[k],1),size=nk.vec[k], replace=TRUE)
        pred[,k] <- diet[k] * apply(preybase.mat[inds,],2,mean)
      }
      
    } # end loop
    
  } # end outer if statement
  
  y.hat <- apply(pred,1,sum)*cal.vec
  
  y.hat <- y.hat/sum(y.hat)
  
  y.hat <- matrix(y.hat,nrow = 1)
  y.hat <- as.data.frame(y.hat)
  dimnames(y.hat)[[2]] <- dimnames(preybase[,-1])[[2]]
  
  return(y.hat)
}
