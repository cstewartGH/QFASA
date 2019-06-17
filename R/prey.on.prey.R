
#'
#' Each prey fatty acid signature is systematically removed from the
#' supplied prey database and its QFASA diet estimate is obtained by
#' treating the individual as a predator.
#'
#' @export
#' 
#' @param preybase first column is name of species and remaining columns are fatty acids.
#' @param dist.meas see help file for \code{\link{p.QFASA}}.
#' @param gamma  see help file for \code{\link{p.QFASA}}.
#' @return diet estimate
#'
#' @examples
#' data(preyFAs)
#' my.preybase <- preyFAs[, -c(1,3)]
#'
#' # Note: uncomment examples to run. CRAN tests fail because execution time > 5 seconds
#' # diets.out <- prey.on.prey(my.preybase, 2)
#' # round(MEANmeth(diets.out), 3)
#'      
prey.on.prey <-function(preybase,dist.meas,gamma=1)
{  
  preybase[, -1] <- preybase[,-1]/apply(preybase[, -1], 1, sum)
  sort.preytype <- order(preybase[, 1])
  preybase <- preybase[sort.preytype,]
  groups <- factor(preybase[, 1], unique(preybase[, 1]))
  
  for( j in 1:nrow(preybase))
  {
      pred.sig <-preybase[j, -1] 
      prey.matrix <- MEANmeth(preybase[-j, ])
      result <- p.QFASA(pred.sig,
                        prey.matrix,
                        cal.mat = rep(1,length(pred.sig)),
                        dist.meas,
                        gamma,
                        ext.fa = colnames(prey.matrix))
      
      diet.est <- as.data.frame(result$`Diet Estimates`)
      
      if(j == 1) {
          DE.df<-diet.est
      }
      else {
          DE.df <- rbind(DE.df, diet.est)
      }
  }
    
  DE.df <- cbind(groups, DE.df)
  dimnames(DE.df)[[2]] <- c("Species", levels(groups))    
  return(DE.df)
}

