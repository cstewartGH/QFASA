#' Predator fatty acid signatures. Each predator signature is a row
#' with fatty acid proportions in columns.
#'
#' Fatty acid signatures are subsetted for the chosen fatty acid set
#' and renormalized during the modelling so there is no need to subset
#' and/or renormalize prior to running p.QFASA. However, make sure that
#' the the same fatty acids appear in the predator and prey files (if a
#' FA appears in one but not the other the code will give you an
#' error).
#'
#' Unlike the original QFASApack code the predator data can
#' contain as much tombstone data in columns as you wish but the
#' predator FA signatures must be extracted as a separate input in
#' order to run in p.QFASA. 
#'
#' @format A data frame with 10 observations and 70 variables:
#' \describe{
#' \item{SampleCode}{TODO}
#' \item{AnimalCode}{TODO}
#' \item{SampleGroup}{TODO}
#' \item{Biopsy}{TODO}
#' \item{fatty acids}{TODO}
#' }
#'
"predatorFAs"
