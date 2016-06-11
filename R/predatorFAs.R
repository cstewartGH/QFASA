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
#' \item{c12.0}{}
#' \item{c13.0}{}
#' \item{Iso14}{}
#' \item{c14.0}{}
#' \item{c14.1w9}{}
#' \item{c14.1w7}{}
#' \item{c14.1w5}{}
#' \item{Iso15}{}
#' \item{Anti15}{}
#' \item{c15.0}{}
#' \item{c15.1w8}{}
#' \item{c15.1w6}{}
#' \item{Iso16}{}
#' \item{c16.0}{}
#' \item{c16.1w11}{}
#' \item{c16.1w9}{}
#' \item{c16.1w7}{}
#' \item{c7Mec16.0}{}
#' \item{c16.1w5}{}
#' \item{c16.2w6}{}
#' \item{Iso17}{}
#' \item{c16.2w4}{}
#' \item{c16.3w6}{}
#' \item{c17.0}{}
#' \item{c16.3w4}{}
#' \item{c17.1}{}
#' \item{c16.4w3}{}
#' \item{c16.4w1}{}
#' \item{c18.0}{}
#' \item{c18.1w13}{}
#' \item{c18.1w11}{}
#' \item{c18.1w9}{}
#' \item{c18.1w7}{}
#' \item{c18.1w5}{}
#' \item{c18.2d5.11}{}
#' \item{c18.2w7}{}
#' \item{c18.2w6}{}
#' \item{c18.2w4}{}
#' \item{c18.3w6}{}
#' \item{c18.3w4}{}
#' \item{c18.3w3}{}
#' \item{c18.3w1}{}
#' \item{c18.4w3}{}
#' \item{c18.4w1}{}
#' \item{c20.0}{}
#' \item{c20.1w11}{}
#' \item{c20.1w9}{}
#' \item{c20.1w7}{}
#' \item{c20.2w9}{}
#' \item{c20.2w6}{}
#' \item{c20.3w6}{}
#' \item{c20.4w6}{}
#' \item{c20.3w3}{}
#' \item{c20.4w3}{}
#' \item{c20.5w3}{}
#' \item{c22.1w11}{}
#' \item{c22.1w9}{}
#' \item{c22.1w7}{}
#' \item{c22.2w6}{}
#' \item{c21.5w3}{}
#' \item{c22.4w6}{}
#' \item{c22.5w6}{}
#' \item{c22.4w3}{}
#' \item{c22.5w3}{}
#' \item{c22.6w3}{}
#' \item{c24.1w9}{}
#' }
#'
"predatorFAs"
