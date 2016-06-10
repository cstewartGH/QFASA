#' Prey fatty acid signatures. Each prey signature is a row with fatty
#' acid proportions in columns.
#'
#' The prey file should contain all of the individual fatty acid signatures of
#' the prey and their lipid contents (where appropriate) - a matrix of
#' the mean values for the FAs (prey.matrix) by the designated prey
#' modelling group is then calculated using the MEANmeth function.
#'
#' Like the predator .csv file you can have as many tombstone data
#' columns as required but there must be at least one column that
#' identifies the modelling group, in this case, Species.
#'
#' Unlike the predator data, the prey data is not subsetted and
#' renomalized during the modelling so the prey file needs to be
#' subsetted for the desired fatty acid set and renormalized to
#' sum to 1 prior to calculating the mean values.
#' 
#' The full FA set is extracted from the data frame
#' (columns 4 onward), subsetted for the FA set in use and then
#' renormalized over 1. The modelling group names (the "Species" column
#' in this case) is then added back to the subsetted and renormalized
#' data (as the first column) and the average values calculated using
#' the MEANmeth function. Note that for the MEANmeth function to work
#' the modelling group name must be in the first column.
#'
#' @format A data frame with 302 observations and 70 variables:
#' \describe{
#' \item{Lab.Code}{TODO}
#' \item{Species}{TODO}
#' \item{lipid}{TODO}
#' \item{fatty acids}{TODO}
#' }
#'
"preyFAs"
