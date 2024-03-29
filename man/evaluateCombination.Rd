% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/All_Selection_Algorithm_Code.R
\name{evaluateCombination}
\alias{evaluateCombination}
\title{Get nll and IC values for a pair of species (this is used in forward
selection if the starting species are not specified; used on a list of all
possible combinations of two species)}
\usage{
evaluateCombination(
  species.index,
  species.list,
  pred.mat,
  prey.mat,
  cal.vec,
  fat.vec,
  ext.fa,
  k = 2
)
}
\description{
Get nll and IC values for a pair of species (this is used in forward
selection if the starting species are not specified; used on a list of all
possible combinations of two species)
}
\keyword{internal}
