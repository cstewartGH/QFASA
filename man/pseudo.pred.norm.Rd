% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/All_MUFASA_Code.R
\name{pseudo.pred.norm}
\alias{pseudo.pred.norm}
\title{Generate a pseudo predator parametrically from multivariate normal
distributions.}
\usage{
pseudo.pred.norm(mu.mat, sigma.pool, diet)
}
\arguments{
\item{mu.mat}{matrix where each row represents the mean transformed FA
signature of each prey type}

\item{sigma.pool}{pooled variance-covariance matrix of the transformed
fatty acid signatures of prey types}

\item{diet}{the "true" or "desired" diet of the pseudo predator with prey
species in the same order as the rows of mu.mat. A compositional vector of
proportions that sums to one with length equal to
the number of prey species.}
}
\value{
A simulated predator FA signature.  See \emph{pseudo.pred} for an
example illustrating how to generate a sample of pseudo predators.
}
\description{
Generate a pseudo predator parametrically from multivariate normal
distributions.
}
\details{
Similar to \emph{pseudo.pred} but instead generates the
pseudo-predators parametrically by assuming ilr transformed FA signatures
have a multivariate normal distribution.
}
