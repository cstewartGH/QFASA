% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Confidence_Interval_Beta_Distribution.R
\name{bisect.beta.lim}
\alias{bisect.beta.lim}
\title{Find simultaneous and individual confidence intervals for diet
proportions of a single prey species i.e. solve f(pio) = PVAL(pio)
= alpha1 and f(pio) = PVAL(pio) = alpha2 using bisection.}
\usage{
bisect.beta.lim(alpha1, alpha2, par.list, R, p.mat, k)
}
\arguments{
\item{alpha1}{simultaneous confidence level}

\item{alpha2}{individual confidence level}

\item{par.list}{a list of R.p lists of I beta distribution parameters phi
and theta that define diet proportion estimates for each of the
prey species. Effectively R.p beta distibutions for each of the
I prey species from which we bootstrap to calculate p-values.}

\item{R}{number of bootstrap replicates to use in p-value estimation.}

\item{p.mat}{predator diet estimates.}

\item{k}{prey species index 1..I}
}
\value{
upper and lower limits for individual and simultaneous
    intervals respectively.
}
\description{
Assumptions: alpha1 < alpha2
}
\details{
Note:  Tried to minimize number of times have to compute a pvalue
since very slow.
}
\keyword{internal}
