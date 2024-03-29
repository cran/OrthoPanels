% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opm.R
\name{longRunEffects}
\alias{longRunEffects}
\title{Long run effects based on the \code{opm} Model Parameters}
\usage{
longRunEffects(opm_obj, parm = NULL, probs = c(0.025, 0.5, 0.975))
}
\arguments{
\item{opm_obj}{an instance of class \code{opm}}

\item{parm}{a specification of which parameters are to be plotted,
a vector of names are the only legal values. If missing, all parameters are considered.}

\item{probs}{a vector of specified quantiles, by default, the c(0.025,0.5,0.975) are ("\code{probs}")}
}
\value{
A matrix with quantiles on the rows, with number of rows specified as length of the \code{probs} vector for the specified quantiles, with covariates on the columns
}
\description{
Computes long run effects and confidence intervals of \code{opm} Model Parameters
}
\examples{
\dontrun{
longRunEffects(opm_obj)
longRunEffects(opm_obj,probs=c(0.975, 0.16, 0.5, 0.84, 0.025))
}

}
