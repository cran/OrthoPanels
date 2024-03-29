% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opm.R
\name{DIC}
\alias{DIC}
\title{Deviance Information Criterion (DIC)}
\usage{
DIC(object, ...)
}
\arguments{
\item{object}{an instance of class \code{opm} whose DIC is wanted.}

\item{...}{further arguments passed to other methods.}
}
\value{
a numeric value with the corresponding DIC
}
\description{
Computes the Deviance Information Criterion (DIC), which is a
generalization of the Akaike Information Criterion. Models with
smaller DIC are considered to fit better than models with larger
DIC.
}
\details{
DIC is defined as \eqn{DIC = 2*\bar{D} - D_\theta}
where:
\eqn{\bar{D} = -2 mean(log-likelihood at parameter samples)}
\eqn{D_\theta = -2 * log(likelihood at expected value of parameters)}

DIC is calculated as: \code{2 * (-2 * mean(log-likelihood at each element of parameter samples)) - (-2 * log(likelihood at mean parameter sample value))}
}
\note{
Note the speed of computation of the DIC in proportional to
    the number of sampled values of the parameters in the opm
    object.
}
