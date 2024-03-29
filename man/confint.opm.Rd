% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opm.R
\name{confint.opm}
\alias{confint.opm}
\title{Credible Intervals for Model Parameters}
\usage{
\method{confint}{opm}(object, parm, level = 0.95, ...)
}
\arguments{
\item{object}{an instance of class \code{opm} whose credible
intervals are wanted}

\item{parm}{a specification of which parameters are to be given
       credible intervals, either a vector of names ("\code{rho}", "\code{sig2}",
and "\code{beta}" are the only legal values) or a vector
       of positional indices. If missing, all parameters are
       considered.}

\item{level}{the size of the interval (e.g., 0.95 for 95\% C.I.)}

\item{...}{additional argument(s) for methods}
}
\value{
A matrix with columns giving lower and upper limits of the
credible interval for each parameter. These will be labeled as (1 -
level/2) and 1 - (1 - level)/2 in \% (by default, "2.5\%" and
"97.5\%").
}
\description{
Computes equal-tailed credible intervals for one or more parameters
in a fitted \code{opm} model. The method used is the quantile
interval of the posterior sample.
}
\seealso{
\code{\link{confint}}
}
