% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OrthoPanels-package.R
\docType{package}
\name{OrthoPanels-package}
\alias{OrthoPanels-package}
\title{OrthoPanels: Orthogonalized Panel Model}
\description{
This package includes the function \code{opm()}, which implements the
orthogonal reparameterization approach recommended by Lancaster
(2002) to estimate dynamic panel models with fixed effects (and
optionally: panel specific intercepts). The OLS estimator for such
models is biased with a fixed (small) \eqn{N} (Nickell 1981). Equivalently,
a maximum likelihood estimation leads to an incidental parameters
problem (Neyman and Scott 1948; Lancaster 2000). The approach by
Lancaster (2002) uses an orthogonal reparameterization of the fixed
effects to produce a likelihood-based estimator of the remaining
parameters that is exact and consistent as \eqn{N} approaches infinity
for \eqn{T} greater than or equal to 2.
}
\details{
Orthopanels can accomodate unbalanced panel data, in that some
respondents may drop out early (attrition) and some respondents may
enter the panel late (refreshment). It is assumed that once
respondents enter the panel, they will have observations up until
they dropout and then NAs in subsequent waves. The estimation is
conducted under the assumption that the data is missing at random
}
\references{
Lancaster, T. (2000) The incidental parameter problem since 1948.
\emph{Journal of Econometrics}, \bold{95}, 391--413.

Lancaster, T. (2002) Orthogonal parameters and panel data.
\emph{Review of Economic Studies}, \bold{69}, 647--666.

Neyman, J. and Scott, E. L. (1948) Consistent estimation from
partially consistent observations. \emph{Econometrica}, \bold{16},
1--32.

Nickell, S. (1981) Biases in dynamic models with fixed effects.
\emph{Econometrica}, \bold{49}, 1417--1426.
}
