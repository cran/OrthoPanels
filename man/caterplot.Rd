% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opm.R
\name{caterplot}
\alias{caterplot}
\title{Caterpillar Plots of \code{opm} Model Parameters}
\usage{
caterplot(
  x,
  parm,
  main = paste("Caterpillar plot of", xname),
  xlab = "Range of parameter samples",
  labels = colnames(ranges)
)
}
\arguments{
\item{x}{an instance of class \code{opm}}

\item{parm}{a specification of which parameters are to be plotted,
       either a vector of names ("\code{rho}", "\code{sig2}"
and "\code{beta}" are the only legal values) or a vector of positional
indices. If missing, all parameters are considered.}

\item{main, xlab}{useful defaults for the plot title and X-axis label}

\item{labels}{labels for each parameter's interval: see \code{\link[graphics]{axis}}}
}
\value{
A matrix of 2.5\%, 5\%, 50\%, 95\%, and 97.5\% quantiles for
each of the desired parameters, with parameters arranged in
columns.
}
\description{
Creates side-by-side plots of equal-tailed credible intervals of \code{opm}
model parameters. The intervals are displayed as horizontal lines,
with 90\% interval using a thicker line width and 95\% interval a
thinner one. The posterior median is indicated with a dot.
}
\examples{
\dontrun{
caterplot(o, main = NULL, labels = expression(alpha, beta, sigma^2))
}

}
