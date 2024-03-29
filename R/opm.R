#' Fitting orthogonal panel models
#'
#' \code{opm} is used to fit orthogonal panel models.
#' 
#' The model can be either specified symbolically with the formula
#' \code{response ~ term1 + term2 ...} or with the terms and response
#' given as a pair of 3- and 2-dimensional arrays, \code{x} and
#' \code{y} respectively. The arrays have to be in the format
#' \code{time x variable x case} for terms and \code{time x case} for
#' the response.
#'
#' The lagged dependent variable does not need to be included in the
#' formula or data, as it is included automatically.
#' 
#' @param x a formula (see description of parameter \code{formula}
#' below) or an array of dimension \code{time x variable x case} of terms.
#' @param n.samp number of samples to use to estimate the parameters.
#' @param ... further arguments passed to other methods.
#'
#' @return An object of class \code{opm} with the following elements:
#' \describe{
#' \item{\code{samples}}{parameter samples used to estimate the model,
#' as a list with following elements:}
#' \describe{
#' \item{\code{rho}}{a vector of \code{n.samp} samples of \eqn{\rho}.}
#' \item{\code{v}}{a vector of \code{n.samp} samples of \eqn{\frac{1}{\sigma^2}}.}
#' \item{\code{beta}}{an \code{n.samp x variable} matrix of samples of \eqn{\beta}.}
#' }
#' \item{\code{call}}{the matched call}
#' \item{\code{index}}{the index variables, when using the formula interface}
#' \item{\code{time.indicators}}{\code{TRUE} if dummy time variables are used (see Notes), \code{FALSE} otherwise}
#' \item{\code{terms}}{the \code{terms} object used}
#' }
#'
#' The function \code{summary} (i.e., \code{summary.opm}) can be used
#' to obtain or print a summary of the results. The generic accessor
#' functions \code{coefficients}, \code{fitted.values},
#' \code{residuals}, \code{logLik}, and \code{df.residual} can be used
#' to extract various useful features of the value returned by \code{opm}.
#' 
#' @examples
#' set.seed(123)
#' N <- 5
#' T <- 2
#' beta <- .5
#' rho <- .5
#' v <- 1
#'
#' f <- runif(N, -2, 2)
#' K <- length(beta)
#' beta <- matrix(beta, K, 1)
#'
#' ## $x_i = 0.75 f + N(0, 1)$:
#' x <- array(.75*f, dim=c(N, K, (T+1))) + rnorm(N*K*(T+1))
#'
#' ## $y_{i,t} = \rho y_{i,t-1} + \beta x_{i,t} + f_i + N(0,1)$:
#' y <- matrix(0, N, T+1)
#' for (t in seq_len(T+1)) {
#'     yy <- if (t>1) y[,t-1] else 0
#'     y[,t] <- rho * yy + f  + x[,,t] %*% beta + rnorm(N, sd = sqrt(1/v))
#' }
#'
#' d <- data.frame(i = rep(seq(N), T+1),
#'                 t = rep(seq(T+1), each = N),
#'                 as.data.frame(matrix(aperm(x, c(1, 3, 2)), N*(T+1), K,
#'                                      dimnames = list(NULL, paste0('x', seq(K))))),
#'                 y = c(y))
#' opm(y~x1, d, n.samp = 10)
#'
#' @export
opm <- function(x, ...) {
    UseMethod('opm')
}


#' @rdname opm
#' @param y a matrix of dimensions \code{time x case} of responses.
#' @param add.time.indicators (logical) if \code{TRUE}, adds dummy
#' variables for time.
#' @note Dummy time variables exist as an additional column for each
#' wave of data, excluding the first and second wave (i.e., at
#' \eqn{t=0} and \eqn{t=1} using the terminology from Lancaster
#' (2000)). The new variables are named \code{tind.}\eqn{t}, where
#' \eqn{t = 2, ...}, and appear as such as elements of the estimated
#' \code{beta} coefficient.
#' @importFrom stats median
#' @export
opm.default <- function(x, y, n.samp, add.time.indicators = FALSE, ...) {
    cl <- match.call()
    ## clean up the call name to refer to the generic
    cl[[1]] <- as.name('opm')
    
    aligned_xy <- align_latecomers(x, y)
    x <- aligned_xy$x
    y <- aligned_xy$y
    
    if (add.time.indicators) {
        x <- with_time_indicators(x)
    }
    
    sample_params <- sample_all(x, y, n.samp)
    
    Ti <- Ti(x, y) - 1                  # x and y are not centered yet
    rho_ <- median(sample_params$rho)
    beta_ <- apply(sample_params$beta, 2, median)
    v_ <- 1/median(sample_params$sig2)
    fitted <- fitted(x, y, rho_, beta_)
    residuals <- center_yt(y, Ti) - fitted
    df.residual <- sum(ifelse(Ti>0, Ti, 0)) - (ncol(x)+1)
    logLik <- log_likelihood(x, y, rho_, matrix(beta_, ncol(x), 1), v_)

    design <- if (all(Ti == nrow(x)-1)) "balanced" else "unbalanced (with dropouts)"
    
    structure(list(samples = sample_params,
                   fitted.values = fitted,
                   residuals = residuals,
                   df.residual = df.residual,
                   logLik = logLik,
                   design = design,
                   call = cl,
                   .Environment = parent.frame(),
                   time.indicators = add.time.indicators && TRUE),
              class = 'opm')
}


#' @rdname opm
#' @param data an optional data frame, list, or environment containing
#' the variables in the model. If not found in \code{data}, the
#' variables are taken from \code{environment(x)}, typically the
#' environment from which \code{opm} is called.
#' @param subset an optional vector specifying a subset of
#' observations to be used in the fitting process.
#' @param index a two-element vector containing the index of the case
#' and time variables, respectively. Variable indices can be specifed
#' by name or position. This argument is ignored if the model is not
#' specified by the formula, because the index is implicit in the
#' organization of the terms and response arrays.
#' @importFrom stats as.formula model.matrix model.response xtabs
#' @export
opm.formula <- function(x, data = environment(x), subset = NULL, index = 1:2, n.samp, ...) {
    cl <- match.call()
    ## clean up the call name to refer to the generic
    cl[[1]] <- as.name('opm')
    
    with(reshape_inputs(x, data, subset, index, parent.frame()), {
        modifyList(opm.default(x, y, n.samp, ...),
                   list(call = cl,
                        index = vapply(index, function(ix) {
                                    if (is.numeric(ix)) names(data)[ix] else ix
                                }, character(1), USE.NAMES = FALSE),
                        terms = mt,
                        .Environment = attr(mt, '.Environment')))
    })
}


## Extend responses with T-2 columns of dummy variables for time
with_time_indicators <- function(x) {
    ## temporary rearrange to make adding dummies easier
    x <- aperm(x, c(1, 3, 2))
    dims <- dim(x)
    T <- dims[1]
    N <- dims[2]
    K <- dims[3]
    K_dum <- T-2
    dims[3] <- K + K_dum
    
    dimnams <- dimnames(x)
    if (is.null(dimnams)) {
        dimnams <- list(NULL, NULL, NULL)
    }
    if (is.null(dimnams[[3]])) {
        dimnams[[3]] <- seq(K)
    }
    dimnams[[3]] <- c(dimnams[[3]], paste('tind', seq(K_dum)+1, sep='.'))
    
    dummys <- rep(c(rbind(matrix(0, 2, K_dum), diag(K_dum))), N)
    x <- array(c(x,
                 dummys),
               dims,
               dimnams)

    ## rearrange back
    x <- aperm(x, c(1, 3, 2))
}


## Convert the input formula and data to a 3-d array for predictors and a matrix for the response
##
## arguments:
## - formula: formula specifying the model
## - data: environment containing the variables in the model
## - subset: subset of observations to use (NULL=all)
## - index: index of the case and time variables
## - envir: environment in which the model.frame should be evaluated
##
## returns a list object with elements:
## - x: 3-d array for predictors
## - y: matrix for responses
## - mt: the 'terms' object used
#' @importFrom stats reshape
reshape_inputs <- function(formula, data, subset, index, envir) {
    mf <- eval(bquote(stats::model.frame(formula = .(formula),
                                         data = .(data),
                                         subset = .(subset),
                                         na.action = na.pass)),
               envir)


    
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0L
    i <- factor(data[[index[1]]])
    t <- factor(data[[index[2]]])
    y <- model.response(mf, "numeric")
    yf <- data.frame(i = i,
                     t = t,
                     y = y)
    y <- xtabs(y~t+i, yf)
    class(y) <- 'matrix'

    yf$i <- as.integer(yf$i)
    yf$t <- as.integer(yf$t)
    y_nas <- as.matrix(yf[is.na(yf$y), 2:1])
    y[y_nas] <- NA

    missing_obs <- reshape(as.data.frame(table(i, t) == 0),
                           direction = 'long',
                           varying=list(levels(t)),
                           timevar = 't',
                           idvar='i',
                           v.names='count')
    missing_obs <- missing_obs[missing_obs[['count']], c('t', 'i')]
    y[as.matrix(missing_obs)] <- NA

    
    x <- model.matrix(mt, mf)
    xf <- data.frame(i = i, t = t,
                     as.data.frame(x))
    ncols=ncol(xf)
    x <- xtabs(as.formula(paste0('cbind(',
                                 paste(colnames(xf)[3:ncols], collapse=','),
                                 ') ~ t+i')),
               xf)
    class(x) <- 'array'
    xf$i <- as.integer(xf$i)
    xf$t <- as.integer(xf$t)
    x_nas <- as.matrix(xf[apply(xf[, -(1:2), drop=FALSE], 1,
                                function(row)any(is.na(row))), 2:1])
    
    if (length(dim(x)) == 2) dim(x) <- c(dim(x), 1)
    for (k in seq_len(dim(x)[3])) {
        x[,,k][x_nas] <- NA
    }
    x[missing_obs[[1]], missing_obs[[2]], ] <- NA
    x <- aperm(x, c(1L, 3L, 2L))
    
    list(x=x, y=y, mt=mt)
}


#' @importFrom stats sd
#' @export
print.opm <- function(x, digits = max(3, getOption("digits") - 2), prefix = "\t", ...) {
    cat(strwrap(paste("Panel design:", x$design), prefix = prefix), sep = "\n")
    cat("\n")
    cat('Call:\n')
    dput(x$call, control = NULL)
    
    cat('\nCoefficients:\n')

    ci <- apply(confint(x), 1, function(c) {
        paste0('(', paste(format(c, digits = 2, trim = TRUE),
                          collapse = ', '), ')')
    })
    
    coefs <- format(data.frame(`mean (SD)` = colMeans(data.frame(x$samples)),
                               med = coef(x),
                               `95-CI` = ci,
                               check.names = FALSE), digits = digits)

    ## concatenate SD in parentheses to the mean column
    coefs[,1] <- paste(coefs[[1]],
                       paste0('(', round(sapply(data.frame(x$samples), sd), 2), ')'))

    print(format(coefs, digits = digits), print.gap = 2, quote = FALSE)
    
    cat('\n')
    invisible(x)
}


#' @export
summary.opm <- function(object, ...) {
    quants <- t(quantile(object, probs=c(.025, .16, .5, .84, .975),
                         names = FALSE))
    colnames(quants) <- c('<--95CI', '<--68CI', 'med',
                          '68CI-->', '95CI-->')
    structure(list(quants = quants,
                   call = object$call),
              class = 'summary.opm')
}


#' @export
print.summary.opm <- function(x, digits = max(3, getOption("digits") - 2), ...) {
    cat('Call:\n')
    dput(x$call, control = NULL)
    
    cat('\nParameter estimates:\n')
    print(x$quants, digits = digits, print.gap = 3, quote = FALSE, ...)
    invisible(x)
}


#' @importFrom stats coef
#' @export
coef.opm <- function(object, ...) {
    if (any(c('probs', 'names') %in% names(list(...)))) {
        stop("Arguments 'probs' and 'names' are not allowed")
    }
    
    quantile(object, probs = .5, names = FALSE, ...)
}


#' Credible Intervals for Model Parameters
#'
#' Computes equal-tailed credible intervals for one or more parameters
#' in a fitted \code{opm} model. The method used is the quantile
#' interval of the posterior sample.
#'
#' @param object an instance of class \code{opm} whose credible
#'        intervals are wanted
#' @param parm a specification of which parameters are to be given
#'        credible intervals, either a vector of names ("\code{rho}", "\code{sig2}",
#' and "\code{beta}" are the only legal values) or a vector
#'        of positional indices. If missing, all parameters are
#'        considered.
#' @param level the size of the interval (e.g., 0.95 for 95\% C.I.)
#' @param ... additional argument(s) for methods
#' 
#' @return A matrix with columns giving lower and upper limits of the
#' credible interval for each parameter. These will be labeled as (1 -
#' level/2) and 1 - (1 - level)/2 in \% (by default, "2.5\%" and
#' "97.5\%").
#' 
#' @seealso \code{\link{confint}}
#' @importFrom stats confint
#' @export
confint.opm <- function(object, parm, level = 0.95, ...) {
    if (missing(parm)) {
        parm <- names(object$samples)
    } else if (is.numeric(parm)) {
        parm <- names(object$samples)[parm]
    }

    a <- (1 - level)/2

    t(quantile(object, parm = parm, probs = c(a, 1-a)))
}


#' Posterior Sample Quantiles
#'
#' Produces quantiles of the posterior samples corresponding to the
#' given probabilities. In other words, it is equivalent to computing
#' "\code{quantile(x, ...)}", where "\code{x}" is the original Monte
#' Carlo sample of the parameter "\code{parm}", as produced by
#' \code{\link{opm}}.
#' 
#' @param x an instance of class \code{opm} whose sample quantiles are
#'        wanted
#' @param parm a specification of which parameters are to be given
#'        quantiles, either a vector of names ("\code{rho}", "\code{sig2}",
#' and "\code{beta}" are the only legal values) or a vector of positional
#' indices. If missing, all parameters are considered.
#' @param ... further arguments passed to the \code{\link{quantile}}
#'        function operating on the individual parameter's samples
#' 
#' @return A matrix of quantiles for each of the desired parameters,
#' with parameters arranged in columns. If arguments include
#' "\code{names = FALSE}", the quantile labels won't be included
#' (i.e., the rownames of the matrix will be \code{NULL}).
#' 
#' @seealso \code{\link{quantile}}
#' @importFrom stats quantile
#' @export
quantile.opm <- function(x, parm, ...) {
    if (missing(parm)) {
        parm <- names(x$samples)
    } else if (is.numeric(parm)) {
        parm <- names(x$samples)[parm]
    }

    sapply(data.frame(x$samples[parm]),
           quantile, ...)
}


#' @importFrom stats logLik
#' @export
logLik.opm <- function(object, ...) {
    structure(object$logLik,
              nobs = length(object$residuals),
              df = object$df.resid,
              class = 'logLik')
}

#' Deviance Information Criterion (DIC)
#' 
#' Computes the Deviance Information Criterion (DIC), which is a
#' generalization of the Akaike Information Criterion. Models with
#' smaller DIC are considered to fit better than models with larger
#' DIC.
#'
#' DIC is defined as \eqn{DIC = 2*\bar{D} - D_\theta}
#' where:
#' \eqn{\bar{D} = -2 mean(log-likelihood at parameter samples)}
#' \eqn{D_\theta = -2 * log(likelihood at expected value of parameters)}
#'
#' DIC is calculated as: \code{2 * (-2 * mean(log-likelihood at each element of parameter samples)) - (-2 * log(likelihood at mean parameter sample value))}
#'
#' @param object an instance of class \code{opm} whose DIC is wanted.
#' @param ... further arguments passed to other methods.
#' 
#' @return a numeric value with the corresponding DIC
#' 
#' @note Note the speed of computation of the DIC in proportional to
#'     the number of sampled values of the parameters in the opm
#'     object.
#' @export
DIC <- function(object, ...) {
    cl <- object$call
    env <- object$.Environment
    if ('terms' %in% names(object)) {
        data <- if (!is.null(cl$data)) eval(cl$data, env) else env
        index <- if (!is.null(cl$index)) eval(cl$index, env) else eval(formals(opm.formula)$index)
        subset <- if (!is.null(cl$subset)) eval(cl$subset, env) else eval(formals(opm.formula)$subset)

        inputs <- reshape_inputs(cl$x, data, subset, index, env)
        x <- inputs$x
        y <- inputs$y
    }
    else {
        env <- object$.Environment
        x <- eval(cl[[2]], env)
        y <- eval(cl[[3]], env)
    }
    if (isTRUE(object$time.indicators)) {
        x <- with_time_indicators(x)
    }
    
    sample_parameters <- object$samples
    K <- ncol(x)
    D1_at_mean_param <- with(sample_parameters,
                             log_likelihood(x, y, mean(rho), matrix(colMeans(beta), K), 1/mean(sig2)))
    mean_D1 <- mean(with(sample_parameters,
                         mapply(function(r, b, v) log_likelihood(x, y, r, matrix(b, K), v),
                                rho,
                                data.frame(t(beta)),
                                1/sig2)))
    2 * (-2 * mean_D1) - 2 * D1_at_mean_param
}


#' Histogram of an \code{opm} Object
#'
#' Method for \code{\link{hist}} applied to \code{\link{opm}} objects.
#' Each parameter will be plotted in a separate figure.
#'
#' @param x an instance of class \code{opm}
#' @param parm a specification of which parameters are to be plotted,
#'        either a vector of names ("\code{rho}", "\code{sig2}"
#' and "\code{beta}" are the only legal values) or a vector of positional
#' indices. If missing, all parameters are considered.
#' @param ask if "\code{TRUE}", and the R session is interactive, the
#' user is asked to press a key before a new figure (i.e., histogram
#' of the next model parameter) is drawn.
#' @param plot if "\code{TRUE}" (default), the resulting object of
#' class "\code{histogram}" is plotted by \code{plot.hist}.
#' @param main,xlab (optional) vector of titles and X-axis labels
#'        for \emph{each} figure.
#' @param ... further arguments passed to the \code{\link{hist}}
#'        function operating on the individual parameter's samples
#'
#' @return A list of objects of class "\code{histogram}", one for each
#' requested model parameter. The elements are named after the
#' parameter.
#' 
#' @importFrom graphics hist par
#' @importFrom grDevices dev.interactive
#' @importFrom stats setNames
#' @export
hist.opm <- function(x, parm, ask = dev.interactive(), plot = TRUE,
                     main = NULL, xlab = NULL, ...) {
    if (missing(parm)) {
        parm <- names(x$samples)
    } else if (is.numeric(parm)) {
        parm <- names(x$samples)[parm]
    }

    old_par <- NULL
    on.exit(par(old_par))
    
    samples <- data.frame(x$samples[parm])
    sample_names <- names(samples)
    
    if (is.null(main)) {
        main <- lapply(sample_names, function(p) {
            paste('Density of samples of', p)
        })
    } else if (length(main) != length(sample_names)) {
        stop('You need to provide titles for each plotted parameter')
    }
    if (is.null(xlab)) {
        xlab <- sample_names
    } else if (length(xlab) != length(sample_names)) {
        stop('You need to provide X-labels for each plotted parameter')
    }
    
    results <- lapply(seq_along(sample_names), function(i) {
        p <- sample_names[i]
        result <- hist(samples[[p]], plot=plot,
                       main = main[i],
                       xlab = xlab[i], ...)
        if (is.null(old_par)) {
            old_par <<- par(ask = ask)
        }
        result
    })
    invisible(setNames(results, names(samples)))
}


#' Plot Method for an \code{opm} Object
#'
#' Method for \code{\link{plot}} applied to \code{\link{opm}} objects.
#' Each parameter will be plotted as a density plot in a separate
#' figure.
#'
#' @param x an instance of class \code{opm}
#' @param parm a specification of which parameters are to be plotted,
#'        either a vector of names ("\code{rho}", "\code{sig2}"
#' and "\code{beta}" are the only legal values) or a vector of positional
#' indices. If missing, all parameters are considered.
#' @param ask if "\code{TRUE}", and the R session is interactive, the
#' user is asked to press a key before a new figure (i.e., histogram
#' of the next model parameter) is drawn.
#' @param main,xlab (optional) vector of titles and X-axis labels
#'        for \emph{each} figure.
#' @param ... further arguments passed to the \code{\link{plot}}
#'        function operating on the individual parameter's samples
#'
#' @return A list of objects of class "\code{density}", one for each
#' requested model parameter. The elements are named after the
#' parameter.
#' 
#' @importFrom graphics plot par
#' @importFrom grDevices dev.interactive
#' @importFrom stats density setNames
#' @export
plot.opm <- function(x, parm, ask = dev.interactive(),
                     main = NULL, xlab = NULL, ...) {
    if (missing(parm)) {
        parm <- names(x$samples)
    } else if (is.numeric(parm)) {
        parm <- names(x$samples)[parm]
    }

    old_par <- NULL
    on.exit(par(old_par))

    samples <- data.frame(x$samples[parm])
    sample_names <- names(samples)
    
    if (is.null(main)) {
        main <- lapply(sample_names, function(p) {
            paste('Density of samples of', p)
        })
    } else if (length(main) != length(sample_names)) {
        stop('You need to provide titles for each plotted parameter')
    }
    if (is.null(xlab)) {
        xlab <- sample_names
    } else if (length(xlab) != length(sample_names)) {
        stop('You need to provide X-labels for each plotted parameter')
    }
    
    results <- lapply(seq_along(sample_names), function(i) {
        p <- sample_names[i]
        dens <- density(samples[[p]])
        plot(dens, ask = FALSE,
             main = main[i],
             xlab = xlab[i], ...)
        if (is.null(old_par)) {
            old_par <- par(ask = ask)
        }
        dens
    })
    invisible(setNames(results, sample_names))
}


#' Caterpillar Plots of \code{opm} Model Parameters
#'
#' Creates side-by-side plots of equal-tailed credible intervals of \code{opm}
#' model parameters. The intervals are displayed as horizontal lines,
#' with 90\% interval using a thicker line width and 95\% interval a
#' thinner one. The posterior median is indicated with a dot.
#' 
#' @param x an instance of class \code{opm}
#' @param parm a specification of which parameters are to be plotted,
#'        either a vector of names ("\code{rho}", "\code{sig2}"
#' and "\code{beta}" are the only legal values) or a vector of positional
#' indices. If missing, all parameters are considered.
#' @param main,xlab useful defaults for the plot title and X-axis label
#' @param labels labels for each parameter's interval: see \code{\link[graphics]{axis}}
#'
#' @return A matrix of 2.5\%, 5\%, 50\%, 95\%, and 97.5\% quantiles for
#' each of the desired parameters, with parameters arranged in
#' columns.
#'
#' @examples
#' \dontrun{
#' caterplot(o, main = NULL, labels = expression(alpha, beta, sigma^2))
#' }
#' 
#' @importFrom graphics axis segments points
#' @export
caterplot <- function(x, parm, main = paste('Caterpillar plot of', xname),
                      xlab = 'Range of parameter samples',
                      labels = colnames(ranges)) {
    xname <- paste(deparse(substitute(x)), collapse='\n')
    
    ranges <- quantile(x, parm, probs=c(.025, 0.05, .5, .95, .975))
    plot(NULL, xlim = range(ranges), ylim = c(ncol(ranges)+1, 0),
         main = main, xlab = xlab, yaxt='n', ylab = '')
    axis(2, at = seq_len(ncol(ranges)), labels = labels, las = 1)
    segments(ranges[1,], seq_len(ncol(ranges)), ranges[5,], seq_len(ncol(ranges)))
    segments(ranges[2,], seq_len(ncol(ranges)), ranges[4,], seq_len(ncol(ranges)), lwd=3)
    points(ranges[3,], seq_len(ncol(ranges)), pch=19)
    invisible(ranges)
}




#' Long run effects based on the \code{opm} Model Parameters
#'
#' Computes long run effects and confidence intervals of \code{opm} Model Parameters
#'
#' @param opm_obj an instance of class \code{opm}
#' @param parm a specification of which parameters are to be plotted,
#'        a vector of names are the only legal values. If missing, all parameters are considered.
#' @param probs a vector of specified quantiles, by default, the c(0.025,0.5,0.975) are ("\code{probs}")
#' @return A matrix with quantiles on the rows, with number of rows specified as length of the \code{probs} vector for the specified quantiles, with covariates on the columns  
#' @examples
#' \dontrun{
#' longRunEffects(opm_obj)
#' longRunEffects(opm_obj,probs=c(0.975, 0.16, 0.5, 0.84, 0.025))
#' }
#' 
#' @importFrom stats quantile
#' @export
longRunEffects <- function(opm_obj,parm=NULL,probs=c(0.025,0.5,0.975)){
    
    if (is.null(parm)){
        param=opm_obj$samples$beta
        VarNames=colnames(param)
    }else{
        param=opm_obj$samples$beta[,parm]
        VarNames=parm
        
    }
    
    longRunEff=list()
    if (is.vector(param)){
          longRunEff[[1]]=quantile(param/(1-opm_obj$samples$rho), probs=probs) 
          
   }else{
    
        for (i in (1:ncol(param))){
          longRunEff[[i]]=quantile(param[,i]/(1-opm_obj$samples$rho), probs=probs) 
        }
    }
    
    longRunEff_new=do.call(rbind,longRunEff)
    rownames(longRunEff_new)=VarNames
    
    return(t(longRunEff_new))
}


#' Caterpillar Plots of long run effects based on \code{opm} Model Parameters
#'
#' Creates side-by-side plots of equal-tailed credible intervals of \code{opm}
#' the long run effects parameters. The intervals are displayed as horizontal lines,
#' with 90\% interval using a thicker line width and 95\% interval a
#' thinner one. The posterior median is indicated with a dot.
#' 
#' @param x an instance of class \code{opm}
#' @param parm a specification of which parameters are to be plotted,
#'        a vector of names are the only legal values. If missing, all parameters are considered.
#' @param main,xlab useful defaults for the plot title and X-axis label
#' @param probs a vector specifying the quantiles, the defaults is 2.5\%, 5\%, 50\%, 95\%, and 97.5\% quantiles
#' @param labels labels for each parameter's interval: see \code{\link[graphics]{axis}}
#'
#' @return A matrix of 2.5\%, 5\%, 50\%, 95\%, and 97.5\% quantiles for
#' each of the desired parameters, with model parameters arranged in
#' columns.
#'
#' @examples
#' \dontrun{
#' caterplot_longRun(o, main = NULL)
#' }
#' 
#' @importFrom graphics axis segments points
#' @export
caterplot_longRun <- function(x, parm=NULL,main = 'Caterpillar plot of long run effects',
                      xlab = 'Range of parameter samples',probs=c(.025, 0.05, .5, .95, .975),labels = colnames(ranges)) {

   # xname <- paste(deparse(substitute(x)), collapse='\n')
    
    ranges <- longRunEffects(x,parm=parm,probs=probs)
    labels = colnames(ranges) 
    ncr=ncol(ranges)
    nProbs=length(probs)
    indMedian=which(probs==0.5)
    # if (!(is.null(parm))){
    #     ranges=data.frame(ranges[,parm])
    #     labels=parm
    #     ncr = length(parm)
    # }
    
    
   # ranges <- quantile(x, parm, probs=c(.025, 0.05, .5, .95, .975))

    plot(NULL, xlim = range(ranges), ylim = c(ncr+1, 0),
         main = main, xlab = xlab, yaxt='n', ylab = '')
    axis(2, at = seq_len(ncr), labels = labels, las = 1)
    segments(ranges[1,], seq_len(ncr), ranges[nProbs,], seq_len(ncr))
    segments(ranges[2,], seq_len(ncr), ranges[nProbs-1,], seq_len(ncr), lwd=3)
    points(ranges[indMedian,], seq_len(ncr), pch=19)
    invisible(ranges)

}


