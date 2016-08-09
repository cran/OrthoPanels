## Draw `n` samples of rho, v, and beta parameters from the
## distribution defined by the data `x` and `y`, with possible values
## of rho defined by the set `pts`
##' @importFrom MASS mvrnorm
##' @importFrom stats rgamma
sample_all <- function(x, y, n, pts = seq(-.995, .995, .001)) {
    T <- nrow(y) - 1                 # initial values are not modelled
    N <- ncol(y)
    K <- dim(x)[2]

    Ti <- Ti(x, y) - 1                  # x and y are not centered yet

    H <- lapply(Ti, h, n=T)

    ## measure X_i from agent-specific means for each column
    x <- center_x(x, Ti)
    
    ## likelihood is calculated for (Y_it - Y_i0)
    y <- center_y(y, Ti)

    density <- numeric(length(pts))
    Q_star <- numeric(length(pts))
    H_ast <- Hstar(x, Ti, H)
    beta_hat <- list(length(pts))
    for (i in seq_along(pts)) {
        W <- w(y, pts[i])
        WHW <- wHw(W, Ti, H)
        WHX <- wHx(x, W, Ti, H)
        XHW <- xHw(x, W, Ti, H)

        beta_hat[[i]] <- solve(H_ast, XHW)
        
        Q_star[i] <- WHW - WHX %*% beta_hat[[i]]

        ## prevent getting infinite/NaN density
        sm <- abs(Q_star[i] - .Machine$double.eps)
        if (is.na(sm) || sm < .Machine$double.eps) Q_star[i] <- .Machine$double.eps
        
        density[i] <- sum(b(pts[i], Ti)) - ((sum(ifelse(Ti >= 1, Ti-1, 0)) - K) / 2) *
            log(Q_star[i])
    }

    density <- exp(density - max(density))
    is <- sample(seq_along(pts), n, replace = TRUE, prob = density / sum(density))
    rhos <- pts[is]
    vs <- numeric(n)
    betas <- matrix(, K, n)
    if (K > 1) {
        rownames(betas) <- colnames(x)
    }
    H_ast_inv <- solve(H_ast)
    for (i in seq(n)) {
        j <- is[i]

        vs[i] <- rgamma(1, shape = (sum(ifelse(Ti >= 1, Ti-1, 0)) - K)/2,
                        rate = Q_star[j]/2)
        betas[, i] <- mvrnorm(1, beta_hat[[j]], H_ast_inv/vs[i])
    }
    
    list(rho = rhos,
         sig2 = 1/vs,
         beta = t(betas))
}


## Calculates the fitted values yf[t] = rho * yc[t-1] + beta * xc[t]
## where xc and yc are x and y, respectively, centered on
## agent-specific means after the first row (t=0) has been dropped
fitted <- function(x, y, rho, beta) {
    T <- nrow(y) - 1                 # initial values are not modelled
    N <- ncol(y)
    K <- dim(x)[2]

    Ti <- Ti(x, y) - 1                  # x and y are not centered yet

    ## measure X_i from agent-specific means for each column
    cx <- center_x(x, Ti)
    
    ## drop first row and subtract from agent-specific means for each column
    cy <- center_yt(y, Ti)

    ## drop the *last* row and subtract from agent-specific means for
    ## each column
    cy_t1 <- center_yt(y, Ti, nrow(y))

    y_f <- matrix(0, T, N)
    
    for (t in seq_len(T)) {
        y_f[t, ] <- rho * cy_t1[t, ] + crossprod(beta, cx[t,,])
    }
    
    y_f
}


## Calculates the log-likelihoood
log_likelihood <- function(x, y, rho, beta, v) {
    T <- nrow(y) - 1                 # initial values are not modelled
    N <- dim(x)[3]

    Ti <- Ti(x, y) - 1                  # x and y are not centered yet

    H <- lapply(Ti, h, n=T)

    ## measure X_i from agent-specific means for each column
    x <- center_x(x, Ti)

    ## likelihood is calculated for (Y_it - Y_i0)
    y <- center_y(y, Ti)

    W <- w(y, rho)

    Q <- 0
    for (i in seq_len(N)) {
        z <- W[,i] - x[,,i] %*% beta
        Q <- Q + t(z) %*% H[[i]] %*% z
    }
    c(sum(Ti)/2 * log(v) + sum(b(rho, Ti)) - Q*v/2)
}


## Draw `n` samples from `rho`s, proportional to the posterior pdf(rho | x, y)
sample_rho <- function(n, x, y, rho) {
    ps <- p_rho(x, y, rho, log.p = TRUE)
    ps <- exp(ps - max(ps))
    sample(rho, n, replace = TRUE, prob = ps / sum(ps))
}


## Sample variance from the gamma distribution with parameters defined
## by the data `x` and `y` and given a sample of `rho`s
##' @importFrom stats rgamma
sample_sig <- function(x, y, rho) {
    T <- nrow(y) - 1
    N <- ncol(y)
    K <- dim(x)[2]

    Ti <- Ti(x, y) - 1                  # x and y are not centered yet
    x <- center_x(x, Ti)
    y <- center_y(y, Ti)
    
    sig <- numeric(length(rho))
    rho_rle <- rle(sort(rho))
    j <- 1
    for (i in seq_along(rho_rle$values)) {
        n <- rho_rle$lengths[i]         # number of draws
        sig[j:(j+n-1)] <- rgamma(n, shape = (N*(T-1) - K)/2,
                                 rate = Q_star(x, w(y, rho_rle$values[i]))/2)
        j <- j + n
    }
    sig[order(order(rho))]
}


## Sample beta from the normal distribution with parameters defined by
## the data `x` and `y` and given a sample of `rho`s and variance `v`s
##' @importFrom MASS mvrnorm
sample_beta <- function(x, y, rho, v) {
    stopifnot(length(rho) == length(v))
    
    Ti <- Ti(x, y) - 1                  # x and y are not centered yet
    x <- center_x(x, Ti)
    y <- center_y(y, Ti)

    K <- dim(x)[2]
    
    h_ast_inv <- solve(Hstar(x))
    beta <- matrix(, K, length(rho))
    
    rho_rle <- rle(sort(rho))
    irho <- order(rho)
    j <- 1
    
    for (i in seq_along(rho_rle$values)) {
        n <- rho_rle$lengths[i]         # number of draws
        beta_hat <- h_ast_inv %*% xHw(x, w(y, rho_rle$values[i]))
        
        for (k in 0:(n-1)) {
            vv <- v[irho[j+k]]
            beta[ , j+k] <- mvrnorm(1, beta_hat, h_ast_inv/vv)
        }
        j <- j + n
    }
    
    t(beta[ , order(irho), drop = FALSE])
}
