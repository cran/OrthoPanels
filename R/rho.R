#' Returns the posterior density \eqn{p(\rho|\Theta)}
#' 
#' @param x an array of dimension \code{time x variable x case} of terms.
#' @param y a matrix of dimension \code{time x case} of responses.
#' @param rho vector of quantiles.
#' @param log.p if \code{TRUE}, probabilities are given as \code{log(p)}.
#' @export
p_rho <- function(x, y, rho, log.p = FALSE) {
    
    T <- nrow(y) - 1                 # initial values are not modelled
    N <- ncol(y)
    K <- dim(x)[2]

    Ti <- Ti(x, y) - 1                  # x and y are not centered yet

    ## measure X_i from agent-specific means for each column
    x <- center_x(x, Ti)
    
    ## likelihood is calculated for (Y_it - Y_i0)
    y <- center_y(y, Ti)

    density <- numeric(length(rho))
    for (i in seq_along(rho)) {
        density[i] <- sum(b(rho[i], Ti)) - ((sum(Ti-1) - K) / 2) *
            log(Q_star(x, w(y, rho[i])))
    
    }

    if (log.p) {
        density
    } else {
        exp(density)
    }
}


## Moves latecomers to still start at T=1
align_latecomers <- function(x, y) {
    T <- nrow(y)
    N <- ncol(y)
    K <- dim(x)[2]
    
    for (i in seq_len(N)) {
        t <- 1L
        while (t <= T &&
               (any(is.na(x[t, , i])) ||
                any(is.na(y[t, i])))) {
            t <- t + 1L
        }
        if (t > 1L && t < T) {
            t_keep <- T-t+1
            x[seq_len(t_keep), , i] <- x[seq(t, len=t_keep), , i]
            y[seq_len(t_keep), i] <- y[seq(t, len=t_keep), i]
            x[seq(t_keep+1, T),, i] <- NA
            y[seq(t_keep+1, T), i] <- NA
        }
    }
    list(x = x, y = y)
}


## Centers X_i from agent-specific means for each column (dropping the first wave)
center_x <- function(x, Ti) {
    T <- nrow(x) - 1                 # initial values are not modelled
    K <- ncol(x)
    N <- dim(x)[3]
    
    ## initial values are not modelled
    x_0 <- matrix(x[1,,], K, N)
    x <- x[-1, , , drop = FALSE]

    ## zero-out values after Ti
    for (i in seq_len(N)) {
        Tix <- seq(to = T, length = T-Ti[i])
        x[Tix, , i] <- 0
    }
    
    ## measure X_i from agent-specific means for each column
    agent_means <- apply(x, 3, colSums) / rep(Ti, each=K)
    agent_means[is.nan(agent_means)] <- 0
    zero_out(sweep(x, 2:3, agent_means), Ti)
}


## Centers Y_i relative to the first wave (dropping it)
center_y <- function(y, Ti) {
    ## initial values are not modelled
    y_0 <- y[1,]
    y <- y[-1,]
    T <- nrow(y)
    
    ## zero-out values after Ti
    y <- zero_out(y, Ti)
    
    ## likelihood is calculated for (Y_it - Y_i0)
    zero_out(sweep(y, 2, y_0), Ti)
}


## Centers Y_i from agent-specific mean after dropping the drop_row wave
center_yt <- function(y, Ti, drop_row = 1) {
    y <- y[-drop_row,]
    T <- nrow(y)
    
    # zero-out values after Ti
    y <- zero_out(y, Ti)
    
    agent_means <- colSums(y) / Ti
    ## now subtract agent-specific means
    y <- sweep(y, 2, agent_means)

    # zero-out values after Ti
    zero_out(y, Ti)
}

## Zero-out all values of `x` for each case `i` after `Ti[i]`-th time step.
##
## Assumes T is the first dimension and cases are the last dimension,
## but allows multi-dimensional values for each case+timestep
## combination
zero_out <- function(x, Ti, val = 0) {
    N <- length(Ti)
    T <- nrow(x)
    
    dims <- lapply(dim(x), seq_len)
    
    ## zero-out values after Ti in the `margin`
    for (i in seq_len(N)) {
        Tix <- seq(to = T, length = T-Ti[i])
        dims[[1]] <- Tix
        dims[[length(dims)]] <- i
        ix <- as.matrix(do.call('expand.grid', dims))
        x[ix] <- val
    }
    x
}

## Returns $w_i = y_i - \rho y_{i-}
w <- function(y, rho) {
    T <- nrow(y)
    N <- ncol(y)

    w <- matrix(0, T, N)
    for (i in seq_len(N)) {
        w[,i] <- y[,i] - rho * y_(y, i)
    }
    w
}


## Returns $y_{i-} = ( y_{i0}, y_{i1}, \dots, y_{i, T-1} ) $
y_ <- function(y, i) {
    T <- nrow(y)
    
    c(0, y[-T, i])
}

## Returns $\frac{1}{T} \sum_{t=1}^{T-1} \frac{T-t}{t} \rho^t$
b <- function(rho, Ti) {
    sapply(Ti, function(T) {
        if (T < 1) return(0)
        t <- seq_len(T-1)
        sum((T - t) * rho^t / t) / T
    })
}


## Returns $\sum_i w'_i H w_i - \sum_i w'_i H x_i \left(
##    \sum_i x'_i H x_i \right)^{-1} \sum_i x'_i H w_i$
Q_star <- function(x, w) {
    wHw(w) - wHx(x, w) %*% solve(Hstar(x)) %*% xHw(x, w)
}


## Returns $H = I_T - (1/T)_{J J'}$
h <- function(n, Ti=n) {
    m <- matrix(0, n, n)
    if (Ti < 1) return(m)
    m[seq_len(Ti), seq_len(Ti)] <- diag(Ti) - matrix(1/Ti, Ti, Ti)
    m
}


## Returns $\sum_i w_i' H w_i$
wHw <- function(w, Ti=rep(T, N), H=lapply(Ti, h, n=T)) {
    T <- nrow(w)
    N <- ncol(w)

    sum <- 0
    for (i in seq_len(N)) {
        sum <- sum + t(w[,i]) %*% H[[i]] %*% w[,i]
    }
    sum
}


## Returns $\sum_i w_i' H x_i$
wHx <- function(x, w, Ti=rep(T, N), H=lapply(Ti, h, n=T)) {
    T <- nrow(w)
    N <- ncol(w)

    sum <- 0
    for (i in seq_len(N)) {
        sum <- sum + t(w[,i]) %*% H[[i]] %*% x[,,i]
    }
    sum
}


## Returns $\sum_i x_i' H w_i$
xHw <- function(x, w, Ti=rep(T, N), H=lapply(Ti, h, n=T)) {
    T <- nrow(w)
    N <- ncol(w)

    sum <- 0
    for (i in seq_len(N)) {
        sum <- sum + t(x[,,i]) %*% H[[i]] %*% w[,i]
    }
    sum
}


## Returns $(\sum_i x_i' H x_i)^{-1}$
Hstar <- function(x, Ti=rep(T, N), H=lapply(Ti, h, n=T)) {
    T <- dim(x)[1]
    N <- dim(x)[3]

    sum <- 0
    for (i in seq_len(N)) {
        sum <- sum + t(x[,,i]) %*% H[[i]] %*% x[,,i]
    }
    sum
}


Ti <- function(x, y) {
    T <- dim(x)[1]
    N <- dim(x)[3]
    
    Ti <- integer(N)

    x <- is.na(x[-1,,, drop = FALSE])
    for (i in seq_len(N)) {
        Xi <- x[,,i, drop=FALSE]
        
        no_response <- which(apply(Xi, 1, any))
        
        Ti[i] <- if (length(no_response)>0) min(no_response) else T
    }
    no_response_y <- apply(is.na(y), 2,
                           function(Y_i) if (any(Y_i)) which.max(Y_i)-1 else T)

    pmin(Ti, no_response_y)
}
