# Functions for constructing covariance with exchangable covariance

exchangable_cov <- function(main_var, re_var, block_sizes) {
  accum <- c(0, cumsum(block_sizes))
  mat <- matrix(0, max(accum), max(accum))
  mat_inv <- matrix(0, max(accum), max(accum))
  diag(mat) <- main_var
  diag(mat_inv) <- 1 / main_var
  for (i in seq_along(block_sizes)) {
    if (block_sizes[i] == 0) next
    mat[(accum[i] + 1):(accum[i + 1]), (accum[i] + 1):(accum[i + 1])] <-
      mat[(accum[i] + 1):(accum[i + 1]), (accum[i] + 1):(accum[i + 1])] + re_var
    mat_inv[(accum[i] + 1):(accum[i + 1]), (accum[i] + 1):(accum[i + 1])] <-
      mat_inv[(accum[i] + 1):(accum[i + 1]), (accum[i] + 1):(accum[i + 1])] -
      re_var / main_var / (main_var + re_var * block_sizes[i])
  }
  return(list(V = mat, Prec = mat_inv))
}
# Supremum difference between two functions
Linf_norm <- function(f1, f2, L = -Inf, R = +Inf) {
  return(optim(0, function(x) abs(f1(x) - f2(x)), lower = L, upper = R, control = list(fnscale = -1)))$value
}

# General package loader; install package if not already installed
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[, 1])) {
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}
# Functions related to umbrella constraints
is.increasing <- function(x) {
  if (length(x) == 1) {
    return(TRUE)
  } else {
    return(all(diff(x) >= 0))
  }
}
is.decreasing <- function(x) {
  if (length(x) == 1) {
    return(TRUE)
  } else {
    return(all(diff(x) <= 0))
  }
}
is.umbrella <- function(x) {
  ind_max <- which.max(x)
  return(is.increasing(x[1:ind_max]) && is.decreasing(x[ind_max:length(x)]))
}
# Sample from umbrella-constrained MVN through rejection sampling
rtMVN <- function(mean, Sigma, posit = 1:length(mean),
                  SShape = FALSE) {
  if (is.null(posit)) {
    return(as.vector(mvtnorm::rmvnorm(1, mean, Sigma)))
  }
  if (!(SShape)) {
    while (TRUE) {
      x <- as.vector(mvtnorm::rmvnorm(1, mean, Sigma))
      if (all(x[posit] >= 0)) {
        return(x)
      }
    }
  } else {
    while (TRUE) {
      x <- as.vector(mvtnorm::rmvnorm(1, mean, Sigma))
      if (is.umbrella(x[posit]) && all(x[posit] >= 0)) {
        return(x)
      }
    }
  }
}
# Generate smooth/empirical quantiles for (0,1) data
betaKDE <- function(z, s, q) {
  # If s is NULL, use empirical quantiles
  # If s is between 0 and 1, use BETA-KDE
  # Otherwise, terminate the program with error
  if (is.null(s)) {
    quan <- quantile(z, q)
    names(quan) <- paste(round(q * 100, 1), "%", sep = "")
    return(list(density = NULL, CDF = NULL, quantile = quan))
  }
  if ((s <= 0) || (s >= 1)) {
    stop("Variance inflation factor should be NULL or between 0 and 1")
  }
  a <- 1 / s
  dens <- function(x) {
    res <- 0
    for (t in z) {
      res <- res + dbeta(x, a * t + 1, a * (1 - t) + 1)
    }
    res <- res / length(z)
    return(res)
  }
  CDF <- function(x) {
    res <- 0
    for (t in z) {
      res <- res + pbeta(x, a * t + 1, a * (1 - t) + 1)
    }
    res <- res / length(z)
    return(res)
  }

  quan <- rep(0, length(q))
  for (i in 1:length(q)) {
    quan[i] <- uniroot(
      function(x) {
        CDF(x) - q[i]
      },
      c(0, 1)
    )$root
  }
  names(quan) <- paste(round(q * 100, 1), "%", sep = "")
  return(list(density = dens, CDF = CDF, quantile = quan))
}
# Create Planck-Taper window
planck_taper <- function(dimension, eps = 0.1) {
  N <- dimension - 1
  w <- rep(0, dimension)
  w[1] <- 0
  if (eps * N >= 1) {
    for (i in seq(1, eps * N, by = 1)) {
      w[i + 1] <- (1 + exp(eps * N / i - eps * N / (eps * N - i)))^(-1)
    }
  }
  for (i in seq(ceiling(eps * N), N / 2)) w[i + 1] <- 1
  for (i in seq(0, N / 2)) w[N - i + 1] <- w[i + 1]
  return(w)
}
# Create random-walk + window prior matrix
penalty_Matrix <- function(dimension, smooth.sigma = Inf, flat.sigma = Inf, weight = NULL) {
  if (is.null(weight)) {
    weight <- (0:(dimension - 1)) * ((dimension - 1):0)
  }
  weight <- weight + 0.001
  weight <- weight / max(weight)
  prec <- matrix(0, dimension, dimension)
  for (j in 1:(dimension - 1)) {
    u <- rep(0, dimension)
    u[j] <- 1
    u[j + 1] <- -1
    prec <- prec + u %*% t(u) / (smooth.sigma / dimension / dimension)
  }
  prec <- prec + diag(1 / flat.sigma / weight)
  V <- chol2inv(chol(prec))
  return(list(V = V, prec = prec))
}
# Create block diagonal matrices
block_Matrix <- function(M1, M2) {
  rbind(
    cbind(M1, matrix(0, nrow = nrow(M1), ncol = ncol(M2))),
    cbind(matrix(0, nrow = nrow(M2), ncol = ncol(M1)), M2)
  )
}
# Generate linear constraint matrix for umbrella region
lincon <- function(l, q) {
  k <- 1:(l - q)
  res <- array(0, dim = c(l + 1, l, l - q)) # c(l,l,l-q))
  for (k0 in k) {
    for (i in 1:(q + 1)) res[i, i, k0] <- 1
    if (q + 2 <= l) {
      for (i in (q + 2):l) {
        res[i, i, k0] <- 1
        res[i, i - 1, k0] <- -1
        if (i > q + k0) res[i, , k0] <- res[i, , k0] * (-1)
      }
    }
    res[l + 1, l, k0] <- 1
  }
  return(res)
}


# Sample from Umbrella-restricted MVN (UMVN) with known inflection index
# using HDTG
ruMVN <- function(n, mu, V, q, k, Ms, burnin = 100) {
  # First q parameter has no restriction
  # The rest length(mu)-q parameters are all positive and peaks at its k-th location
  require(hdtg)
  stopifnot(c(length(mu) == dim(V)[1], dim(V)[1] == dim(V)[2], q + k <= length(mu)))
  l <- length(mu)
  M <- Ms[(q + 1):(l + 1), , k]
  init0 <- c(runif(q, -0.5, 0.5))
  tt <- c(rnorm(1, sd = 0.1), runif(k - 1, 0, 0.1), runif(l - q - k, -0.1, 0))
  init0 <- c(init0, exp(cumsum(tt)))
  res <- hdtg::harmonicHMC(n, burnin, mu, cholesky(V),
    M, rep(0, l - q + 1), init0,
    precFlg = FALSE
  )
  return(res)
}
# Sample from UMVN with UNKNOWN inflection index using HDTG
hdtg_S <- function(n, mu, sigma, free = NULL, burnin = 5) {
  p <- length(mu)
  if (is.null(free)) {
    res_index <- 1:p
  } else {
    res_index <- (1:p)[-free]
  }
  p_res <- length(res_index)
  FF <- matrix(0, nrow = p_res + 1, ncol = p)
  for (i in 1:p_res) {
    FF[, res_index][i, i] <- 1
    if (i > 1) FF[, res_index][i, i - 1] <- -1
  }
  FF[, res_index][p_res + 1, p_res] <- 1
  logprob <- rep(0, p_res)
  for (i in 1:p_res) {
    Fmat <- FF
    if (i < p_res) Fmat[(i + 1):p_res, ] <- -Fmat[(i + 1):p_res, ]
    g <- rep(0, p_res + 1)
    new.S <- Fmat %*% sigma %*% t(Fmat)
    new.S[p_res + 1, p_res + 1] <- new.S[p_res + 1, p_res + 1] * (1 + 1e-3)
    new.mu <- Fmat %*% mu
    logprob[i] <- logpmvnorm(rep(0, p_res + 1), rep(+Inf, p_res + 1), new.mu, new.S)
  }
  logprob <- logprob - max(logprob)
  inflex.points <- sample(1:p_res, n, replace = TRUE, prob = exp(logprob))
  result <- matrix(0, nrow = n, ncol = p)
  for (i in 1:p_res) {
    pnums <- sum(inflex.points == i)
    if (pnums == 0) next
    Fmat <- FF
    if (i < p_res) Fmat[(i + 1):p_res, ] <- -Fmat[(i + 1):p_res, ]
    g <- rep(0, p_res + 1)
    init <- rep(0, p)
    init[res_index][i] <- 0.1
    init[res_index][-i] <- 0.1 * runif(p_res - 1)
    if (i > 1) init[res_index][1:(i - 1)] <- sort(init[res_index][1:(i - 1)], decreasing = FALSE)
    if (i < p_res) init[res_index][(i + 1):p_res] <- sort(init[res_index][(i + 1):p_res], decreasing = TRUE)
    result[inflex.points == i, ] <- hdtg::harmonicHMC(
      n = pnums, burnin = burnin,
      mean = mu, choleskyFactor = chol(sigma),
      F = Fmat, g = g, init = init, precFlg = FALSE
    )
  }
  if (n == 1) {
    return(as.vector(result))
  } else {
    return(result)
  }
}
# Calculate log of normalizing constant of UMVN
puMVN <- function(mu, V, q, Ms, k = NULL, eps = 1e-2, log = T, verbose = FALSE) {
  stopifnot(c(length(mu) == dim(V)[1], dim(V)[1] == dim(V)[2], q <= length(mu)))
  require(mvtnorm)
  l <- length(mu)
  if (is.null(k)) k <- 1:(l - q)
  res <- rep(0, length(k))
  for (ki in 1:length(k)) {
    k0 <- k[ki]
    M <- Ms[, , k0]
    new.mu <- drop(M %*% mu)
    new.V <- M %*% V %*% t(M)
    new.V[l + 1, l + 1] <- new.V[l + 1, l + 1] * (1 + eps)
    # print(new.V[l+1,l+1]-eps)
    lb <- c(rep(-Inf, q), rep(0, l - q + 1))
    ub <- rep(+Inf, l + 1)
    pp <- logpmvnorm(lb, ub, new.mu, new.V)
    res[ki] <- pp
  }
  if (log) {
    return((res))
  } else {
    return(exp(res))
  }
}

# Functions for Updating Parameters
# Update BETA and GAMMA together
# RE is the replicated random effects, and should always be the same dimension as Y
update_coef <- function(covars.list, nX, Y, RE, sy, sw, id, prior.mean, prior.precision, Ms, verbose = FALSE,
                        samples = 1, burnin = 5) {
  require(matrixStats)
  res <- array(0, c(samples, ncol(covars.list[[1]]), ncol(Y)))
  unique_id <- unique(id)
  for (k in 1:ncol(Y)) {
    non_mis <- !is.na(Y[, k])
    CC <- covars.list[[k]][non_mis, ]
    comp_id <- id[non_mis]
    block_size <- integer(length(unique_id))
    for (j in 1:length(block_size)) {
      block_size[j] <- sum(comp_id == j)
    }
    lik_prec <- exchangable_cov(sy, sw, block_size)$Prec
    precision <- prior.precision + t(CC) %*% lik_prec %*% CC
    variance <- solve(precision)

    mu <- as.vector(precision %*% prior.mean)
    mu <- mu + t(CC) %*% lik_prec %*% Y[, k][non_mis]

    mu <- as.vector(variance %*% mu)
    if (nX > 0) {
      free_indice <- 1:nX
    } else {
      free_indice <- NULL
    }
    res[, , k] <- hdtg_S(samples, mu, variance, free_indice, burnin)
  }
  return(list(res = res))
}
# Update SIGMA_Y----
update_sigmay <- function(covars.list, Y, RE, coefs, prior.shape, prior.scale) {
  residual <- array(0, dim(Y))
  for (k in 1:ncol(Y)) {
    residual[, k] <- Y[, k] - RE[, k] - drop(covars.list[[k]] %*% coefs[, k])
  }
  return(1 / rgamma(1,
    shape = prior.shape + sum(!is.na(Y)) / 2,
    rate = prior.scale + sum(residual^2, na.rm = TRUE) / 2
  ))
}
# Update PENALTIES----
# Old version: NORMAL transition kernel
update_pens0 <- function(gamma, mu, pen, lpd, ls, weight) {
  p <- dim(gamma)[1]
  V <- penalty_Matrix(p, pen[1], pen[2], weight)$V
  ll <- sum(dtmvnorm(t(gamma), mu, V, lb = rep(0, p), log = TRUE)) +
    lpd(pen)
  pen_pro <- pen + rnorm(2, 0, exp(ls))
  V2 <- penalty_Matrix(p, pen_pro[1], pen_pro[2], weight)$V
  if (pen_pro[1] > 0 & pen_pro[2] > 0) {
    ll2 <- sum(dtmvnorm(t(gamma), mu, V2, lb = rep(0, p), log = TRUE)) +
      lpd(pen_pro)
  } else {
    ll2 <- -Inf
  }
  diff <- ll - ll2
  temp <- rexp(1)
  if (temp > diff) {
    acc_status <- 1
    new <- pen_pro
  } else {
    acc_status <- 0
    new <- pen
  }
  return(list(
    acc_status = acc_status,
    new = new
  ))
}
# New Version: Log-normal transition kernel, recommended
update_pens <- function(
    gamma, # Matrix containing K gamma coefficients
    mu, # Prior Mean of Gamma
    lambda, # Current value of penalty parameters, 1 for smooth, 2 for flat
    lpd, # Log posterior density function of lambda
    ls, # Log of standard deviation of proposal distribution
    weight, # Window Function; NULL means default quadratic window,
    Ms,
    verbose) {
  require(mvtnorm)
  p <- dim(gamma)[1]
  V <- penalty_Matrix(p, lambda[1], lambda[2], weight)$V # Current variance
  # Log-density at current value
  ll <- sum(dmvnorm(t(gamma), mu, V, log = TRUE) - logSumExp(puMVN(mu, V, 0, Ms, verbose = verbose))) + lpd(lambda)
  # Propose new state
  # Since lambda is positive, we use log-normal jumps
  new <- exp(log(lambda) + rnorm(2, sd = exp(ls)))

  V2 <- penalty_Matrix(p, new[1], new[2], weight)$V # New variance
  # Log-density at new value
  ll2 <- sum(dmvnorm(t(gamma), mu, V2, log = TRUE) - logSumExp(puMVN(mu, V2, 0, Ms, verbose = verbose))) + lpd(new)

  diff <- ll - sum(log(new)) - ll2 + sum(log(lambda))
  # diff <- ll - ll2
  temp <- rexp(1)
  if (temp > diff) {
    acc_status <- 1
    state <- new
  } else {
    acc_status <- 0
    state <- lambda
  }
  return(list(
    acc_status = acc_status,
    new = state
  ))
}

# Update Random Effect VARIANCE----
# W is the non-replicated version of RE
# W should contain N*K elements, where N is the number of individuals
# and K is the number of biomarkers
update_sigmaw <- function(W, prior.shape, prior.scale) {
  post.shape <- prior.shape + length(W) / 2
  post.scale <- prior.scale + sum(W^2) / 2
  return(1 / rgamma(1,
    shape = post.shape,
    rate = post.scale
  ))
}
# Update Individual Random Intercept----
update_W <- function(covars.list, Y, coefs, long_ss, ID,
                     sigmay, sigmaw) {
  residual <- array(0, dim(Y))
  for (k in 1:ncol(Y)) {
    residual[, k] <- Y[, k] - drop(covars.list[[k]] %*% coefs[, k])
  }
  sds <- (1 / sigmaw + long_ss / sigmay)^(-1 / 2)
  means <- matrix(0, nrow = nrow(long_ss), ncol = ncol(long_ss))
  for (k in 1:ncol(long_ss)) {
    means[, k] <- tapply(residual[, k], ID, sum, na.rm = TRUE)
  }
  means <- means / (sigmay / sigmaw + long_ss)
  norms <- array(rnorm(prod(dim(long_ss))),
    dim = dim(long_ss)
  )
  return(norms * sds + means)
}
# Calculate w*a+(1-w)*b in log scale; a and b are provided in log scale
lwmean <- function(a, b, w = 1 / 2) {
  require(matrixStats)
  # require(Rmpfr)
  a0 <- a + log(w)
  b0 <- b + log(1 - w)
  mat <- rbind(a0, b0)
  return(colLogSumExps(mat))
}
# Adaptation of PMVNORM function to calculate log probability
# Uses method from Genz et al
logpmvnorm <- function(lb, ub, mu, Sigma, Nmax = 1e3) {
  require(matrixStats)
  require(VGAM)
  precision <- 100
  a <- lb - mu
  b <- ub - mu
  C <- t(chol(Sigma))
  m <- ncol(C)
  d <- array(0, dim = c(Nmax, m))
  e <- array(0, dim = c(Nmax, m))
  f <- array(0, dim = c(Nmax, m))
  y <- array(0, dim = c(Nmax, m - 1))
  temp <- array(0, dim = c(Nmax, m - 1))
  d[, 1] <- ifelse(a[1] == -Inf, -Inf, pnorm(a[1] / C[1, 1], log = T))
  e[, 1] <- ifelse(b[1] == +Inf, 0, pnorm(b[1] / C[1, 1], log = T))
  f[, 1] <- e[1, 1] + log1mexp(e[1, 1] - d[1, 1])

  w <- array(runif((m - 1) * Nmax), dim = c(Nmax, m - 1))
  for (i in 2:m) {
    temp[, i - 1] <- lwmean(d[, i - 1], e[, i - 1], w[, i - 1])
    temp[temp[, i - 1] >= 0, i - 1] <- 0
    y[, i - 1] <- qnorm(temp[, i - 1], log = T)
    y[y[, i - 1] == +Inf, i - 1] <- 1e3
    y[y[, i - 1] == -Inf, i - 1] <- -1e3
    if (a[i] == -Inf) {
      d[, i] <- -Inf
    } else {
      d[, i] <- pnorm((a[i] - colSums(C[i, 1:(i - 1)] * t(y[, 1:(i - 1)]))) / C[i, i], log = T)
    }
    if (b[i] == +Inf) {
      e[, i] <- 0
    } else {
      e[, i] <- pnorm((b[i] - colSums(C[i, 1:(i - 1)] * t(y[, 1:(i - 1)]))) / C[i, i], log = T)
    }
    diff <- log1mexp(e[, i] - d[, i])
    diff <- ifelse(diff != -Inf, diff, log1mexp(.Machine$double.xmin))
    f[, i] <- e[, i] + diff + f[, i - 1]
  }
  res <- logSumExp(f[, m], na.rm = TRUE) - log(sum(!is.na(f[, m])))
  if (is.na(res) || res == -Inf) {
    res <- (-Inf)
    View(f)
    View(d)
    View(e)
    View(y)
    View(temp)
    stop("Error")
    attr(res, "rel_error") <- 0
    return(res)
  } else {
    v <- cbind(f[, m], res)
    idx <- v[, 1] < v[, 2]
    idx[is.na(idx)] <- FALSE
    v[idx, ] <- v[idx, c(2, 1)]
    ad <- v[, 1] + log1mexp(v[, 1] - v[, 2])
    ad <- ad * 2
    var0 <- logSumExp(ad, na.rm = TRUE) - log(sum(!is.na(ad)))
    var0 <- var0 - log(sum(!is.na(ad)))
    sd0 <- var0 / 2
    attr(res, "rel_error") <- exp(sd0 - res)
  }
  return(res)
}
# Sample coefficients with consideration of group structures
update_coef_grouped <- function(covars.list, nX, Y, RE, sy, sw, id, prior.mean, prior.precision, Ms, verbose = FALSE, group = NULL,
                                samples = 1, burnin = 5) {
  require(matrixStats)
  res <- array(0, c(samples, ncol(covars.list[[1]]), ncol(Y)))
  unique_id <- unique(id)
  mu_list <- list()
  var_list <- list()
  for (k in 1:ncol(Y)) {
    non_mis <- !is.na(Y[, k])
    CC <- covars.list[[k]][non_mis, ]
    comp_id <- id[non_mis]
    block_size <- integer(length(unique_id))
    for (j in 1:length(block_size)) {
      block_size[j] <- sum(comp_id == j)
    }
    lik_prec <- exchangable_cov(sy, sw, block_size)$Prec
    precision <- prior.precision + t(CC) %*% lik_prec %*% CC
    variance <- solve(precision)

    mu <- as.vector(precision %*% prior.mean)
    mu <- mu + t(CC) %*% lik_prec %*% Y[, k][non_mis]

    mu <- as.vector(variance %*% mu)
    mu_list[[k]] <- mu
    var_list[[k]] <- variance
  }
  if (nX > 0) {
    free_indice <- 1:nX
  } else {
    free_indice <- NULL
  }
  res_list <- hdtg_S_grouped(samples, mu_list, var_list, free_indice, group, burnin)
  for (k in 1:ncol(Y)) {
    res[, , k] <- res_list[[k]]
  }
  return(list(res = res))
}
# Sample multiple UMVNs with grouped inflection indices
hdtg_S_grouped <- function(n, mu_list, sigma_list, free = NULL, group = NULL, burnin = 5) {
  p <- length(mu_list[[1]])
  K <- length(mu_list)
  if (is.null(group)) {
    group <- 1:K
  } else {
    group <- match(group, unique(group))
  }
  if (is.null(free)) {
    res_index <- 1:p
  } else {
    res_index <- (1:p)[-free]
  }
  p_res <- length(res_index)
  FF <- matrix(0, nrow = p_res + 1, ncol = p)
  for (i in 1:p_res) {
    FF[, res_index][i, i] <- 1
    if (i > 1) FF[, res_index][i, i - 1] <- -1
  }
  FF[, res_index][p_res + 1, p_res] <- 1
  logprob <- list()
  for (k0 in 1:K) logprob[[k0]] <- rep(0, p_res)
  for (k0 in 1:K) {
    for (i in 1:p_res) {
      Fmat <- FF
      if (i < p_res) Fmat[(i + 1):p_res, ] <- -Fmat[(i + 1):p_res, ]
      g <- rep(0, p_res + 1)
      new.S <- Fmat %*% sigma_list[[k0]] %*% t(Fmat)
      new.S[p_res + 1, p_res + 1] <- new.S[p_res + 1, p_res + 1] * (1 + 1e-3)
      new.mu <- Fmat %*% mu_list[[k0]]
      logprob[[k0]][i] <- logpmvnorm(rep(0, p_res + 1), rep(+Inf, p_res + 1), new.mu, new.S)
    }
    logprob[[k0]] <- logprob[[k0]] - max(logprob[[k0]])
  }
  joint_logprob <- list()
  for (g0 in 1:length(unique(group))) {
    joint_logprob[[g0]] <- rep(0, p_res)
  }
  for (k0 in 1:K) {
    joint_logprob[[group[k0]]] <-
      joint_logprob[[group[k0]]] + logprob[[k0]]
  }
  inflex.points <- list()
  for (g0 in 1:length(unique(group))) {
    inflex.points[[g0]] <- sample(1:p_res, n, replace = TRUE, prob = exp(joint_logprob[[g0]]))
  }
  result <- list()
  for (k0 in 1:K) {
    result[[k0]] <- matrix(0, nrow = n, ncol = p)
  }
  for (k0 in 1:K) {
    for (i in 1:p_res) {
      pnums <- sum(inflex.points[[group[k0]]] == i)
      if (pnums == 0) next
      Fmat <- FF
      if (i < p_res) Fmat[(i + 1):p_res, ] <- -Fmat[(i + 1):p_res, ]
      g <- rep(0, p_res + 1)
      init <- rep(0, p)
      init[res_index][i] <- 0.1
      init[res_index][-i] <- 0.1 * runif(p_res - 1)
      if (i > 1) init[res_index][1:(i - 1)] <- sort(init[res_index][1:(i - 1)], decreasing = FALSE)
      if (i < p_res) init[res_index][(i + 1):p_res] <- sort(init[res_index][(i + 1):p_res], decreasing = TRUE)
      result[[k0]][inflex.points[[group[k0]]] == i, ] <- hdtg::harmonicHMC(
        n = pnums, burnin = burnin,
        mean = mu_list[[k0]], choleskyFactor = chol(sigma_list[[k0]]),
        F = Fmat, g = g, init = init, precFlg = FALSE
      )
    }
  }
  if (n == 1) {
    return(lapply(result, as.vector))
  } else {
    return(result)
  }
}
# Sample penalty parameters in model with grouped inflection indices
update_pens_grouped <- function(
    gamma, # Matrix containing K gamma coefficients
    mu, # Prior Mean of Gamma
    lambda, # Current value of penalty parameters, 1 for smooth, 2 for flat
    lpd, # Log posterior density function of lambda
    ls, # Log of standard deviation of proposal distribution
    weight, # Window Function; NULL means default quadratic window,
    Ms,
    group,
    verbose) {
  require(mvtnorm)
  p <- dim(gamma)[1]
  K <- dim(gamma)[2]
  mat <- foreach(i = 1:K, .packages = c("matrixStats", "mvtnorm"), .combine = rbind) %dorng% {
    source("temp.R")
    Q <- penalty_Matrix(p, lambda[1], lambda[2], weight)
    n <- puMVN(mu, Q$V, 0, Ms, verbose = verbose)
    l <- dmvnorm(t(gamma)[i, ], mu, Q$V, log = TRUE)
    c(l, n)
  }
  norm.const <- mat[, -1]
  lpdf <- mat[, 1]
  group.id <- sort(unique(group))
  nC <- 0
  for (i in 1:length(group.id)) {
    temp <- colSums(norm.const[group == group.id[i], ])
    nC <- nC + logSumExp(temp)
  }
  ll <- sum(lpdf) - nC + lpd(lambda)
  # Propose new state
  # Since lambda is positive, we use log-normal jumps
  new <- exp(log(lambda) + rnorm(2, sd = exp(ls)))
  mat <- foreach(i = 1:K, .packages = c("matrixStats", "mvtnorm"), .combine = rbind) %dorng% {
    source("temp.R")
    Q <- penalty_Matrix(p, new[1], new[2], weight)
    n <- puMVN(mu, Q$V, 0, Ms, verbose = verbose)
    l <- dmvnorm(t(gamma)[i, ], mu, Q$V, log = TRUE)
    c(l, n)
  }
  norm.const <- mat[, -1]
  lpdf <- mat[, 1]
  group.id <- sort(unique(group))
  nC <- 0
  for (i in 1:length(group.id)) {
    temp <- colSums(norm.const[group == group.id[i], ])
    nC <- nC + logSumExp(temp)
  }
  ll2 <- sum(lpdf) - nC + lpd(new)

  diff <- ll - sum(log(new)) - ll2 + sum(log(lambda))
  # diff <- ll - ll2
  temp <- rexp(1)
  if (temp > diff) {
    acc_status <- 1
    state <- new
  } else {
    acc_status <- 0
    state <- lambda
  }
  return(list(acc_status = acc_status, new = state))
}
