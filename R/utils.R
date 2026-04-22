# -- Model 1 (ordinal + probit) Helpers ----
## -- ACML Helpers ----
# Spread collapsed cutpoints for initialization.
# When some Y levels are unobserved, polr collapses their cutpoints to nearly
# identical values. In the log-spacing parameterization, this makes exp(gamma) ~ 0,
# causing vanishing gradients. Fix: evenly space all cutpoints between the first
# and last (which polr identifies correctly as boundaries of observed groups).
# called by orm_acml2()
spread_collapsed_cutpoints <- function(alpha) {
  q <- length(alpha)
  if (q <= 1) return(alpha)
  # Check if any spacings are tiny relative to the total range
  total_range <- alpha[q] - alpha[1]
  min_spacing <- min(diff(alpha))
  if (min_spacing > total_range / (10 * q)) return(alpha)  # Already well-spread
  seq(alpha[1], alpha[q], length.out = q)
}

# Negative ascertainment-corrected log-likelihood
# called by orm_acml()
acml_probit_nll <- function(theta, yvec, Xmat, pi_perY, family) {
  K <- length(pi_perY)
  p <- ncol(Xmat)

  # Extract parameters
  beta <- theta[1:p]
  alpha <- theta[(p + 1):(p + K - 1)]
  if (any(diff(alpha) <= 0)) return(1e10)  # Constraint: Cutpoints must be strictly increasing
  mu <- as.vector(Xmat %*% beta)  # Linear predictor

  # Calculate category probabilities P(Y=m | X)
  p_m <- matrix(0, nrow = nrow(Xmat), ncol = K)
  if (family == "probit") {
    F_link <- stats::pnorm
  } else if (family == "logistic") {
    F_link <- stats::plogis
  }
  p_m[, 1] <- F_link(alpha[1] - mu)
  if (K > 2) {
    for (m in 2:(K-1)) {
      p_m[, m] <- F_link(alpha[m] - mu) - F_link(alpha[m-1] - mu)
    }
  }
  p_m[, K] <- 1 - F_link(alpha[K-1] - mu)

  # Ascertainment-corrected likelihood
  # Numerator: Probability of the observed outcome
  yobs_idx <- as.numeric(yvec)
  p_yobs <- p_m[cbind(1:nrow(Xmat), yobs_idx)]
  p_yobs <- pmax(p_yobs, 1e-16)  # Prevent zero or negative probabilities
  # Denominator: Weighted sum of probabilities (Probability of selection)
  denom <- pmax(1e-16, as.vector(p_m %*% pi_perY))  # Numerical stability

  return(-sum(log(p_yobs) - log(denom)))
}

# Analytical derivatives for ACML optimization
# called by orm_acml()
acml_probit_grad <- function(theta, yvec, Xmat, pi_perY, family) {
  K <- length(pi_perY)
  n <- nrow(Xmat)
  p <- ncol(Xmat)

  # Extract parameters
  beta <- theta[1:p]
  alpha <- theta[(p + 1):(p + K - 1)]
  mu <- as.vector(Xmat %*% beta)
  yobs <- as.numeric(yvec)

  # Pre-calculate normal PDF (phi) and CDF (Phi) at each cutpoint
  phi_mat <- matrix(0, nrow = n, ncol = K+1)
  Phi_mat <- matrix(0, nrow = n, ncol = K+1)
  Phi_mat[, K + 1] <- 1  # Phi(Inf) = 1
  if (family == "probit") {
    f_link <- stats::dnorm
    F_link <- stats::pnorm
  } else if (family == "logistic") {
    f_link <- stats::dlogis
    F_link <- stats::plogis
  }
  for(j in 1:(K - 1)) {
    z <- alpha[j] - mu
    Phi_mat[, j + 1] <- F_link(z)
    phi_mat[, j + 1] <- f_link(z)
  }
  p_m <- Phi_mat[, 2:(K+1)] - Phi_mat[, 1:K]
  D <- as.vector(p_m %*% pi_perY)

  # --- Gradient for beta ---
  # d(p_im)/d(beta) = X * (phi_{m-1} - phi_m)
  dp_dbeta_obs <- Xmat * (phi_mat[cbind(1:n, yobs)] - phi_mat[cbind(1:n, yobs+1)])
  dD_dbeta <- matrix(0, nrow = n, ncol = p)
  for(m in 1:K) {
    dD_dbeta <- dD_dbeta + pi_perY[m] * (Xmat * (phi_mat[, m] - phi_mat[, m+1]))
  }
  grad_beta <- -colSums(dp_dbeta_obs / p_m[cbind(1:n, yobs)] - dD_dbeta / D)

  # --- Gradient for alpha ---
  grad_alpha <- numeric(K - 1)
  for(j in 1:(K-1)) {
    # Derivative of observed part w.r.t alpha_j
    d_obs_alpha_j <- (yobs == j) * phi_mat[, j+1] - (yobs == j+1) * phi_mat[, j+1]
    # Derivative of selection part w.r.t alpha_j
    term2_alpha <- (pi_perY[j] - pi_perY[j+1]) * phi_mat[, j+1] / D
    grad_alpha[j] <- -sum(d_obs_alpha_j / p_m[cbind(1:n, yobs)] - term2_alpha)
  }

  return(c(grad_beta, grad_alpha))
}

## -- SMLE Helpers ----
# Compute P(Y_i | x_v, Z_i; theta) for S=0 across x_v in x_support
# called by orm_smle(), em_joint()
compute_py_given_xv_s0 <- function(theta, yvec_s0, Xmat_s0, x_support, x_colname, family = "probit",
                                   x_inter_colname = NULL, z_inter_colname = NULL) {
  pcovs <- ncol(Xmat_s0)
  beta <- theta[1:pcovs]
  alpha <- theta[(pcovs + 1):length(theta)]

  # Pre-calculate the linear predictor for Z (excluding all X-dependent columns)
  all_xcol_idx <- which(colnames(Xmat_s0) %in% c(x_colname, x_inter_colname))
  beta_Z <- beta[-all_xcol_idx]
  Zmat_s0 <- Xmat_s0[, -all_xcol_idx, drop = FALSE]
  # Constant part of mu for each subject: mu_fixed = Z %*% beta_z
  mu_fixed <- as.vector(Zmat_s0 %*% beta_Z)

  # Vectorized calculation of mu for all support points v
  # Main X effect: mu_iv += x_support[v,] %*% beta_x_main
  xcol_idx <- which(colnames(Xmat_s0) %in% x_colname)
  beta_x <- beta[xcol_idx]
  mu_mat <- outer(mu_fixed, as.vector(x_support %*% beta_x), "+")  # (n0, d_size)

  # Interaction contribution: mu_iv += Z_inter[i,] . (x_support[v,] %*% beta_inter reshaped)
  if (length(x_inter_colname) > 0) {
    inter_idx <- which(colnames(Xmat_s0) %in% x_inter_colname)
    beta_inter <- beta[inter_idx]
    Z_inter <- Xmat_s0[, z_inter_colname, drop = FALSE]  # (n0, n_z)
    n_x <- ncol(x_support)
    n_z <- ncol(Z_inter)
    # R model.matrix orders interaction columns with X varying fastest within Z
    beta_inter_mat <- matrix(beta_inter, nrow = n_x, ncol = n_z, byrow = FALSE)
    # C[v, j] = x_support[v,] %*% beta_inter_mat[, j]; C is (d_size, n_z)
    # inter_contrib[i, v] = Z_inter[i,] . C[v,] = Z_inter %*% t(C)
    C <- x_support %*% beta_inter_mat  # (d_size, n_z)
    mu_mat <- mu_mat + Z_inter %*% t(C)  # (n0, d_size)
  }

  # Calculate ordered probit probabilities for all i, v simultaneously
  alpha_ext <- c(-Inf, alpha, Inf)
  if (family == "probit") {
    F_link <- stats::pnorm
  } else if (family == "logistic") {
    F_link <- stats::plogis
  }
  probs <- F_link(alpha_ext[yvec_s0 + 1] - mu_mat) - F_link(alpha_ext[yvec_s0] - mu_mat)
  probs[probs < 1e-16] <- 1e-16  # Prevent exactly zero probabilities

  return(probs)
}

## -- E-Step in EM Algorithm ----
# Compute w_iv = E[I_iv | Y_i, Z_i, theta(m), p_vl] for S=0
# called by orm_smle(), em_joint()
compute_w_iv_s0 <- function(prob_y_given_xv_s0, p_vl, Bbasis_s0) {
  # P(X = x_v | Z) = sum_l { p_vl * B_l(Z) }
  prob_xv_s0 <- Bbasis_s0 %*% t(p_vl)
  w_iv_num <- prob_y_given_xv_s0 * prob_xv_s0
  w_denom <- pmax(1e-16, rowSums(w_iv_num))

  return(w_iv_num / w_denom)
}

## -- M-Step in EM Algorithm ----
# Update theta in M-step
# called by orm_smle()
update_theta_smle <- function(theta, w_iv, yvec_s1, yvec_s0, Xmat_s1, Xmat_s0, x_support, x_colname, family,
                              x_inter_colname = NULL, z_inter_colname = NULL) {
  n1 <- length(yvec_s1)
  n0 <- length(yvec_s0)
  d_size <- nrow(x_support)

  # Prepare S=1 data: Weights are simply 1

  # Prepare S=0 data: Repeat each person in S=0 for every support point x_v
  yvec_s0_long <- rep(yvec_s0, each = d_size)
  w_iv_s0_long <- as.vector(t(w_iv))  # (n0 * d_size)

  Xmat_s0_long <- Xmat_s0[rep(1:n0, each = d_size), , drop = FALSE]
  Xmat_s0_long[, x_colname] <- x_support[rep(1:nrow(x_support), n0), , drop = FALSE]

  # Fill in interaction columns: each is the product of corresponding X main and Z columns
  if (length(x_inter_colname) > 0) {
    n_x <- length(x_colname)
    n_z <- length(z_inter_colname)
    for (j in seq_len(n_z)) {
      for (k in seq_len(n_x)) {
        inter_col <- x_inter_colname[(j - 1) * n_x + k]
        Xmat_s0_long[, inter_col] <- Xmat_s0_long[, x_colname[k]] * Xmat_s0_long[, z_inter_colname[j]]
      }
    }
  }

  # Combine S=1 and S=0
  yvec_full <- c(yvec_s1, yvec_s0_long)
  Xmat_full <- rbind(Xmat_s1, Xmat_s0_long)
  w_iv_full <- c(rep(1, n1), w_iv_s0_long)

  # Optimize
  fit <- stats::optim(
    par = theta, fn = weighted_nll, gr = weighted_grad,
    yvec = yvec_full, Xmat = Xmat_full, w_iv = w_iv_full, family = family, method = "BFGS"
  )

  return(fit$par)
}

# Weighted negative log-likelihood for SMLE optimization
# called by update_theta_smle()
weighted_nll <- function(theta, yvec, Xmat, w_iv, family) {
  pcovs <- ncol(Xmat)
  beta <- theta[1:pcovs]
  alpha <- theta[(pcovs + 1):length(theta)]

  # Constraint: Cutpoints must be strictly increasing
  if (any(diff(alpha) <= 0)) return(1e10)

  # Probability Calculation: P(Y=k) = Phi(alpha_k - mu) - Phi(alpha_{k-1} - mu)
  mu <- as.vector(Xmat %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  if (family == "probit") {
    F_link <- stats::pnorm
  } else if (family == "logistic") {
    F_link <- stats::plogis
  }
  p_y <- F_link(alpha_ext[yvec + 1] - mu) - F_link(alpha_ext[yvec] - mu)
  p_y <- pmax(1e-16, p_y)  # Numerical stability

  return(-sum(w_iv * log(p_y)))
}

# Analytical derivatives for SMLE optimization
# called by update_theta_smle()
weighted_grad <- function(theta, yvec, Xmat, w_iv, family) {
  pcovs <- ncol(Xmat)
  beta <- theta[1:pcovs]
  alpha <- theta[(pcovs + 1):length(theta)]

  if (any(diff(alpha) <= 0)) return(rep(0, length(theta)))

  mu <- as.vector(Xmat %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  z_u <- alpha_ext[yvec+1] - mu; z_l <- alpha_ext[yvec] - mu
  if (family == "probit") {
    f_link <- stats::dnorm
    F_link <- stats::pnorm
  } else if (family == "logistic") {
    f_link <- stats::dlogis
    F_link <- stats::plogis
  }
  phi_u <- f_link(z_u); phi_l <- f_link(z_l)
  p_y <- pmax(1e-16, F_link(z_u) - F_link(z_l))  # Numerical Stability

  # Gradient for beta
  d_mu <- (phi_l - phi_u) / p_y
  grad_beta <- -colSums(w_iv * d_mu * Xmat)

  # Gradient for alpha
  K_minus_1 <- length(alpha)
  grad_alpha <- numeric(K_minus_1)
  for (j in 1:K_minus_1) {
    # Indicator if alpha_j is upper bound vs lower bound
    d_alpha_j <- (yvec == j) * (phi_u / p_y) - (yvec == j+1) * (phi_l / p_y)
    grad_alpha[j] <- -sum(w_iv * d_alpha_j)
  }

  return(c(grad_beta, grad_alpha))
}

## -- SE Estimation ----
# Estimate SE via profile likelihood
# called by orm_smle()
# method = "forward":  manual second-order forward-difference Hessian with h_n = n^(-1/2) * h_n_scale
# method = "numDeriv": Richardson-extrapolated Hessian via numDeriv::hessian()
estimate_se_smle <- function(theta, p_vl, p_vl_s1, yvec_s1, yvec_s0, Xmat_s1, Xmat_s0, Bbasis_s0, Bbasis_s1,
                             s1_support_idx, x_support, x_colname, family,
                             method = "forward", h_n_scale = 0.1, hessian_method_args = list(), max_iter, tol,
                             verbose = FALSE, x_inter_colname = NULL, z_inter_colname = NULL) {
  se_pll_max_iters <- integer(0)
  theta_perturb <- numeric(0)
  p_vl_moves <- numeric(0)
  pll_func <- function(t) {
    result <- smle_pll(
      theta = t, p_vl = p_vl, p_vl_s1 = p_vl_s1, yvec_s1 = yvec_s1, yvec_s0 = yvec_s0, Xmat_s1 = Xmat_s1, Xmat_s0 = Xmat_s0,
      Bbasis_s0 = Bbasis_s0, Bbasis_s1 = Bbasis_s1, s1_support_idx = s1_support_idx,
      x_support = x_support, x_colname = x_colname, family = family, max_iter = max_iter, tol = tol,
      x_inter_colname = x_inter_colname, z_inter_colname = z_inter_colname)
    se_pll_max_iters[[length(se_pll_max_iters) + 1]] <<- attr(result, "pll_inner_iter")
    theta_perturb[[length(theta_perturb) + 1]] <<- max(abs(t - theta))
    p_vl_moves[[length(p_vl_moves) + 1]] <<- attr(result, "p_vl_move")
    result
  }

  if (method == "numDeriv") {
    hess <- numDeriv::hessian(pll_func, theta, method.args = hessian_method_args)
    hn <- NA_real_
  } else {
    # Manual second-order forward-difference Hessian:
    #   H[i,j] = (f(x+h*ei+h*ej) - f(x+h*ei) - f(x+h*ej) + f(x)) / h^2, h = h_n_scale * n^(-1/2)
    nparams <- length(theta)
    n <- length(yvec_s1) + length(yvec_s0)
    hn <- h_n_scale * n^(-1/2)
    e_mat <- diag(hn, nparams)

    pl_0d <- pll_func(theta)
    pl_1d <- numeric(nparams)
    for (i in seq_len(nparams)) {
      pl_1d[i] <- pll_func(theta + e_mat[i, ])
    }
    pl_2d <- matrix(NA, nparams, nparams)
    for (i in seq_len(nparams)) {
      for (j in i:nparams) {
        pl_2d[i, j] <- pl_2d[j, i] <- pll_func(theta + e_mat[i, ] + e_mat[j, ])
      }
    }

    hess <- matrix(NA, nparams, nparams)
    for (i in seq_len(nparams)) {
      for (j in i:nparams) {
        hess[i, j] <- hess[j, i] <- (pl_2d[i, j] - pl_1d[i] - pl_1d[j] + pl_0d) / (hn^2)
      }
    }
  }

  if (verbose && length(se_pll_max_iters) > 0) {
    n_evals <- length(se_pll_max_iters)
    n_hit_max <- sum(se_pll_max_iters >= max_iter)
    cat(sprintf("SE (%s): %d PLL evals | inner iters %d-%d (max=%d, hit-max=%d)\n",
                method, n_evals, min(se_pll_max_iters), max(se_pll_max_iters), max_iter, n_hit_max))
    if (method == "numDeriv") {
      cat(sprintf("  theta perturb |h|:   min=%.2e  med=%.2e  max=%.2e\n",
                  min(theta_perturb), stats::median(theta_perturb), max(theta_perturb)))
    } else {
      cat(sprintf("  theta perturb h_n=%.2e\n", hn))
    }
    cat(sprintf("  p_vl move per eval:  min=%.2e  med=%.2e  max=%.2e  (se_tol=%.2e)\n",
                min(p_vl_moves), stats::median(p_vl_moves), max(p_vl_moves), tol))
  }

  diagnostics <- list(method = method, hn = hn, hessian_method_args = hessian_method_args,
                      inner_iters = se_pll_max_iters, theta_perturb = theta_perturb, p_vl_moves = p_vl_moves,
                      hess = hess)

  # Estimate Covariance and SE
  vcov_mat <- try(solve(-hess), silent = TRUE)
  if (inherits(vcov_mat, "try-error")) {
    warning("Hessian inversion failed; standard errors cannot be computed.")
    return(list(se = rep(NA, length(theta)), vcov = NULL, se_max_iterations = NULL,
                diagnostics = diagnostics))
  }

  return(list(se = sqrt(diag(vcov_mat)), vcov = vcov_mat,
              se_max_iterations = max(se_pll_max_iters), diagnostics = diagnostics))
}

# Profile log-likelihood
# called by estimate_se_smle()
smle_pll <- function(theta, p_vl, p_vl_s1, yvec_s1, yvec_s0, Xmat_s1, Xmat_s0, Bbasis_s0, Bbasis_s1, s1_support_idx,
                         x_support, x_colname, family, max_iter, tol, x_inter_colname = NULL, z_inter_colname = NULL) {
  # Pre-compute prob_y0_v for FIXED theta
  prob_y_given_xv_s0 <- compute_py_given_xv_s0(
    theta = theta, yvec_s0 = yvec_s0, Xmat_s0 = Xmat_s0, x_support = x_support, x_colname = x_colname, family = family,
    x_inter_colname = x_inter_colname, z_inter_colname = z_inter_colname)

  # EM Loop: Update only p_vl
  p_vl_opt <- update_p_vl_cpp(prob_y_given_xv_s0, Bbasis_s0, p_vl, p_vl_s1, max_iter, tol)
  p_vl_move <- max(abs(as.vector(p_vl_opt) - as.vector(p_vl)))

  # Pass pre-computed prob_y_given_xv_s0 to avoid redundant computation
  ll <- smle_ll(
    theta = theta, p_vl = p_vl_opt, yvec_s1 = yvec_s1, yvec_s0 = yvec_s0, Xmat_s1 = Xmat_s1, Xmat_s0 = Xmat_s0,
    Bbasis_s0 = Bbasis_s0, Bbasis_s1 = Bbasis_s1, s1_support_idx = s1_support_idx,
    x_support = x_support, x_colname = x_colname, family = family, prob_y_given_xv_s0 = prob_y_given_xv_s0,
    x_inter_colname = x_inter_colname, z_inter_colname = z_inter_colname)
  attr(ll, "pll_inner_iter") <- as.integer(attr(p_vl_opt, "iterations"))
  attr(ll, "p_vl_move") <- p_vl_move
  return(ll)
}

# Calculate log-likelihood at fixed theta, p_vl
# called by estimate_se_smle()
# The SMLE objective has three additive pieces:
#   S=1:  log P(Y|X,Z; theta) (in prob_s1)
#   S=1:  log p(X|Z; p_vl) (in prob_x_s1)
#   S=0:  log sum_v P(Y|x_v,Z) * sum_l B_l(Z) p_vl (in prob_s0)
smle_ll <- function(theta, p_vl, yvec_s1, yvec_s0, Xmat_s1, Xmat_s0, Bbasis_s0, Bbasis_s1, s1_support_idx,
                        x_support, x_colname, family, prob_y_given_xv_s0 = NULL,
                        x_inter_colname = NULL, z_inter_colname = NULL) {
  # For S=1, we need P(Y|X,Z)
  pcovs <- ncol(Xmat_s1)
  beta <- theta[1:pcovs]
  alpha <- theta[(pcovs + 1):length(theta)]
  mu_s1 <- as.vector(Xmat_s1 %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  if (family == "probit") {
    F_link <- stats::pnorm
  } else if (family == "logistic") {
    F_link <- stats::plogis
  }
  prob_s1 <- pmax(1e-16, F_link(alpha_ext[yvec_s1 + 1] - mu_s1) - F_link(alpha_ext[yvec_s1] - mu_s1))

  # For S=1, sieve density p(X_i | Z_i) = sum_l B_l(Z_i) * p_{v_i, l}
  prob_xv_s1 <- Bbasis_s1 %*% t(p_vl)  # (n1, d)
  prob_x_s1 <- pmax(1e-16, prob_xv_s1[cbind(seq_along(s1_support_idx), s1_support_idx)])

  # For S=0, we need P(Y|Z) = sum_v sum_l P(Y|x_v, Z) * B_l(Z) * p_vl
  if (is.null(prob_y_given_xv_s0)) {
    prob_y_given_xv_s0 <- compute_py_given_xv_s0(
      theta = theta, yvec_s0 = yvec_s0, Xmat_s0 = Xmat_s0, x_support = x_support, x_colname = x_colname, family = family,
      x_inter_colname = x_inter_colname, z_inter_colname = z_inter_colname)
  }
  prob_xv_s0 <- Bbasis_s0 %*% t(p_vl)
  prob_s0 <- pmax(1e-16, rowSums(prob_y_given_xv_s0 * prob_xv_s0))

  return(sum(log(prob_s1)) + sum(log(prob_x_s1)) + sum(log(prob_s0)))
}

## -- MLE0 Helpers ----
# Compute P(Y_i = y_i | support_j; theta) for all (i, j) combinations
# called by orm_mle0()
compute_py_given_support <- function(theta, yvec, support_mat, family) {
  pcovs <- ncol(support_mat)
  beta <- theta[1:pcovs]
  alpha <- theta[(pcovs + 1):length(theta)]

  mu_support <- as.vector(support_mat %*% beta)  # (m,)
  alpha_ext <- c(-Inf, alpha, Inf)

  if (family == "probit") {
    F_link <- stats::pnorm
  } else if (family == "logistic") {
    F_link <- stats::plogis
  }

  # P(Y_i = y_i | support_j) = F(alpha_{y_i} - mu_j) - F(alpha_{y_i-1} - mu_j)
  upper <- outer(alpha_ext[yvec + 1], mu_support, "-")  # (n, m)
  lower <- outer(alpha_ext[yvec], mu_support, "-")  # (n, m)
  probs <- F_link(upper) - F_link(lower)
  probs[probs < 1e-16] <- 1e-16

  return(probs)
}

# Compute psi_ij = E[I{covariates_i = support_j} | Y_i; theta, q_j] for S=0
# called by orm_mle0()
compute_psi_mle0 <- function(prob_y_given_support, q_j) {
  psi_num <- sweep(prob_y_given_support, 2, q_j, "*")  # (n0, m)
  psi_denom <- pmax(1e-16, rowSums(psi_num))
  return(psi_num / psi_denom)
}

# M-step: update theta by maximizing the weighted log-likelihood
# called by orm_mle0()
update_theta_mle0 <- function(theta, psi_ij, yvec_s1, yvec_s0, Xmat_s1, support_mat, family) {
  n1 <- length(yvec_s1)
  n0 <- length(yvec_s0)
  m <- nrow(support_mat)

  # S=0: expand each subject across all m support points
  yvec_s0_long <- rep(yvec_s0, each = m)
  psi_s0_long <- as.vector(t(psi_ij))  # (n0 * m,)
  Xmat_s0_long <- support_mat[rep(1:m, n0), , drop = FALSE]

  # Combine S=1 and S=0
  yvec_full <- c(yvec_s1, yvec_s0_long)
  Xmat_full <- rbind(Xmat_s1, Xmat_s0_long)
  w_full <- c(rep(1, n1), psi_s0_long)

  fit <- stats::optim(
    par = theta, fn = weighted_nll, gr = weighted_grad,
    yvec = yvec_full, Xmat = Xmat_full, w_iv = w_full, family, method = "BFGS"
  )

  return(fit$par)
}

# M-step: update q_j = (1/n) * sum_i psi_ij
# called by orm_mle0()
update_q_mle0 <- function(psi_ij, s1_support_idx, n) {
  m <- ncol(psi_ij)
  s1_counts <- tabulate(s1_support_idx, nbins = m)
  s0_sums <- colSums(psi_ij)
  q_j <- (s1_counts + s0_sums) / n
  return(q_j)
}

# Log-likelihood for MLE0
# called by mle0_pll()
mle0_ll <- function(theta, q_j, yvec_s1, yvec_s0, Xmat_s1, support_mat, s1_support_idx, family,
                    prob_y_given_support_s0 = NULL) {
  pcovs <- ncol(Xmat_s1)
  beta <- theta[1:pcovs]
  alpha <- theta[(pcovs + 1):length(theta)]

  # S=1: log P(Y_i | X_i; theta) + log q_{j(i)}
  mu_s1 <- as.vector(Xmat_s1 %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  if (family == "probit") {
    F_link <- stats::pnorm
  } else if (family == "logistic") {
    F_link <- stats::plogis
  }
  prob_s1 <- pmax(1e-16, F_link(alpha_ext[yvec_s1 + 1] - mu_s1) - F_link(alpha_ext[yvec_s1] - mu_s1))
  ll_s1 <- sum(log(prob_s1)) + sum(log(pmax(1e-16, q_j[s1_support_idx])))

  # S=0: log sum_j P(Y_i | support_j; theta) * q_j
  if (is.null(prob_y_given_support_s0)) {
    prob_y_given_support_s0 <- compute_py_given_support(theta, yvec_s0, support_mat, family)
  }
  prob_s0 <- pmax(1e-16, as.vector(prob_y_given_support_s0 %*% q_j))
  ll_s0 <- sum(log(prob_s0))

  return(ll_s1 + ll_s0)
}

# Profile log-likelihood for MLE0: hold theta fixed, optimize q_j
# called by estimate_se_mle0()
mle0_pll <- function(theta, q_j, yvec_s1, yvec_s0, Xmat_s1, support_mat, s1_support_idx, family, max_iter, tol) {
  n <- length(yvec_s1) + length(yvec_s0)

  # Pre-compute P(Y_i | support_j; theta) for S=0 (fixed theta)
  prob_y_s0 <- compute_py_given_support(theta, yvec_s0, support_mat, family)

  # EM loop: update only q_j
  for (iter in 1:max_iter) {
    q_old <- q_j
    psi_ij <- compute_psi_mle0(prob_y_s0, q_j)
    q_j <- update_q_mle0(psi_ij, s1_support_idx, n)
    if (max(abs(q_j - q_old)) < tol) break
  }

  ll <- mle0_ll(theta, q_j, yvec_s1, yvec_s0, Xmat_s1, support_mat, s1_support_idx, family, prob_y_s0)
  attr(ll, "pll_inner_iter") <- as.integer(iter)
  return(ll)
}

# Estimate SE for MLE0 via profile likelihood Hessian
# called by orm_mle0()
estimate_se_mle0 <- function(theta, q_j, yvec_s1, yvec_s0, Xmat_s1, support_mat, s1_support_idx, family, max_iter, tol,
                             method = "forward", verbose = FALSE) {
  se_pll_max_iters <- integer(0)
  pll_func <- function(t) {
    result <- mle0_pll(t, q_j, yvec_s1, yvec_s0, Xmat_s1, support_mat, s1_support_idx, family, max_iter, tol)
    se_pll_max_iters[[length(se_pll_max_iters) + 1]] <<- attr(result, "pll_inner_iter")
    result
  }

  if (method == "numDeriv") {
    hess <- numDeriv::hessian(pll_func, theta)
  } else {
    # Second-order forward-difference Hessian
    nparams <- length(theta)
    n <- length(yvec_s1) + length(yvec_s0)
    hn <- n^(-1/2)
    e_mat <- diag(hn, nparams)

    pl_0d <- pll_func(theta)
    pl_1d <- numeric(nparams)
    for (i in seq_len(nparams)) {
      pl_1d[i] <- pll_func(theta + e_mat[i, ])
    }
    pl_2d <- matrix(NA, nparams, nparams)
    for (i in seq_len(nparams)) {
      for (j in i:nparams) {
        pl_2d[i, j] <- pl_2d[j, i] <- pll_func(theta + e_mat[i, ] + e_mat[j, ])
      }
    }

    hess <- matrix(NA, nparams, nparams)
    for (i in seq_len(nparams)) {
      for (j in i:nparams) {
        hess[i, j] <- hess[j, i] <- (pl_2d[i, j] - pl_1d[i] - pl_1d[j] + pl_0d) / (hn^2)
      }
    }
  }

  if (verbose && length(se_pll_max_iters) > 0) {
    n_evals <- length(se_pll_max_iters)
    n_hit_max <- sum(se_pll_max_iters >= max_iter)
    cat(sprintf("SE: %d PLL evaluations, inner q_j EM used %d-%d iters (max_iter = %d, %d hit max)\n",
                n_evals, min(se_pll_max_iters), max(se_pll_max_iters), max_iter, n_hit_max))
  }

  vcov_mat <- try(solve(-hess), silent = TRUE)
  if (inherits(vcov_mat, "try-error")) {
    warning("Hessian inversion failed; standard errors cannot be computed.")
    return(list(se = rep(NA, length(theta)), vcov = NULL))
  }

  return(list(se = sqrt(diag(vcov_mat)), vcov = vcov_mat, se_max_iterations = max(se_pll_max_iters)))
}

# -- Model 2 (continuous) Helpers ----
## -- E-Step in EM Algorithm ----
# Compute moments of [Y1* | Y2, X, Z, Y1, theta1(m)] for S=1
# called by em_conditional(), em_joint()
compute_y1star_s1 <- function(theta1, theta2, y1vec_s1, y2vec, Xmat_m1_s1, Xmat_m2, mu1 = NULL) {
  pcovs_m1 <- ncol(Xmat_m1_s1)
  beta <- theta1[1:pcovs_m1]
  alpha <- theta1[(pcovs_m1 + 1):length(theta1)]

  pcovs_m2 <- ncol(Xmat_m2)
  gamma <- theta2[1:pcovs_m2]
  sigma12 <- theta2[pcovs_m2 + 1]
  sigma22 <- theta2[pcovs_m2 + 2]

  # [Y1* | Y2, X, Z] ~ N(mu_cond, sd_cond)
  if (is.null(mu1)) mu1 <- as.vector(Xmat_m1_s1 %*% beta)
  mu2 <- as.vector(Xmat_m2 %*% gamma)
  mu_cond <- mu1 + (sigma12 / sigma22) * (y2vec - mu2)
  sd_cond <- sqrt(1 - (sigma12^2 / sigma22))

  alpha_ext <- c(-Inf, alpha, Inf)
  idx1 <- which(y1vec_s1 == 1); idxmax <- which(y1vec_s1 == max(y1vec_s1))
  # [Y1* | Y2, X, Z, Y1] ~ truncated normal
  a <- alpha_ext[y1vec_s1]; b <- alpha_ext[y1vec_s1 + 1]
  a_std <- (a-mu_cond) / sd_cond; b_std <- (b-mu_cond) / sd_cond
  varphi_a <- stats::dnorm(a_std); varphi_b <- stats::dnorm(b_std)
  Phi_a <- stats::pnorm(a_std); Phi_b <- stats::pnorm(b_std)
  a_std[idx1] <- 0; b_std[idxmax] <- 0
  tmp_denom <- pmax(1e-16, Phi_b - Phi_a)  # Numerical stability
  tmp1 <- (varphi_b - varphi_a) / tmp_denom
  tmp2 <- (b_std * varphi_b - a_std * varphi_a) / tmp_denom
  e_y1star <- mu_cond - sd_cond * tmp1
  var_y1star <- sd_cond^2 * (1 - tmp2 - tmp1^2)

  return(list(mean = e_y1star, var = var_y1star))
}

# Compute E[Y1* | x_v, Z, Y1, theta1(m)] for S=0
# called by em_joint()
compute_y1star_s0 <- function(theta1, y1vec_s0, Xmat_m1_s0, x_support, x_colname) {
  pcovs_m1 <- ncol(Xmat_m1_s0)
  beta <- theta1[1:pcovs_m1]
  alpha <- theta1[(pcovs_m1+1):length(theta1)]

  # Pre-calculate the linear predictor for Z (excluding X)
  xcol_idx <- which(colnames(Xmat_m1_s0) %in% x_colname)
  beta_Z <- beta[-xcol_idx]
  Zmat_m1_s0 <- Xmat_m1_s0[, -xcol_idx, drop = FALSE]
  # Constant part of mu for each subject: mu_fixed = Z %*% beta_z
  mu_fixed <- as.vector(Zmat_m1_s0 %*% beta_Z)

  # Vectorized calculation of mu for all support points v
  # mu_iv = mu_fixed_i + beta_x * x_v
  beta_x <- beta[xcol_idx]
  mu_mat <- outer(mu_fixed, as.vector(x_support %*% beta_x), "+")  # (n0, d_size)

  # Calculate truncated normal mean
  alpha_ext <- c(-Inf, alpha, Inf)
  a_mat <- alpha_ext[y1vec_s0] - mu_mat
  b_mat <- alpha_ext[y1vec_s0+1] - mu_mat

  e_y1star_s0 <- mu_mat -
    (stats::dnorm(b_mat)-stats::dnorm(a_mat)) / pmax(1e-16, stats::pnorm(b_mat)-stats::pnorm(a_mat))

  return(e_y1star_s0)
}

## -- M-Step in EM Algorithm ----
# Update theta in M-step
# called by em_conditional(), em_joint()
update_theta2 <- function(theta1, theta2, y1vec, y2vec, Xmat_m1, Xmat_m2, e_y1star, var_y1star) {
  pcovs_m1 <- ncol(Xmat_m1)
  beta <- theta1[1:pcovs_m1]
  alpha <- theta1[(pcovs_m1+1):length(theta1)]
  alpha_ext <- c(-Inf, alpha, Inf)

  pcovs_m2 <- ncol(Xmat_m2)

  # Update theta2 via regression: Y2 ~ mu2 + sigma12(Y* - mu_cond)
  mu1 <- as.vector(Xmat_m1 %*% beta)
  latent_resid <- e_y1star - mu1  # Residual of latent Y1* from its mean
  Xaug_m2 <- cbind(Xmat_m2, latent_resid)
  XtX <- t(Xaug_m2) %*% Xaug_m2
  XtX[pcovs_m2+1, pcovs_m2+1] <- sum(var_y1star + latent_resid^2)
  Xty <- t(Xaug_m2) %*% y2vec
  theta2_new <- as.vector(solve(XtX) %*% Xty)

  gamma_curr <- theta2_new[1:pcovs_m2]
  sigma12_curr <- theta2_new[pcovs_m2 + 1]
  sigma22_curr <- sum((y2vec - as.vector(Xaug_m2 %*% theta2_new))^2 +
                        sigma12_curr^2 * var_y1star) / length(y1vec) + sigma12_curr^2

  return(c(gamma_curr, sigma12_curr, sigma22_curr))
}

# Update theta in M-step
# called by em_joint()
update_theta1 <- function(theta1, y1vec_s1, y1vec_s0, Xmat_m1_s1, Xmat_m1_s0,
                          w_iv, e_y1star_s1, e_y1star_s0, x_support, x_colname) {
  n0 <- length(y1vec_s0)
  d_size <- nrow(x_support)

  # Prepare S=1 data: Weights are simply 1

  # Prepare S=0 data: each person in S=0 for every support point x_v
  y1vec_s0_long <- rep(y1vec_s0, each = d_size)
  w_iv_s0_long <- as.vector(t(w_iv))

  Xmat_m1_s0_long <- Xmat_m1_s0[rep(1:n0, each = d_size), , drop = FALSE]
  Xmat_m1_s0_long[, x_colname] <- x_support[rep(1:nrow(x_support), n0), , drop = FALSE]

  # Combine S=1 and S=0
  y1vec_full <- c(y1vec_s1, y1vec_s0_long)
  e_y1star_full <- c(e_y1star_s1, as.vector(t(e_y1star_s0)))
  Xmat_full <- rbind(Xmat_m1_s1, Xmat_m1_s0_long)
  w_iv_full <- c(rep(1, length(y1vec_s1)), w_iv_s0_long)

  # Update beta via WLS
  beta <- as.vector(solve(t(Xmat_full) %*% (w_iv_full * Xmat_full)) %*% (t(Xmat_full) %*% (w_iv_full * e_y1star_full)))
  alpha <- theta1[(length(beta) + 1):length(theta1)]

  # Optimize
  fit <- stats::optim(
    par = alpha, fn = alpha_weighted_nll, gr = alpha_weighted_grad,
    beta = beta, y1vec = y1vec_full, Xmat = Xmat_full, w_iv = w_iv_full, method = "BFGS"
  )

  return(c(beta, fit$par))
}

# Weighted negative log-likelihood for alpha optimization
# called by update_theta1()
alpha_weighted_nll <- function(alpha, beta, y1vec, Xmat, w_iv) {
  # Constraint: Cutpoints must be strictly increasing
  if (any(diff(alpha) <= 0)) return(1e10)

  # Probability Calculation: P(Y=k) = Phi(alpha_k - mu) - Phi(alpha_{k-1} - mu)
  mu <- as.vector(Xmat %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  p_y <- stats::pnorm(alpha_ext[y1vec+1] - mu) - stats::pnorm(alpha_ext[y1vec] - mu)
  p_y <- pmax(1e-16, p_y)  # Numerical Stability

  return(-sum(w_iv * log(p_y)))
}

# Analytical derivatives for alpha optimization
# called by update_theta1()
alpha_weighted_grad <- function(alpha, beta, y1vec, Xmat, w_iv) {
  if (any(diff(alpha) <= 0)) return(rep(0, length(alpha)))

  mu <- as.vector(Xmat %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  z_u <- alpha_ext[y1vec+1] - mu; z_l <- alpha_ext[y1vec] - mu
  phi_u <- stats::dnorm(z_u); phi_l <- stats::dnorm(z_l)
  p_y <- pmax(1e-16, stats::pnorm(z_u) - stats::pnorm(z_l))  # Numerical Stability

  # Gradient for alpha
  K_minus_1 <- length(alpha)
  grad_alpha <- numeric(K_minus_1)
  for (j in 1:K_minus_1) {
    # Indicator if alpha_j is upper bound vs lower bound
    d_alpha_j <- (y1vec == j) * (phi_u / p_y) - (y1vec == j+1) * (phi_l / p_y)
    grad_alpha[j] <- -sum(w_iv * d_alpha_j)
  }

  return(grad_alpha)
}

## -- SE Estimation ----
# Estimate SE
# called by em_conditional()
estimate_se_em_cond <- function(theta2, y2vec, Xmat_m2, theta_cond, vcov_cond, y1vec, Xmat_m1, mu1, pcovs_m1, alpha_ext) {
  pcovs_m2 <- ncol(Xmat_m2)
  hess <- numDeriv::hessian(
    func = function(t) {
      m2_ll(theta = t, y2vec, y1vec, Xmat_m2, mu1, alpha_ext, pcovs_m2)
    }, x = theta2)
  # Estimate the conditional variability
  V22_inv <- try(solve(-hess), silent = TRUE)

  if (inherits(V22_inv, "try-error")) {
    warning("Hessian inversion failed; standard errors cannot be computed.")
    return(list(se = rep(NA, length(theta2)), vcov = NULL))
  }

  # Calculate the Jacobian of the score function w.r.t theta1
  V21 <- numDeriv::jacobian(
    func = function(t1) {
      beta_tmp <- t1[1:pcovs_m1]
      alpha_tmp <- t1[(pcovs_m1 + 1):length(t1)]
      alpha_ext_tmp <- c(-Inf, alpha_tmp, Inf)
      mu1_tmp <- as.vector(Xmat_m1 %*% beta_tmp)
      m2_score(theta2, y2vec, y1vec, Xmat_m2, mu1_tmp, alpha_ext_tmp, pcovs_m2)
    }, x = theta_cond)

  # Law of total variance
  vcov <- V22_inv + (V22_inv %*% V21 %*% vcov_cond %*% t(V21) %*% V22_inv)

  return(list(se = sqrt(diag(vcov)), vcov = vcov))
}

m2_ll <- function(theta2, y2vec, y1vec, Xmat_m2, mu1, alpha_ext, pcovs_m2) {
  # Parameter mapping
  gamma <- theta2[1:pcovs_m2]
  sigma12 <- theta2[pcovs_m2 + 1]
  sigma22 <- theta2[pcovs_m2 + 2]

  mu2 <- as.vector(Xmat_m2 %*% gamma)
  # [Y1* | Y2, X, Z] ~ N(mu_cond, sd_cond)
  mu_cond <- mu1 + (sigma12 / sigma22) * (y2vec - mu2)
  sd_cond <- sqrt(1 - (sigma12^2 / sigma22))
  # Calculate P(Y1 | Y2, X, Z) from [Y1* | Y2, X, Z] ~ N(mu_cond, sd_cond)
  a <- alpha_ext[y1vec]; b <- alpha_ext[y1vec + 1]
  prob_y1_given_y2 <- pmax(1e-16, stats::pnorm(b, mean = mu_cond, sd = sd_cond) -
                             stats::pnorm(a, mean = mu_cond, sd = sd_cond))

  # Y2 ~ N(mu2, sigma22)
  prob_y2 <- pmax(1e-16, stats::dnorm(y2vec, mean = mu2, sd = sqrt(sigma22)))

  # log[P(Y1, Y2)] = log[P(Y1 | Y2) * P(Y2)]
  return(sum(log(prob_y1_given_y2) + log(prob_y2)))
}

m2_score <- function(theta2, y2vec, y1vec, Xmat_m2, mu1, alpha_ext, pcovs_m2) {
  # Parameter mapping
  gamma <- theta2[1:pcovs_m2]
  sigma12 <- theta2[pcovs_m2 + 1]
  sigma22 <- theta2[pcovs_m2 + 2]

  mu2 <- as.vector(Xmat_m2 %*% gamma)
  # [Y1* | Y2, X, Z] ~ N(mu_cond, sd_cond)
  mu_cond <- mu1 + (sigma12 / sigma22) * (y2vec - mu2)
  sd_cond <- sqrt(1 - (sigma12^2 / sigma22))

  # Standardized Bounds
  a_std <- (alpha_ext[y1vec] - mu_cond) / sd_cond
  b_std <- (alpha_ext[y1vec + 1] - mu_cond) / sd_cond
  pdf_a <- stats::dnorm(a_std); pdf_b <- stats::dnorm(b_std)
  cdf_a <- stats::pnorm(a_std); cdf_b <- stats::pnorm(b_std)
  idx1 <- which(y1vec == 1); a_std[idx1] <- 0
  idxmax <- which(y1vec == max(y1vec)); b_std[idxmax] <- 0
  cdf_diff <- pmax(1e-16, cdf_b - cdf_a)

  # Partial Derivatives
  # d_loglik / d_mu2
  d_mu2 <- ((pdf_b - pdf_a) / cdf_diff) * (sigma12 / (sigma22 * sd_cond)) +
    (y2vec - mu2) / sigma22
  # d_loglik / d_gamma
  score_gamma <- as.vector(t(Xmat_m2) %*% d_mu2)

  # d_loglik / d_sigma12
  da_sigma12 <- -((y2vec - mu2) / sigma22 - (sigma12 * a_std) / (sigma22 * sd_cond)) / sd_cond
  db_sigma12 <- -((y2vec - mu2) / sigma22 - (sigma12 * b_std) / (sigma22 * sd_cond)) / sd_cond
  d_sigma12 <- (pdf_b * db_sigma12 - pdf_a * da_sigma12) / cdf_diff
  score_sigma12 <- sum(d_sigma12)

  # d_loglik / d_sigma22
  db_sigma22 <- ((sigma12 * (y2vec - mu2)) / sigma22^2 - (sigma12^2 * b_std) / (2 * sigma22^2 * sd_cond)) / sd_cond
  da_sigma22 <- ((sigma12 * (y2vec - mu2)) / sigma22^2 - (sigma12^2 * a_std) / (2 * sigma22^2 * sd_cond)) / sd_cond
  d_sigma22 <- -1/(2 * sigma22) + ((y2vec - mu2)^2) / (2 * sigma22^2) +
    (pdf_b * db_sigma22 - pdf_a * da_sigma22) / cdf_diff
  score_sigma22 <- sum(d_sigma22)

  return(c(score_gamma, score_sigma12, score_sigma22))
}
