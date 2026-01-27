# -- [E-Step] Compute P(Y | x_v, Z) for S=0 ----
compute_prob_matrix <- function(theta, y0, X0, x_support, x_name) {
  p_covs <- ncol(X0)
  beta   <- theta[1:p_covs]
  alpha  <- theta[(p_covs + 1):length(theta)]

  # Pre-calculate the linear predictor for Z (excluding X)
  x_col_idx <- which(colnames(X0) == x_name)
  beta_no_x <- beta[-x_col_idx]
  X0_no_x   <- X0[, -x_col_idx, drop = FALSE]
  # Constant part of mu for each subject: mu_fixed = Z * gamma
  mu_fixed <- as.vector(X0_no_x %*% beta_no_x)

  # Vectorized calculation of mu for all support points v
  # mu_iv = mu_fixed_i + beta_x * x_v
  beta_x <- beta[x_col_idx]
  mu_mat <- outer(mu_fixed, beta_x * x_support, "+")  # (n0, d_size)

  # Calculate Ordered Probit probabilities for all i, v simultaneously
  alpha_ext <- c(-Inf, alpha, Inf)
  probs <- pnorm(alpha_ext[y0 + 1] - mu_mat) - pnorm(alpha_ext[y0] - mu_mat)
  probs[probs < 1e-12] <- 1e-12  # Prevent exactly zero probabilities

  return(probs)
}

# -- [M-Step] Weighted Negative Log-Likelihood for Probit ----
weighted_probit_nll <- function(theta, y, X, weights) {
  p_covs <- ncol(X)
  beta   <- theta[1:p_covs]
  alpha  <- theta[(p_covs + 1):length(theta)]

  # Constraint: Cutpoints must be strictly increasing
  if (any(diff(alpha) <= 0)) return(1e10)

  # Probability Calculation: P(Y=k) = Phi(alpha_k - mu) - Phi(alpha_{k-1} - mu)
  mu <- as.vector(X %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  p_y <- pnorm(alpha_ext[y + 1] - mu) - pnorm(alpha_ext[y] - mu)
  p_y <- pmax(p_y, 1e-16)  # Numerical Stability
  log_lik <- sum(weights * log(p_y))
  if (!is.finite(log_lik)) return(1e15)  # Prevent non-finite returns to optim

  return(-log_lik)
}

# -- [M-Step] Gradient for Weighted Probit ----
weighted_probit_grad <- function(theta, y, X, weights) {
  p_covs <- ncol(X)
  beta   <- theta[1:p_covs]
  alpha  <- theta[(p_covs + 1):length(theta)]
  K_minus_1 <- length(alpha)

  if (any(diff(alpha) <= 0)) return(rep(0, length(theta)))

  mu <- as.vector(X %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  z_u <- alpha_ext[y + 1] - mu
  z_l <- alpha_ext[y] - mu
  phi_u <- dnorm(z_u)
  phi_l <- dnorm(z_l)
  p_y <- pmax(pnorm(z_u) - pnorm(z_l), 1e-16)  # Numerical Stability

  # Gradient for Beta
  d_mu <- (phi_l - phi_u) / p_y
  grad_beta <- -colSums(weights * d_mu * X)

  # Gradient for Alpha
  grad_alpha <- numeric(K_minus_1)
  for (j in 1:K_minus_1) {
    # Indicator if alpha_j is upper bound vs lower bound
    d_alpha_j <- (y == j) * (phi_u / p_y) - (y == j + 1) * (phi_l / p_y)
    grad_alpha[j] <- -sum(weights * d_alpha_j)
  }

  return(c(grad_beta, grad_alpha))
}

# -- [M-Step] Update theta in the M-step ----
update_theta_weighted <- function(start_theta, q_iv, y1, X1, y0, X0, x_support, x_name) {
  n1 <- length(y1)
  n0 <- length(y0)
  d_size <- length(x_support)

  # Prepare S=1 data: Weights are simply 1

  # Prepare S=0 data: Repeat each person in S=0 for every support point x_v
  y0_long <- rep(y0, each = d_size)
  weights0_long <- as.vector(t(q_iv))  # (n0, d_size)

  X0_long <- X0[rep(1:n0, each = d_size), , drop = FALSE]
  X0_long[, x_name] <- rep(x_support, n0)

  # Combine S=1 and S=0
  y_full <- c(y1, y0_long)
  X_full <- rbind(X1, X0_long)
  weights_full <- c(rep(1, n1), weights0_long)

  # Optimize
  fit <- optim(
    par = start_theta,
    fn = weighted_probit_nll,
    gr = weighted_probit_grad,
    y = y_full,
    X = X_full,
    weights = weights_full,
    method = "BFGS"
  )

  return(fit$par)
}

# -- Calculate log-likelihood at fixed theta, p_vl ----
calculate_full_loglik <- function(theta, p_vl, y1, X1, y0, X0, B0, x_support, x_name) {
  # --- Phase 2 (S=1): Full Data ---
  p_covs <- ncol(X1)
  beta   <- theta[1:p_covs]
  alpha  <- theta[(p_covs + 1):length(theta)]

  mu1 <- as.vector(X1 %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  prob1 <- pnorm(alpha_ext[y1 + 1] - mu1) - pnorm(alpha_ext[y1] - mu1)

  # --- Phase 1 (S=0): Marginal Data ---
  # We need P(Y|Z) = sum_v sum_l P(Y|x_v, Z) * B_l(Z) * p_vl
  # Use your existing compute_prob_matrix for P(Y|x_v, Z)
  prob_y0_v <- compute_prob_matrix(theta, y0, X0, x_support, x_name)
  # marginal_prob_i = sum_l [ B_l(Z_i) * sum_v (P(Y_i|x_v, Z_i) * p_vl) ]
  # Inner sum over v: [n0 x s_n]
  inner_sum <- prob_y0_v %*% p_vl
  # Outer sum over l: multiply by B-spline basis values
  prob0 <- rowSums(inner_sum * B0)

  return(sum(log(prob1 + 1e-16)) + sum(log(prob0 + 1e-16)))
}

# -- Profile log-likelihood ----
pl_theta <- function(theta, p_vl_init, p_vl_R1, y1, X1, y0, X0, B0, x_support,
                     x_name, max_iter = 100, tol = 1e-8) {
  # Pre-compute prob_y0_v for FIXED theta
  prob_y0_v <- compute_prob_matrix(theta, y0, X0, x_support, x_name)

  # EM Loop: Update only p_vl
  p_vl_curr <- update_p_vl_cpp(prob_y0_v, B0, p_vl_init, p_vl_R1, max_iter, tol)

  return(calculate_full_loglik(theta, p_vl_curr, y1, X1, y0, X0, B0, x_support, x_name))
}

# -- Estimate SE via profile likelihood ----
estimate_smle_se <- function(theta_hat, p_vl_hat, p_vl_R1,
                             y1, X1, y0, X0, B0, x_support, x_name) {
  # Define the function to differentiate
  obj_func <- function(t) {
    pl_theta(theta = t, p_vl_init = p_vl_hat, p_vl_R1,
             y1, X1, y0, X0, B0, x_support, x_name)
  }

  # Calculate Hessian numerically
  hess <- numDeriv::hessian(func = obj_func, x = theta_hat)

  # Estimate Covariance and SE
  vcov_mat <- try(solve(-hess), silent = TRUE)
  if (inherits(vcov_mat, "try-error")) {
    return(list(se = rep(NA, length(theta_hat)), vcov = NULL))
  }
  se <- sqrt(diag(vcov_mat))

  return(list(se = se, vcov = vcov_mat))
}
