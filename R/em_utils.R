obs_loglik <- function(theta2, p_covs_m2, y2_vec, y1_vec, X_mat_m2, mu1, alpha_ext) {
  # Parameter mapping
  gamma   <- theta2[1:p_covs_m2]
  sigma12 <- theta2[p_covs_m2 + 1]
  sigma22 <- theta2[p_covs_m2 + 2]
  
  # Conditional moments of Y1* | Y2
  mu2 <- as.vector(X_mat_m2 %*% gamma)
  mu_cond <- mu1 + (sigma12 / sigma22) * (y2_vec - mu2)
  sd_cond <- sqrt(1 - (sigma12^2 / sigma22))
  # Calculate P(Y1 | Y2) from the normal distribution of Y1* | Y2
  a <- alpha_ext[y1_vec]
  b <- alpha_ext[y1_vec + 1]
  prob_y1_given_y2 <- pmax(1e-16, stats::pnorm(b, mean = mu_cond, sd = sd_cond) - 
                             stats::pnorm(a, mean = mu_cond, sd = sd_cond))
  
  # Marginal density of Y2: Y2 ~ N(mu2, sigma22)
  dens_y2 <- stats::dnorm(y2_vec, mean = mu2, sd = sqrt(sigma22))
  
  # Observed-data log-likelihood: log[P(Y1, Y2)] = log[P(Y1 | Y2) * P(Y2)]
  return(sum(log(prob_y1_given_y2) + log(dens_y2)))
}

score_theta2_function <- function(theta2, p_covs_m2, y2_vec, y1_vec, X_mat_m2, mu1, alpha_ext) {
  # Parameter mapping
  gamma   <- theta2[1:p_covs_m2]
  sigma12 <- theta2[p_covs_m2 + 1]
  sigma22 <- theta2[p_covs_m2 + 2]
  
  # Conditional moments of Y1* | Y2
  mu2 <- as.vector(X_mat_m2 %*% gamma)
  mu_cond  <- mu1 + (sigma12 / sigma22) * (y2_vec - mu2)
  sd_cond  <- sqrt(1 - (sigma12^2 / sigma22))
  
  # Standardized Bounds
  a_std <- (alpha_ext[y1_vec] - mu_cond) / sd_cond
  b_std <- (alpha_ext[y1_vec + 1] - mu_cond) / sd_cond
  pdf_a <- stats::dnorm(a_std); pdf_b <- stats::dnorm(b_std)
  cdf_a <- stats::pnorm(a_std); cdf_b <- stats::pnorm(b_std)
  idx1 <- which(y1_vec == 1); a_std[idx1] <- 0
  idxmax <- which(y1_vec == max(y1_vec)); b_std[idxmax] <- 0
  cdf_diff <- pmax(1e-16, cdf_b - cdf_a)
  
  # Partial Derivatives
  # d_loglik / d_mu2
  d_mu2 <- ((pdf_b - pdf_a) / cdf_diff) * (sigma12 / (sigma22 * sd_cond)) +
    (y2_vec - mu2) / sigma22
  # d_loglik / d_gamma
  score_gamma <- as.vector(t(X_mat_m2) %*% d_mu2)
  
  # d_loglik / d_sigma12
  da_sigma12 <- -((y2_vec - mu2) / sigma22 - (sigma12 * a_std) / (sigma22 * sd_cond)) / sd_cond
  db_sigma12 <- -((y2_vec - mu2) / sigma22 - (sigma12 * b_std) / (sigma22 * sd_cond)) / sd_cond
  d_sigma12 <- (pdf_b * db_sigma12 - pdf_a * da_sigma12) / cdf_diff
  score_sigma12 <- sum(d_sigma12)
  
  # d_loglik / d_sigma22
  db_sigma22 <- ((sigma12 * (y2_vec - mu2)) / sigma22^2 - (sigma12^2 * b_std) / (2 * sigma22^2 * sd_cond)) / sd_cond
  da_sigma22 <- ((sigma12 * (y2_vec - mu2)) / sigma22^2 - (sigma12^2 * a_std) / (2 * sigma22^2 * sd_cond)) / sd_cond
  d_sigma22 <- -1/(2 * sigma22) + ((y2_vec - mu2)^2) / (2 * sigma22^2) +
    (pdf_b * db_sigma22 - pdf_a * da_sigma22) / cdf_diff
  score_sigma22 <- sum(d_sigma22)
  
  return(c(score_gamma, score_sigma12, score_sigma22))
}

compute_ey1_star_s0 <- function(theta, y0, X0, x_support, x_name) {
  p_covs <- ncol(X0)
  beta   <- theta[1:p_covs]
  alpha  <- theta[(p_covs + 1):length(theta)]
  
  # Pre-calculate the linear predictor for Z (excluding X)
  x_col_idx <- which(colnames(X0) == x_name)
  beta_no_x <- beta[-x_col_idx]
  X0_no_x   <- X0[, -x_col_idx, drop = FALSE]
  # Constant part of mu for each subject: mu_fixed = Z * beta_z
  mu_fixed <- as.vector(X0_no_x %*% beta_no_x)
  
  # Vectorized calculation of mu for all support points v
  # mu_iv = mu_fixed_i + beta_x * x_v
  beta_x <- beta[x_col_idx]
  mu_mat <- outer(mu_fixed, beta_x * x_support, "+")  # (n0, d_size)
  
  # Calculate truncated normal mean
  alpha_ext <- c(-Inf, alpha, Inf)
  a_mat <- alpha_ext[y0] - mu_mat
  b_mat <- alpha_ext[y0 + 1] - mu_mat
  
  ey1_star_mat <- mu_mat - 
    (stats::dnorm(b_mat) - stats::dnorm(a_mat)) / pmax(1e-16, stats::pnorm(b_mat) - stats::pnorm(a_mat))
  
  return(ey1_star_mat)
}

compute_ey1_star_s1 <- function(theta1, theta2, y1, y2, Xmat_m1, Xmat_m2) {
  p_covs1 <- ncol(Xmat_m1)
  beta    <- theta1[1:p_covs1]
  alpha   <- theta1[(p_covs1 + 1):length(theta1)]
  alpha_ext <- c(-Inf, alpha, Inf)
  
  p_covs2 <- ncol(Xmat_m2)
  gamma   <- theta2[1:p_covs2]
  sigma12 <- theta2[p_covs2 + 1]
  sigma22 <- theta2[p_covs2 + 2]
  if (length(theta2) != p_covs2 + 2) {
    stop(paste0("Dimension mismatch in 'theta2': expected ", p_covs2 + 2, 
                " elements (", p_covs2, " for gamma + 2 for sigma12/sigma22), but found ", 
                length(theta2), "."))
  }
  
  mu1 <- as.vector(Xmat_m1 %*% beta)
  mu2 <- as.vector(Xmat_m2 %*% gamma)
  
  mu_cond <- mu1 + (sigma12 / sigma22) * (y2 - mu2)
  sd_cond <- sqrt(1 - (sigma12^2 / sigma22))
  
  # Calculate E[Y2* | Y2, Y1, X, Z1, Z2] using Truncated Normal Mean
  a <- alpha_ext[y1]
  b <- alpha_ext[y1 + 1]
  idx1 <- which(y1 == 1)
  idxmax <- which(y1 == max(y1))
  a_std <- (a - mu_cond) / sd_cond
  b_std <- (b - mu_cond) / sd_cond
  varphi_a <- stats::dnorm(a_std); varphi_b <- stats::dnorm(b_std)
  Phi_a <- stats::pnorm(a_std); Phi_b <- stats::pnorm(b_std)
  a_std[idx1] <- 0; b_std[idxmax] <- 0
  tmp1 <- (varphi_b - varphi_a) / pmax(1e-16, Phi_b - Phi_a)
  tmp2 <- (b_std * varphi_b - a_std * varphi_a) / pmax(1e-16, Phi_b - Phi_a)
  Ey1_star <- mu_cond - sd_cond * tmp1
  Vary1_star <- sd_cond^2 * (1 - tmp2 - tmp1^2)
  
  return(list(E = Ey1_star, Var = Vary1_star))
}

weighted_alpha_nll <- function(alpha, beta, y, X, weights) {
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

weighted_alpha_grad <- function(alpha, beta, y, X, weights) {
  if (any(diff(alpha) <= 0)) return(rep(0, length(alpha)))
  
  K_minus_1 <- length(alpha)
  mu <- as.vector(X %*% beta)
  alpha_ext <- c(-Inf, alpha, Inf)
  z_u <- alpha_ext[y + 1] - mu
  z_l <- alpha_ext[y] - mu
  phi_u <- dnorm(z_u)
  phi_l <- dnorm(z_l)
  p_y <- pmax(pnorm(z_u) - pnorm(z_l), 1e-16)  # Numerical Stability
  
  # Gradient for Alpha
  grad_alpha <- numeric(K_minus_1)
  for (j in 1:K_minus_1) {
    # Indicator if alpha_j is upper bound vs lower bound
    d_alpha_j <- (y == j) * (phi_u / p_y) - (y == j + 1) * (phi_l / p_y)
    grad_alpha[j] <- -sum(weights * d_alpha_j)
  }
  
  return(grad_alpha)
}

update_theta1_joint <- function(theta1, w_iv, y1_s0, y1_s1, Ey1_star_s0, Ey1_star_s1, Xmat_s0, Xmat_s1, x_support, x_name) {
  n1 <- length(y1_s1)
  n0 <- length(y1_s0)
  d_size <- length(x_support)
  
  # Prepare S=1 data: Weights are simply 1
  
  # Prepare S=0 data: each person in S=0 for every support point x_v
  y0_long <- rep(y1_s0, each = d_size)
  weights0_long <- as.vector(t(w_iv)) # Row-by-row weights from E-step
  
  X0_long <- Xmat_s0[rep(1:n0, each = d_size), , drop = FALSE]
  X0_long[, x_name] <- rep(x_support, n0)
  
  # Combine S=1 and S=0
  y_combined <- c(y1_s1, y0_long)
  Ey1_star_combined <- c(Ey1_star_s1, as.vector(t(Ey1_star_s0)))
  X_combined <- rbind(Xmat_s1, X0_long)
  weights_combined <- c(rep(1, n1), weights0_long)
  
  # Update beta via WLS
  beta <- as.vector(solve(t(X_combined) %*% (weights_combined * X_combined)) %*% (t(X_combined) %*% (weights_combined * Ey1_star_combined)))
  alpha <- theta1[(length(beta) + 1):length(theta1)]
  
  # Optimize
  fit <- optim(
    par = alpha,
    fn = weighted_alpha_nll,
    gr = weighted_alpha_grad,
    beta = beta,
    y = y_combined,
    X = X_combined,
    weights = weights_combined,
    method = "BFGS"
  )
  
  return(c(beta, fit$par))
}

update_theta2_joint <- function(theta1, theta2, y1, y2, Xmat_m1, Xmat_m2, Ey1_star, Vary1_star) {
  p_covs1 <- ncol(Xmat_m1)
  beta    <- theta1[1:p_covs1]
  alpha   <- theta1[(p_covs1 + 1):length(theta1)]
  alpha_ext <- c(-Inf, alpha, Inf)
  mu1 <- as.vector(Xmat_m1 %*% beta)
  p_covs2 <- ncol(Xmat_m2)
  n1 <- length(y1)
  # Residual of latent Y1* from its mean
  latent_resid <- Ey1_star - mu1
  # Update theta and rho via regression: Y2 = mu2 + rho*(Y* - mu_cond)
  # We use (Y* - mu_cond) as a covariate to account for the correlation
  Xaug_m2 <- cbind(Xmat_m2, latent_resid)
  XtX <- t(Xaug_m2) %*% Xaug_m2
  XtX[p_covs2 + 1, p_covs2 + 1] <- sum(Vary1_star + latent_resid^2)
  Xty <- t(Xaug_m2) %*% y2
  theta2_new <- as.vector(solve(XtX) %*% Xty)
  gamma_curr <- theta2_new[1:p_covs2]
  sigma12_curr <- theta2_new[p_covs2 + 1]
  sigma22_curr <- sum((y2 - as.vector(Xaug_m2 %*% theta2_new))^2 +
                        sigma12_curr^2 * Vary1_star) / n1 + sigma12_curr^2
  
  return(c(gamma_curr, sigma12_curr, sigma22_curr))
}
