# -- Objective Function: Negative Ascertainment-Corrected Log-Likelihood ----
acml_probit_nll <- function(params, X_mat, y_vec, pi_perY) {
  K <- length(pi_perY)
  p <- ncol(X_mat)
  
  # Extract parameters
  beta <- params[1:p]
  alpha <- params[(p + 1):(p + K - 1)]
  if (any(diff(alpha) <= 0)) return(1e10)  # Constraint: Thresholds must be strictly increasing
  mu <- as.vector(X_mat %*% beta)  # Linear predictor
  
  # Calculate category probabilities P(Y=m | X)
  p_m <- matrix(0, nrow = nrow(X_mat), ncol = K)
  p_m[, 1] <- pnorm(alpha[1] - mu)
  if (K > 2) {
    for (m in 2:(K - 1)) {
      p_m[, m] <- pnorm(alpha[m] - mu) - pnorm(alpha[m - 1] - mu)
    }
  }
  p_m[, K] <- 1 - pnorm(alpha[K - 1] - mu)
  
  # Ascertainment-corrected likelihood components
  # Numerator: Probability of the observed outcome
  y_obs_idx <- as.numeric(y_vec)
  p_y_obs <- p_m[cbind(1:nrow(X_mat), y_obs_idx)]
  # Denominator: Weighted sum of probabilities (Probability of selection)
  denom <- as.vector(p_m %*% pi_perY)
  
  return(-sum(log(p_y_obs + 1e-10) - log(denom + 1e-10)))
}

# -- Gradient Function: Analytical Derivatives for BFGS ----
acml_probit_grad <- function(params, X_mat, y_vec, pi_perY) {
  K <- length(pi_perY)
  n <- nrow(X_mat)
  p <- ncol(X_mat)
  
  beta <- params[1:p]
  alpha <- params[(p + 1):(p + K - 1)]
  mu <- as.vector(X_mat %*% beta)
  y_obs <- as.numeric(y_vec)
  
  # Pre-calculate normal PDF (phi) and CDF (Phi) at each threshold
  phi_mat <- matrix(0, nrow = n, ncol = K + 1)
  Phi_mat <- matrix(0, nrow = n, ncol = K + 1)
  Phi_mat[, K + 1] <- 1  # Phi(Inf) = 1
  
  for(j in 1:(K-1)) {
    z <- alpha[j] - mu
    Phi_mat[, j+1] <- pnorm(z)
    phi_mat[, j+1] <- dnorm(z)
  }
  
  p_m <- Phi_mat[, 2:(K+1)] - Phi_mat[, 1:K]
  D <- as.vector(p_m %*% pi_perY)
  
  # --- Gradient for beta ---
  # d(p_im)/d(beta) = X * (phi_{m-1} - phi_m)
  dp_dbeta_obs <- X_mat * (phi_mat[cbind(1:n, y_obs)] - phi_mat[cbind(1:n, y_obs + 1)])
  
  dD_dbeta <- matrix(0, nrow = n, ncol = p)
  for(m in 1:K) {
    dD_dbeta <- dD_dbeta + pi_perY[m] * (X_mat * (phi_mat[, m] - phi_mat[, m+1]))
  }
  grad_beta <- -colSums(dp_dbeta_obs / p_m[cbind(1:n, y_obs)] - dD_dbeta / D)
  
  # --- Gradient for alpha ---
  grad_alpha <- numeric(K - 1)
  for(j in 1:(K-1)) {
    # Derivative of observed part w.r.t alpha_j
    d_obs_alpha_j <- (y_obs == j) * phi_mat[, j+1] - (y_obs == j+1) * phi_mat[, j+1]
    # Derivative of selection part w.r.t alpha_j
    term2_alpha <- (pi_perY[j] - pi_perY[j+1]) * phi_mat[, j+1] / D
    grad_alpha[j] <- -sum(d_obs_alpha_j / p_m[cbind(1:n, y_obs)] - term2_alpha)
  }
  
  return(c(grad_beta, grad_alpha))
}