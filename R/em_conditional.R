#' Estimates parameters for a continuous outcome model conditioned on a previously
#' estimated ordinal probit model. This function implements an EM algorithm where
#' the latent continuous variable \eqn{Y_1^*} is treated as missing data.
#'
#' @param formula A \code{formula} object for the continuous outcome (\eqn{Y_2}).
#' @param data A data frame containing the variables in the model.
#' @param formula_cond The \code{formula} used for the ordinal outcome model (\eqn{Y_1}).
#' @param theta_cond Vector of estimated parameters from the ordinal model (\eqn{\theta_1}).
#' @param vcov_cond Covariance matrix of \eqn{\theta_1}. Optional; only required if \code{se_calc = TRUE}.
#' @param theta2_init Optional numeric vector of initial values for theta2.
#' @param se_calc Logical; whether to compute standard errors for Model 2 parameters. Defaults to TRUE.
#' @param max_iter Maximum number of EM iterations. Defaults to 100.
#' @param tol Convergence tolerance.
#'
#' @return A list containing estimates and convergence info.
#'
#' @importFrom stats lm coef var residuals dnorm pnorm model.frame model.response model.matrix
#' @importFrom numDeriv hessian jacobian
#' @export
em_conditional <- function(formula, data, formula_cond, theta_cond, vcov_cond = NULL,
                           theta2_init = NULL, se_calc = TRUE, max_iter = 100, tol = 1e-6) {
  # -- 0. Validate ----
  if (se_calc && is.null(vcov_cond)) {
    stop("vcov_cond must be provided when se_calc is TRUE to account for Model 1 uncertainty.")
  }

  cols2check <- unique(c(all.vars(formula), all.vars(formula_cond)))
  na_counts <- colSums(is.na(data[, cols2check, drop = FALSE]))
  if (any(na_counts > 0)) {
    missing_cols <- names(na_counts[na_counts > 0])
    stop(paste("Missing values detected in the following columns:",
               paste(missing_cols, collapse = ", ")))
  }

  # -- 1. Setup Parameters from Model 1 ----
  mf1 <- model.frame(formula_cond, data)
  y1_vec <- as.numeric(model.response(mf1))
  n1 <- length(y1_vec)
  Xmat_m1 <- model.matrix(formula_cond, mf1)[, -1, drop = FALSE]  # remove the intercept (column 1)
  p_covs_m1 <- ncol(Xmat_m1)
  beta <- theta_cond[1:p_covs_m1]
  alpha <- theta_cond[(p_covs_m1 + 1):length(theta_cond)]
  alpha_ext <- c(-Inf, alpha, Inf)

  # -- 2. Initialize Model 2 Parameters ----
  mf2 <- model.frame(formula, data)
  y2_vec <- as.numeric(model.response(mf2))
  Xmat_m2 <- model.matrix(formula, mf2)
  p_covs_m2 <- ncol(Xmat_m2)
  use_defaults <- TRUE
  if (!is.null(theta2_init)) {
    theta2_len <- p_covs_m2 + 2
    if (length(theta2_init) == theta2_len) {
      use_defaults <- FALSE
      gamma_curr <- theta2_init[1:p_covs_m2]
      sigma12_curr <- theta2_init[p_covs_m2 + 1]
      sigma22_curr <- theta2_init[p_covs_m2 + 2]
    } else {
      warning(paste("theta2_init must have length", theta2_len))
    }
  }
  if (use_defaults) {
    # Initial OLS guess for theta2
    fit_init <- stats::lm(formula, data)
    gamma_curr <- as.vector(stats::coef(fit_init))
    sigma22_curr <- stats::var(stats::residuals(fit_init))
    sigma12_curr <- 0.1  # Initial correlation guess
  }

  # theta1: fixed in the EM algorithm
  # Truncation bounds for Y1* based on observed category Y1
  a <- alpha_ext[y1_vec]
  b <- alpha_ext[y1_vec + 1]
  # Linear predictor
  mu1 <- as.vector(Xmat_m1 %*% beta)
  # Initialize vectors for the calculation using truncated normal distribution
  idx1 <- which(y1_vec == 1)
  idxmax <- which(y1_vec == max(y1_vec))
  Ey1_star <- numeric(n1); Vary1_star <- numeric(n1)
  # -- 3. EM Loop ----
  for (iter in 1:max_iter) {
    gamma_old <- gamma_curr
    # Linear predictor
    mu2 <- as.vector(Xmat_m2 %*% gamma_curr)

    ## -- E-STEP: Conditional Expectation of Latent Y* ----
    # Distribution of Y1* | Y2, X, Z1, Z2 is Normal with:
    # Mean = mu1 + (sigma12/sigma22)*(Y2 - mu2)
    # Var  = 1 - (sigma12^2/sigma22)
    mu_tmp <- mu1 + (sigma12_curr / sigma22_curr) * (y2_vec - mu2)
    sigma_tmp <- sqrt(1 - (sigma12_curr^2 / sigma22_curr))

    # Calculate E[Y2* | Y2, Y1, X, Z1, Z2] using Truncated Normal Mean
    a_std <- (a - mu_tmp) / sigma_tmp
    b_std <- (b - mu_tmp) / sigma_tmp
    varphi_a <- stats::dnorm(a_std); varphi_b <- stats::dnorm(b_std)
    Phi_a <- stats::pnorm(a_std); Phi_b <- stats::pnorm(b_std)
    a_std[idx1] <- 0; b_std[idxmax] <- 0
    tmp1 <- (varphi_b - varphi_a) / pmax(1e-16, Phi_b - Phi_a)
    tmp2 <- (b_std * varphi_b - a_std * varphi_a) / pmax(1e-16, Phi_b - Phi_a)
    Ey1_star <- mu_tmp - sigma_tmp * tmp1
    Vary1_star <- sigma_tmp^2 * (1 - tmp2 - tmp1^2)

    # -- M-STEP: Update theta, sigma, rho ----
    # Residual of latent Y1* from its mean
    latent_resid <- Ey1_star - mu1

    # Update theta and rho via regression: Y2 = mu2 + rho*(Y* - mu_cond)
    # We use (Y* - mu_cond) as a covariate to account for the correlation
    Xaug_m2 <- cbind(Xmat_m2, latent_resid)
    XtX <- t(Xaug_m2) %*% Xaug_m2
    XtX[p_covs_m2 + 1, p_covs_m2 + 1] <- sum(Vary1_star + latent_resid^2)
    Xty <- t(Xaug_m2) %*% y2_vec
    theta2_new <- as.vector(solve(XtX) %*% Xty)
    gamma_curr <- theta2_new[1:p_covs_m2]
    sigma12_curr <- theta2_new[p_covs_m2 + 1]
    sigma22_curr <- sum((y2_vec - as.vector(Xaug_m2 %*% theta2_new))^2 +
                          sigma12_curr^2 * Vary1_star) / n1 + sigma12_curr^2

    # Convergence Check
    if (max(abs(gamma_curr - gamma_old)) < tol) break
  }

  # -- 4. Standard Error Estimation ----
  se <- NULL
  vcov_mat <- NULL
  if (se_calc) {
    # Calculate Hessian numerically for conditional variability
    theta2_hat <- c(gamma_curr, sigma12_curr, sigma22_curr)
    hess <- numDeriv::hessian(
      func = function(t) {
        obs_loglik(theta = t, p_covs_m2, y2_vec, y1_vec, Xmat_m2, mu1, alpha_ext)
      }, x = theta2_hat)
    # Estimate the conditional variability
    V22_inv <- try(solve(-hess), silent = TRUE)

    if (inherits(V22_inv, "try-error")) {
      warning("Hessian inversion failed; standard errors cannot be computed.")
    } else {
      # Calculate the Jacobian of the score function w.r.t theta1
      V21 <- numDeriv::jacobian(
        func = function(t1) {
          beta_tmp <- t1[1:p_covs_m1]
          alpha_tmp <- t1[(p_covs_m1 + 1):length(t1)]
          alpha_ext_tmp <- c(-Inf, alpha_tmp, Inf)
          mu1_tmp <- as.vector(Xmat_m1 %*% beta_tmp)
          score_theta2_function(theta2_hat, p_covs_m2, y2_vec, y1_vec, Xmat_m2, mu1_tmp, alpha_ext_tmp)
        },
        x = theta_cond)

      # Law of total variance
      vcov_mat <- V22_inv + (V22_inv %*% V21 %*% vcov_cond %*% t(V21) %*% V22_inv)
      se <- sqrt(diag(vcov_mat))
    }
  }

  return(list(
    gamma = gamma_curr,
    sigma12 = sigma12_curr,
    sigma22 = sigma22_curr,
    se = se,
    vcov_mat = vcov_mat,
    iterations = iter
  ))
}
