#' EM Algorithm for Conditional Modeling (Step 2 of Two-Step Estimator)
#'
#' @description
#' Estimates parameters for a continuous outcome (Y2) conditioned on a previously
#' estimated ordinal probit model (Y1), treating the latent variable Y1* as missing.
#'
#' @param formula Formula for the continuous outcome Y2 (Model 2).
#' @param data Data frame containing all variables in \code{formula} and \code{formula_cond}. No NA allowed.
#' @param formula_cond Formula for the ordinal outcome Y1 (Model 1).
#' @param theta_cond Vector of estimated parameters from Model 1 (betas, cutpoints).
#' @param vcov_cond Covariance matrix of \code{theta_cond}; required if \code{se_calc} is TRUE.
#' @param theta2_init Optional initial vector for theta2 (gamma, sigma12, sigma22).
#' @param se_calc Logical; whether to compute standard errors. Defaults to TRUE.
#' @param max_iter Maximum number of EM iterations. Defaults to 100.
#' @param tol Convergence tolerance for the optimizer.
#'
#' @return A list containing estimates and convergence info.
#' @export
em_conditional <- function(formula, data, formula_cond, theta_cond, vcov_cond = NULL,
                           theta2_init = NULL, se_calc = TRUE, max_iter = 100, tol = 1e-6) {
  # -- 0. Validate ----
  if (se_calc && is.null(vcov_cond)) {
    stop("vcov_cond must be provided when se_calc is TRUE to account for Model 1 uncertainty.")
  }

  cols2check <- intersect(unique(c(all.vars(formula), all.vars(formula_cond))), colnames(data))
  na_counts <- colSums(is.na(data[, cols2check, drop = FALSE]))
  if (any(na_counts > 0)) {
    missing_cols <- names(na_counts[na_counts > 0])
    stop(paste("Missing values detected in the following columns:",
               paste(missing_cols, collapse = ", ")))
  }

  # -- 1. Construct Data Matrices ----
  mf2 <- model.frame(formula, data)
  y2vec <- as.numeric(model.response(mf2))
  Xmat_m2 <- model.matrix(formula, mf2)
  pcovs_m2 <- ncol(Xmat_m2)

  # -- 2. Initialize Parameters ----
  # theta2: (gamma, sigma12, sigma22)
  use_defaults <- TRUE
  if (!is.null(theta2_init)) {
    theta2_len <- pcovs_m2 + 2
    if (length(theta2_init) == theta2_len) {
      use_defaults <- FALSE
      gamma_curr <- theta2_init[1:pcovs_m2]
      sigma12_curr <- theta2_init[pcovs_m2 + 1]
      sigma22_curr <- theta2_init[pcovs_m2 + 2]
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
    theta2_curr <- c(gamma_curr, sigma12_curr, sigma22_curr)
  }

  # theta1: (beta, cutpoints) - fixed in em_conditional()
  mf1 <- model.frame(formula_cond, data)
  y1vec <- as.numeric(model.response(mf1))
  n1 <- length(y1vec)
  Xmat_m1 <- model.matrix(formula_cond, mf1)[, -1, drop = FALSE]  # Remove intercept
  pcovs_m1 <- ncol(Xmat_m1)
  beta <- theta_cond[1:pcovs_m1]
  alpha <- theta_cond[(pcovs_m1+1):length(theta_cond)]
  alpha_ext <- c(-Inf, alpha, Inf)
  a <- alpha_ext[y1vec]; b <- alpha_ext[y1vec + 1]
  mu1 <- as.vector(Xmat_m1 %*% beta)

  # -- 3. EM Loop ----
  converged <- FALSE
  for (iter in 1:max_iter) {
    gamma_old <- gamma_curr

    ## -- E-STEP ----
    # Compute moments of [Y1* | Y2, X, Z, Y1, theta1(m)]
    y1star <- compute_y1star_s1(theta_cond, theta2_curr, y1vec, y2vec, Xmat_m1, Xmat_m2, mu1)
    e_y1star <- y1star$mean
    var_y1star <- y1star$var

    # -- M-STEP ----
    # Update theta2(m+1)
    theta2_curr <- update_theta2(theta_cond, theta2_curr, y1vec, y2vec, Xmat_m1, Xmat_m2, e_y1star, var_y1star)
    gamma_curr <- theta2_curr[1:pcovs_m2]

    # Check convergence
    if (max(abs(gamma_curr - gamma_old)) < tol) {
      converged <- TRUE
      break
    }
  }
  sigma12_curr <- theta2_curr[pcovs_m2 + 1]
  sigma22_curr <- theta2_curr[pcovs_m2 + 2]
  names(gamma_curr) <- colnames(Xmat_m2)

  if (!converged) {
    warning("EM algorithm reached max_iter without converging to tolerance.")
  }

  # -- 4. Standard Error Estimation ----
  se <- NULL
  vcov <- NULL
  if (se_calc) {
    se_results <- estimate_se_em_cond(
      theta2 = theta2_curr, y2vec, Xmat_m2, theta_cond, vcov_cond, y1vec, Xmat_m1, mu1, pcovs_m1, alpha_ext
    )
    se <- se_results$se[1:pcovs_m2]
    vcov <- se_results$vcov[1:pcovs_m2, 1:pcovs_m2]
    names(se) <- rownames(vcov) <- colnames(vcov) <- names(gamma_curr)
  }

  return(list(gamma = gamma_curr, sigma12 = sigma12_curr, sigma22 = sigma22_curr,
              se = se, vcov = vcov, iterations = iter))
}
