#' EM Algorithm for Joint Modeling - assuming both X and Y2 are available only for subjects selected in Phase 2
#'
#' @param formula1 Formula for the ordinal outcome Y1 (probit).
#' @param formula2 Formula for the continuous outcome Y2 (Model 2).
#' @param data Full dataset containing Phase 1 and Phase 2.
#' @param Bbasis Basis matrix for B-splines of Z1.
#' @param x_name Name of the expensive covariate X.
#' @param theta1_init Optional numeric vector of initial values for theta1.
#' @param theta2_init Optional numeric vector of initial values for theta2.
#' @param se_calc Logical; if TRUE, computes SEs via profile likelihood.
#' @param max_iter Maximum number of EM iterations. Defaults to 500.
#' @param tol Convergence tolerance.
#'
#' @return A list containing estimates and convergence info.
#' @export
em_joint <- function(formula1, formula2, data, Bbasis, x_name,
                     theta1_init = NULL, theta2_init = NULL, se_calc = TRUE, max_iter = 500, tol = 1e-6) {

  # -- 0. Validate ----
  if (any(is.na(Bbasis))) {
    stop("The 'Bbasis' matrix contains missing values. B-spline basis must be fully computed.")
  }
  cols2check <- setdiff(all.vars(formula1), x_name)
  na_counts <- colSums(is.na(data[, cols2check, drop = FALSE]))
  if (any(na_counts > 0)) {
    missing_cols <- names(na_counts[na_counts > 0])
    stop(paste("Missing values detected in the following columns:",
               paste(missing_cols, collapse = ", "),
               "- Only the column specified in 'x_name' can have missingness."))
  }
  if (nrow(Bbasis) != nrow(data)) {
    stop("Bbasis must have the same number of rows as data.")
  }

  # Identify selection indicator (assuming S=1 if x_name is NOT NA, and S=0 if x_name is NA)
  s_indicator <- as.numeric(!is.na(data[[x_name]]))
  n1_idx <- which(s_indicator == 1)
  n0_idx <- which(s_indicator == 0)
  n0 <- length(n0_idx)
  x_support <- sort(unique(data[[x_name]][n1_idx]))
  d_size <- length(x_support)
  data_p2 <- data[n1_idx, ]

  # -- 1. Construct Data Matrices ----
  # Model 1 (formula1)
  mf1 <- model.frame(formula1, data, na.action = na.pass)
  y1vec <- as.numeric(model.response(mf1))
  Xmat_m1 <- model.matrix(formula1, mf1)[, -1, drop = FALSE]  # remove the intercept (column 1)
  y1_s1 <- y1vec[n1_idx]; y1_s0 <- y1vec[n0_idx]
  Xmat_m1_s1 <- Xmat_m1[n1_idx, , drop = FALSE]; Xmat_m1_s0 <- Xmat_m1[n0_idx, , drop = FALSE]
  B1 <- Bbasis[n1_idx, , drop = FALSE]; B0 <- Bbasis[n0_idx, , drop = FALSE]
  sn_size <- ncol(B0)

  # Model 2 (formula2)
  mf2 <- model.frame(formula2, data_p2)
  y2_vec <- as.numeric(model.response(mf2))
  Xmat_m2 <- model.matrix(formula2, mf2) # Includes intercept
  p_covs2 <- ncol(Xmat_m2)

  # -- 2. Initialize Parameters ----
  # theta1: (beta, alpha/cutpoints)
  use_defaults <- TRUE
  if (!is.null(theta1_init)) {
    theta1_len <- ncol(Xmat_m1) + (length(unique(y1vec)) - 1)
    if (length(theta1_init) == theta1_len) {
      use_defaults <- FALSE
      theta1_curr <- theta1_init
    } else {
      warning(paste("theta1_init must have length", theta1_len))
    }
  }
  if (use_defaults) {
    # Initial probit guess for theta1
    fit_init <- MASS::polr(formula1, data_p2, method = "probit")
    theta1_curr <- as.vector(c(fit_init$coefficients, fit_init$zeta))
  }

  # theta2: (gamma, sigma12, sigma22)
  if (!is.null(theta2_init)) {
    theta2_len <- p_covs2 + 2
    if (length(theta2_init) == theta2_len) {
      use_defaults <- FALSE
      theta2_curr <- theta2_init
      gamma_curr <- theta2_curr[1:p_covs2]
      sigma12_curr <- theta2_curr[p_covs2 + 1]
      sigma22_curr <- theta2_curr[p_covs2 + 2]
    } else {
      warning(paste("theta2_init must have length", theta2_len))
    }
  }
  if (use_defaults) {
    # Initial OLS guess for theta2
    fit_init <- stats::lm(formula2, data_p2)
    gamma_curr <- as.vector(stats::coef(fit_init))
    sigma22_curr <- stats::var(stats::residuals(fit_init))
    sigma12_curr <- 0.1  # Initial correlation guess
    theta2_curr <- c(gamma_curr, sigma12_curr, sigma22_curr)
  }
  # p_vl: normalized sieve coefficients
  p_vl_num_init <- t(outer(data[[x_name]][n1_idx], x_support, "==")) %*% B1  # (d, s_n): numerator for initial p_vl
  p_vl_curr <- sweep(p_vl_num_init, 2, pmax(colSums(p_vl_num_init), 1e-16), FUN = "/")  # Column-wise normalization

  # -- 3. EM Loop ----
  for (iter in 1:max_iter) {
    theta_old <- c(theta1_curr, gamma_curr)

    ## -- E-STEP ----
    # for S=0 (i=n1+1,...,n)
    ## Compute w_iv
    ### P(Y1 | x_v, Z)
    prob_y1_given_xz <- compute_prob_matrix(theta1_curr, y1_s0, Xmat_m1_s0, x_support, x_name)  # (n0, d)
    ### P(X = x_v | Z) = sum_l { p_vl * B_l(Z) }
    prob_x_given_z <- B0 %*% t(p_vl_curr)  # (n0, d)
    w_iv_num <- prob_y1_given_xz * prob_x_given_z  # (n0, d)
    w_denom <- pmax(1e-16, rowSums(w_iv_num))  # (n0)
    w_iv <- w_iv_num / w_denom  # (n0, d)
    ## Compute conditional expectation of Y1* | Y1, X, Z
    Ey1_star_s0 <- compute_ey1_star_s0(theta1_curr, y1_s0, Xmat_m1_s0, x_support, x_name)  # (n0, d)

    # for S=1 (i=1,...,n1)
    ## Compute conditional moments of Y1* | Y2, Y1, X, Z
    y1_star_moments_s1 <- compute_ey1_star_s1(theta1_curr, theta2_curr, y1_s1, y2_vec, Xmat_m1_s1, Xmat_m2)
    Ey1_star_s1 <- y1_star_moments_s1$E; Vary1_star_s1 <- y1_star_moments_s1$Var

    ## -- M-STEP ----
    # Update {p_vl}
    p_vl_curr <- update_p_vl_cpp(prob_y1_given_xz, B0, p_vl_curr, p_vl_num_init, max_iter, tol)

    # Update theta2
    theta2_curr <- update_theta2_joint(theta1_curr, theta2_curr, y1_s1, y2_vec, Xmat_m1_s1, Xmat_m2, Ey1_star_s1, Vary1_star_s1)
    gamma_curr <- theta2_curr[1:p_covs2]

    # Update theta1
    theta1_curr <- update_theta1_joint(theta1_curr, w_iv, y1_s0, y1_s1, Ey1_star_s0, Ey1_star_s1, Xmat_m1_s0, Xmat_m1_s1, x_support, x_name)

    if (max(abs(c(theta1_curr, gamma_curr) - theta_old)) < tol) break
  }

  # -- 4. Standard Error Estimation ----
  se_theta1 <- NULL
  se_theta2 <- NULL
  vcov_mat <- NULL
  if (se_calc) {
    se_res <- estimate_em_joint_se(theta1_curr, theta2_curr, p_vl_curr, p_vl_num_init,
                                   y1_s0, y1_s1, y2_vec, Xmat_m1_s0, Xmat_m1_s1, Xmat_m2,
                                   B0, x_support, x_name)
    p_covs1 <- length(theta1_curr)
    se_theta1 <- se_res$se[1:p_covs1]
    se_theta2 <- se_res$se[(p_covs1 + 1):length(se_res$se)]
    vcov_mat <- se_res$vcov
  }

  return(list(theta1 = theta1_curr, se_theta1 = se_theta1,
              gamma = theta2_curr[1:p_covs2],
              sigma12 = theta2_curr[p_covs2 + 1],
              sigma22 = theta2_curr[p_covs2 + 2],
              se_theta2 = se_theta2, iterations = iter))
}
