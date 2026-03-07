#' EM Algorithm for Joint Modeling (Missing X and Y2)
#'
#' @param formula1 Formula for the ordinal outcome Y1 (Model 1).
#' @param formula2 Formula for the continuous outcome Y2 (Model 2).
#' @param data Data frame. NA allowed only in \code{x_name} and Y2 for S=0 subjects.
#' @param Bbasis Basis matrix for B-splines of Z1.
#' @param x_name Name of the expensive covariate X.
#' @param theta1_init Optional initial vector for theta1 (beta, cutpoints)
#' @param theta2_init Optional initial vector for theta2 (gamma, sigma12, sigma22).
#' @param se_calc Logical; whether to compute standard errors. Defaults to TRUE.
#' @param max_iter Maximum number of EM iterations. Defaults to 500.
#' @param tol Convergence tolerance for the optimizer.
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
  n1_idx <- which(s_indicator == 1); n0_idx <- which(s_indicator == 0)

  # -- 1. Construct Data Matrices ----
  # Model 1 (formula1)
  mf1 <- stats::model.frame(formula1, data, na.action = na.pass)
  y1vec <- as.numeric(stats::model.response(mf1))
  Xmat_m1 <- stats::model.matrix(formula1, mf1)[, -1, drop = FALSE]  # Remove intercept
  # Split by selection indicator S
  y1vec_s1 <- y1vec[n1_idx]; y1vec_s0 <- y1vec[n0_idx]
  Xmat_m1_s1 <- Xmat_m1[n1_idx, , drop = FALSE]; Xmat_m1_s0 <- Xmat_m1[n0_idx, , drop = FALSE]
  Bbasis_s1 <- Bbasis[n1_idx, , drop = FALSE]; Bbasis_s0 <- Bbasis[n0_idx, , drop = FALSE]
  # Support points for X: all unique values observed in Phase 2
  x_name <- colnames(Xmat)[grepl(x_name, colnames(Xmat), fixed = TRUE)]
  x_support <- unique(Xmat_s1[, x_name, drop = FALSE])

  # Model 2 (formula2)
  mf2 <- stats::model.frame(formula2, data)
  y2vec <- as.numeric(stats::model.response(mf2))
  Xmat_m2 <- stats::model.matrix(formula2, mf2)  # Include intercept
  pcovs_m2 <- ncol(Xmat_m2)

  # -- 2. Initialize Parameters ----
  # theta1: (beta, cutpoints)
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
    fit_init <- MASS::polr(formula1, data, method = "probit")
    theta1_curr <- as.vector(c(fit_init$coefficients, fit_init$zeta))
  }

  # theta2: (gamma, sigma12, sigma22)
  use_defaults <- TRUE
  if (!is.null(theta2_init)) {
    theta2_len <- pcovs_m2 + 2
    if (length(theta2_init) == theta2_len) {
      use_defaults <- FALSE
      theta2_curr <- theta2_init
      gamma_curr <- theta2_curr[1:pcovs_m2]
      sigma12_curr <- theta2_curr[pcovs_m2 + 1]
      sigma22_curr <- theta2_curr[pcovs_m2 + 2]
    } else {
      warning(paste("theta2_init must have length", theta2_len))
    }
  }
  if (use_defaults) {
    # Initial OLS guess for theta2
    fit_init <- stats::lm(formula2, data)
    gamma_curr <- as.vector(stats::coef(fit_init))
    sigma22_curr <- stats::var(stats::residuals(fit_init))
    sigma12_curr <- 0.1  # Initial correlation guess
    theta2_curr <- c(gamma_curr, sigma12_curr, sigma22_curr)
  }

  # p_vl: normalized sieve coefficients
  s1_keys <- do.call(paste, as.data.frame(Xmat_s1[, x_name, drop = FALSE]))
  support_keys <- do.call(paste, as.data.frame(x_support))
  p_vl_num_init <- t(outer(s1_keys, support_keys, "==")) %*% Bbasis_s1  # (d, s_n)
  p_vl_curr <- sweep(p_vl_num_init, 2, pmax(1e-16, colSums(p_vl_num_init)), FUN = "/")

  # -- 3. EM Loop ----
  converged <- FALSE
  for (iter in 1:max_iter) {
    theta_old <- c(theta1_curr, gamma_curr)

    ## -- E-STEP ----
    # For S=0, compute w_v = E[Iv | Y1, Z, theta1(m), p_vl(m)]
    prob_y1_given_xv_s0 <- compute_py_given_xv_s0(theta1_curr, y1vec_s0, Xmat_m1_s0, x_support, x_name)  # (n0, d)
    w_iv <- compute_w_iv_s0(prob_y1_given_xv_s0, p_vl_curr, Bbasis_s0)  # (n0, d)
    # For S=0, compute E[Y1* | x_v, Z, Y1, theta1(m)]
    e_y1star_s0 <- compute_y1star_s0(theta1_curr, y1vec_s0, Xmat_m1_s0, x_support, x_name)  # (n0, d)
    # For S=1, compute moments of [Y1* | Y2, X, Z, Y1, theta1(m)]
    y1star_s1 <- compute_y1star_s1(theta1_curr, theta2_curr, y1vec_s1, y2vec, Xmat_m1_s1, Xmat_m2)
    e_y1star_s1 <- y1star_s1$mean
    var_y1star_s1 <- y1star_s1$var

    # -- M-STEP ----
    # Update p_vl(m+1)
    p_vl_curr <- update_p_vl_cpp(prob_y1_given_xv_s0, Bbasis_s0, p_vl_curr, p_vl_num_init, max_iter, tol)

    # Update theta2(m+1)
    theta2_curr <- update_theta2(theta1_curr, theta2_curr, y1vec_s1, y2vec, Xmat_m1_s1, Xmat_m2, e_y1star_s1, var_y1star_s1)
    gamma_curr <- theta2_curr[1:pcovs_m2]

    # Update theta1(m+1)
    theta1_curr <- update_theta1(theta1_curr, y1vec_s1, y1vec_s0, Xmat_m1_s1, Xmat_m1_s0,
                                 w_iv, e_y1star_s1, e_y1star_s0, x_support, x_name)

    # Check convergence
    if (max(abs(c(theta1_curr, gamma_curr) - theta_old)) < tol) {
      converged <- TRUE
      break
    }
  }
  sigma12_curr <- theta2_curr[pcovs_m2 + 1]
  sigma22_curr <- theta2_curr[pcovs_m2 + 2]

  if (!converged) {
    warning("EM algorithm reached max_iter without converging to tolerance.")
  }

  return(list(theta1 = theta1_curr, gamma = gamma_curr, sigma12 = sigma12_curr, sigma22 = sigma22_curr, iterations = iter))
}
