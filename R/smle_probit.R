#' SMLE under Probit Link
#'
#' @param formula Formula for the ordinal outcome Y.
#' @param data Data frame. NA allowed only in \code{x_name} for S=0 subjects.
#' @param Bbasis Basis matrix for the spline basis.
#' @param x_name Name of the expensive covariate X.
#' @param theta_init Optional initial vector for theta (beta, cutpoints).
#' @param se_calc Logical; whether to compute standard errors. Defaults to TRUE.
#' @param max_iter Maximum number of EM iterations. Defaults to 500.
#' @param tol Convergence tolerance for the optimizer.
#'
#' @return A list containing estimates and convergence info.
#' @export
smle_probit <- function(formula, data, Bbasis, x_name,
                        theta_init = NULL, se_calc = TRUE, max_iter = 500, tol = 1e-6) {
  # -- 0. Validate ----
  if (any(is.na(Bbasis))) {
    stop("The 'Bbasis' matrix contains missing values. B-spline basis must be fully computed.")
  }
  cols2check <- setdiff(all.vars(formula), x_name)
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
  mf <- stats::model.frame(formula, data, na.action = na.pass)
  yvec <- as.numeric(stats::model.response(mf))
  Xmat <- stats::model.matrix(formula, mf)[, -1, drop = FALSE]  # Remove intercept
  # Split by selection indicator S
  yvec_s1 <- yvec[n1_idx]; yvec_s0 <- yvec[n0_idx]
  Xmat_s1 <- Xmat[n1_idx, , drop = FALSE]; Xmat_s0 <- Xmat[n0_idx, , drop = FALSE]
  Bbasis_s1 <- Bbasis[n1_idx, , drop = FALSE]; Bbasis_s0 <- Bbasis[n0_idx, , drop = FALSE]
  # Support points for X: all unique values observed in Phase 2
  x_name <- colnames(Xmat)[grepl(x_name, colnames(Xmat), fixed = TRUE)]
  x_support <- unique(Xmat_s1[, x_name, drop = FALSE])

  # -- 2. Initialize Parameters ----
  # theta: (beta, cutpoints)
  use_defaults <- TRUE
  if (!is.null(theta_init)) {
    theta_len <- ncol(Xmat) + (length(unique(yvec)) - 1)
    if (length(theta_init) == theta_len) {
      use_defaults <- FALSE
      theta_curr <- theta_init
    } else {
      warning(paste("theta_init must have length", theta_len))
    }
  }
  if (use_defaults) {
    fit_init <- MASS::polr(formula, data, method = "probit")
    theta_curr <- as.vector(c(fit_init$coefficients, fit_init$zeta))
  }

  # p_vl: normalized sieve coefficients
  s1_keys <- do.call(paste, as.data.frame(Xmat_s1[, x_name, drop = FALSE]))
  support_keys <- do.call(paste, as.data.frame(x_support))
  p_vl_num_init <- t(outer(s1_keys, support_keys, "==")) %*% Bbasis_s1  # (d, s_n)
  p_vl_curr <- sweep(p_vl_num_init, 2, pmax(1e-16, colSums(p_vl_num_init)), FUN = "/")

  # -- 3. EM Loop ----
  converged <- FALSE
  for (iter in 1:max_iter) {
    theta_old <- theta_curr

    ## -- E-STEP ----
    # For S=0, compute w_iv = E[I_iv | Y_i, Z_i, theta(m), p_vl(m)]
    prob_y_given_xv_s0 <- compute_py_given_xv_s0(theta_curr, yvec_s0, Xmat_s0, x_support, x_name)  # (n0, d)
    w_iv <- compute_w_iv_s0(prob_y_given_xv_s0, p_vl_curr, Bbasis_s0)  # (n0, d)

    # -- M-STEP ----
    # Update theta(m+1)
    theta_curr <- update_theta_smle(theta_curr, w_iv, yvec_s1, yvec_s0, Xmat_s1, Xmat_s0, x_support, x_name)
    # Update p_vl(m+1)
    p_vl_curr <- update_p_vl_cpp(prob_y_given_xv_s0, Bbasis_s0, p_vl_curr, p_vl_num_init, max_iter, tol)

    # Check convergence
    if (max(abs(theta_curr - theta_old)) < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning("EM algorithm reached max_iter without converging to tolerance.")
  }

  # -- 4. Standard Error Estimation ----
  se <- NULL
  vcov <- NULL
  if (se_calc) {
    se_results <- estimate_se_smle(
      theta = theta_curr, p_vl = p_vl_curr, p_vl_s1 = p_vl_num_init,
      yvec_s1, yvec_s0, Xmat_s1, Xmat_s0, Bbasis_s0, x_support, x_name
    )
    se <- se_results$se
    vcov <- se_results$vcov
  }

  return(list(est = theta_curr, se = se, vcov = vcov, p_vl = p_vl_curr, iterations = iter))
}
