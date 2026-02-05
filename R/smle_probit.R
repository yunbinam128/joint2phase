#' SMLE under probit
#'
#' @param formula A model formula for the probit specification.
#' @param data A data frame containing the variables in the formula.
#' @param B_basis Basis matrix for the spline basis.
#' @param x_name Character name of the expensive covariate (missingness).
#' @param se_calc Logical; whether to compute standard errors.
#' @param max_iter Maximum number of optimization iterations.
#' @param tol Convergence tolerance for the optimizer.
#'
#' @return A list containing estimates and convergence info.
#' @export
smle_probit <- function(formula, data, B_basis, x_name, se_calc = TRUE, max_iter = 100, tol = 1e-6) {
  # -- 0. Validate ----
  if (any(is.na(B_basis))) {
    stop("The 'B_basis' matrix contains missing values. B-spline basis must be fully computed.")
  }
  cols2check <- setdiff(all.vars(formula), x_name)
  na_counts <- colSums(is.na(data[, cols2check, drop = FALSE]))
  if (any(na_counts > 0)) {
    missing_cols <- names(na_counts[na_counts > 0])
    stop(paste("Missing values detected in the following columns:",
               paste(missing_cols, collapse = ", "),
               "- Only the column specified in 'x_name' can have missingness."))
  }
  if (nrow(B_basis) != nrow(data)) {
    stop("B_basis must have the same number of rows as data.")
  }

  # -- 1. Construct Data Matrices ----
  mf <- model.frame(formula, data, na.action = na.pass)
  y_vec <- as.numeric(model.response(mf))
  X_mat <- model.matrix(formula, mf)[, -1, drop = FALSE]  # remove the intercept (column 1)

  # Identify selection indicator (assuming S=1 if x_name is NOT NA, and S=0 if x_name is NA)
  s_indicator <- as.numeric(!is.na(data[[x_name]]))
  # Split by selection indicator S
  y1 <- y_vec[s_indicator == 1]; y0 <- y_vec[s_indicator == 0]
  X1 <- X_mat[s_indicator == 1, , drop = FALSE]; X0 <- X_mat[s_indicator == 0, , drop = FALSE]
  B1 <- B_basis[s_indicator == 1, , drop = FALSE]; B0 <- B_basis[s_indicator == 0, , drop = FALSE]

  # -- 2. Setup Support ----
  # Support points for X: all unique values observed in Phase 2
  x_support <- sort(unique(data[[x_name]][s_indicator == 1]))

  # -- 3. Initialize Parameters ----
  # Initialize p betas at 0 and K-1 alphas spread from -1 to 1
  theta_curr <- c(numeric(ncol(X_mat)), seq(-1, 1, length.out = length(unique(y_vec)) - 1))
  # Initialize p_vl
  p_vl_num_init <- t(outer(data[[x_name]][s_indicator == 1], x_support, "==")) %*% B1  # (d_size, s_n): numerator for initial p_vl
  p_vl_curr <- sweep(p_vl_num_init, 2, colSums(p_vl_num_init) + 1e-12, FUN = "/")      # Column-wise normalization

  # -- 4. EM Algorithm Loop ----
  for (iter in 1:max_iter) {
    theta_old <- theta_curr

    ## -- E-STEP: Compute expectations for S=0 ----
    # Calculate P(Y_i | x_v, Z_i; theta) for all i in S=0 and all v in support
    prob_y0_v <- compute_prob_matrix(theta_curr, y0, X0, x_support, x_name)  # (n0, d_size)
    # Profile out p_vl
    p_vl_curr <- update_p_vl_cpp(prob_y0_v, B0, p_vl_curr, p_vl_num_init, max_iter, tol)
    # Calculate weights
    # Marginal P(X=v | Z) = sum_l { B_l(Z) * p_vl }
    prob_x_v_given_z0 <- B0 %*% t(p_vl_curr)
    # Compute q_iv weights: P(X=v | Y, Z; theta, p_vl)
    joint_y_x_v <- prob_y0_v * prob_x_v_given_z0
    denom_y0    <- rowSums(joint_y_x_v) + 1e-12
    q_iv        <- joint_y_x_v / denom_y0

    # -- M-STEP: Update parameters ----
    # 1. Update theta (theta = argmax weighted log-likelihood)
    # weights=1 for S=1 and weights=q_iv for S=0 (each subject i is repeated d times for S=0)
    theta_curr <- update_theta_weighted(theta_curr, q_iv, y1, X1, y0, X0, x_support, x_name)

    # Check Convergence
    if (max(abs(theta_curr - theta_old)) < tol) break
  }

  # -- 5. Standard Error Estimation via profile log-likelihood (se_calc = TRUE) ----
  se <- NULL
  vcov <- NULL
  if (se_calc) {
    se_results <- estimate_smle_se(
      theta_hat = theta_curr,
      p_vl_hat = p_vl_curr,
      p_vl_R1 = p_vl_num_init,
      y1, X1, y0, X0, B0, x_support, x_name
    )
    se <- se_results$se
    vcov <- se_results$vcov
  }

  return(list(
    theta = theta_curr,
    p_vl = p_vl_curr,
    se = se,
    vcov = vcov,
    iterations = iter
  ))
}
