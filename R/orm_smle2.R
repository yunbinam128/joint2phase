#' SMLE for ordinal outcomes with envelope theorem SE
#'
#' Same estimator as \code{\link{orm_smle}}, but computes standard errors using
#' the envelope theorem: the analytical profile score is differentiated numerically
#' via \code{numDeriv::jacobian} (first-order), rather than differentiating the
#' profile log-likelihood via \code{numDeriv::hessian} (second-order).
#'
#' @inheritParams orm_smle
#' @param jacobian_method Passed to \code{numDeriv::jacobian}. \code{"Richardson"}
#'   (default) performs extrapolation across multiple small steps; \code{"simple"}
#'   is a forward difference with a single user-specified step.
#' @param jacobian_method_args Named list forwarded as \code{method.args} to
#'   \code{numDeriv::jacobian}. For \code{"simple"} this controls \code{eps}
#'   (step size). To mirror the PLL step, use
#'   \code{list(eps = length(data)^(-1/2))}.
#'
#' @return A list containing estimates, SE, and (if \code{se_calc=TRUE}) a
#'   \code{se_diagnostics} element with per-evaluation \code{p_vl_move},
#'   \code{score_mags}, \code{inner_iters}, \code{theta_perturb}, the raw
#'   (unsymmetrized) Hessian, and the numDeriv settings used.
#' @export
orm_smle2 <- function(formula, data, Bbasis, x_name, family = "probit",
                      theta_init = NULL,
                      se_calc = TRUE,
                      verbose = TRUE, max_iter = 500, tol = 1e-6,
                      se_max_iter = 1000, se_tol = 1e-8,
                      p_vl_init = NULL,
                      jacobian_method = "Richardson",
                      jacobian_method_args = list()) {
  # -- 0. Validate ----
  if (!(family %in% c("probit", "logistic"))) {
    stop("The 'family' must be \"probit\" or \"logistic\".")
  }
  if (any(is.na(Bbasis))) {
    stop("The 'Bbasis' matrix contains missing values. B-spline basis must be fully computed.")
  }
  cols2check <- intersect(setdiff(all.vars(formula), x_name), colnames(data))
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
  Xmat_tmp <- stats::model.matrix(formula, mf)
  Xmat_assign <- attr(Xmat_tmp, "assign")[-1]
  Xmat <- Xmat_tmp[, -1, drop = FALSE]  # Remove intercept
  # Split by selection indicator S
  yvec_s1 <- yvec[n1_idx]; yvec_s0 <- yvec[n0_idx]
  Xmat_s1 <- Xmat[n1_idx, , drop = FALSE]; Xmat_s0 <- Xmat[n0_idx, , drop = FALSE]
  Bbasis_s1 <- Bbasis[n1_idx, , drop = FALSE]; Bbasis_s0 <- Bbasis[n0_idx, , drop = FALSE]

  # Identify model matrix columns associated with x_name
  terms_mat <- attr(terms(formula), "factors")
  x_rowidx <- which(sapply(rownames(terms_mat), function(r) x_name %in% all.vars(parse(text = r))))
  x_colidx <- which(terms_mat[x_rowidx, , drop = FALSE] > 0)
  xonly_colidx <- x_colidx[colSums(terms_mat[, x_colidx, drop = FALSE]) == 1]
  xonly_colname <- colnames(Xmat)[Xmat_assign %in% xonly_colidx]
  # Interaction term with X?
  xinterz_colname <- character(0)
  xinterz_z_colname <- character(0)
  xinterz_colidx <- x_colidx[colSums(terms_mat[, x_colidx, drop = FALSE]) > 1]
  if (length(xinterz_colidx) > 0) {
    xinterz_colname <- colnames(Xmat)[Xmat_assign %in% xinterz_colidx]
    xinterz_z_rowidx <- setdiff(which(terms_mat[, xinterz_colidx] > 0), x_rowidx)
    if (length(xinterz_z_rowidx) > 1) {
      stop("Higher-order interactions (e.g., X:Z1:Z2) are not supported. Only two-way interactions allowed.")
    }
    xinterz_z_colname <- colnames(Xmat)[Xmat_assign == xinterz_z_rowidx]
  }

  # Support points for X: unique raw values observed in Phase 2
  x_raw_s1 <- data[[x_name]][n1_idx]
  x_raw_unique <- sort(unique(x_raw_s1))
  # Map each Phase 2 subject to their support point index
  s1_raw_match <- match(x_raw_s1, x_raw_unique)
  # Build x_support as model matrix columns for each unique raw X value
  # Create a template data frame with one row per support point
  template_df <- data[1:length(x_raw_unique), , drop = FALSE]
  template_df[[x_name]] <- x_raw_unique
  mf_template <- stats::model.frame(formula, template_df, na.action = na.pass)
  Xmat_template <- stats::model.matrix(formula, mf_template)[, -1, drop = FALSE]
  x_support <- Xmat_template[, xonly_colname, drop = FALSE]

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
    fit_init <- tryCatch(
      MASS::polr(formula, data, method = family),
      error = function(e) NULL
    )
    if (!is.null(fit_init)) {
      theta_curr <- as.vector(c(fit_init$coefficients, fit_init$zeta))
    } else {
      theta_curr <- c(rep(1, ncol(Xmat_tmp) - 1), 1:(length(unique(yvec)) - 1))
    }
  }

  # p_vl: normalized sieve coefficients
  # Use support point index matching (exact, not string-based)
  d_size <- nrow(x_support)
  indicator_mat <- matrix(0, nrow = length(n1_idx), ncol = d_size)
  indicator_mat[cbind(seq_along(s1_raw_match), s1_raw_match)] <- 1
  p_vl_num_init <- t(indicator_mat) %*% Bbasis_s1  # (d, s_n)
  p_vl_curr <- sweep(p_vl_num_init, 2, pmax(1e-16, colSums(p_vl_num_init)), FUN = "/")
  if (!is.null(p_vl_init)) {
    if (!is.matrix(p_vl_init) || !identical(dim(p_vl_init), dim(p_vl_curr))) {
      warning(sprintf("p_vl_init must be a %s matrix; got %s. Ignoring warm start.",
                      paste(dim(p_vl_curr), collapse = "x"),
                      if (is.matrix(p_vl_init)) paste(dim(p_vl_init), collapse = "x")
                      else paste0(class(p_vl_init)[1], " of length ", length(p_vl_init))))
    } else {
      p_vl_curr <- p_vl_init
    }
  }

  # -- 3. EM Loop ----
  converged <- FALSE
  for (iter in 1:max_iter) {
    theta_old <- theta_curr
    p_vl_old <- p_vl_curr

    ## -- E-STEP ----
    # For S=0, compute w_iv = E[I_iv | Y_i, Z_i, theta(m), p_vl(m)]
    prob_y_given_xv_s0 <- compute_py_given_xv_s0(theta_curr, yvec_s0, Xmat_s0, x_support, xonly_colname, family,
                                                 x_inter_colname = xinterz_colname, z_inter_colname = xinterz_z_colname)  # (n0, d)
    w_iv <- compute_w_iv_s0(prob_y_given_xv_s0, p_vl_curr, Bbasis_s0)  # (n0, d)

    # -- M-STEP ----
    # Update theta(m+1)
    theta_curr <- update_theta_smle(theta_curr, w_iv, yvec_s1, yvec_s0, Xmat_s1, Xmat_s0, x_support, xonly_colname, family,
                                    x_inter_colname = xinterz_colname, z_inter_colname = xinterz_z_colname)
    # Update p_vl(m+1)
    p_vl_curr <- update_p_vl_cpp(prob_y_given_xv_s0, Bbasis_s0, p_vl_curr, p_vl_num_init, max_iter, tol)

    # Check convergence on both theta and p_vl
    if (max(abs(theta_curr - theta_old)) < tol && max(abs(p_vl_curr - p_vl_old)) < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    message(sprintf("EM algorithm reached max_iter without converging to tolerance.\n
                    Current max(abs(theta_curr - theta_old)) = %.6f\n", max(abs(theta_curr - theta_old))))
  }
  if (verbose && converged) {
    cat(sprintf("EM converged at iter = %d (tol = %g)\n", iter, tol))
  }

  names(theta_curr)[1:ncol(Xmat)] <- colnames(Xmat)
  names(theta_curr)[(ncol(Xmat) + 1):length(theta_curr)] <- paste0("alpha", 1:(length(theta_curr) - ncol(Xmat)))

  # -- 4. Standard Error Estimation (Envelope Theorem) ----
  se <- NULL
  vcov <- NULL
  se_max_iterations <- NULL
  se_diagnostics <- NULL
  if (se_calc) {
    se_results <- estimate_se_smle_envelope(
      theta = theta_curr, p_vl = p_vl_curr, p_vl_s1 = p_vl_num_init,
      yvec_s1, yvec_s0, Xmat_s1, Xmat_s0, Bbasis_s0, x_support, xonly_colname, family, se_max_iter, se_tol,
      verbose = verbose, x_inter_colname = xinterz_colname, z_inter_colname = xinterz_z_colname,
      jacobian_method = jacobian_method, jacobian_method_args = jacobian_method_args
    )
    se <- se_results$se
    vcov <- se_results$vcov
    se_max_iterations <- se_results$se_max_iterations
    se_diagnostics <- se_results$diagnostics
    names(se) <- rownames(vcov) <- colnames(vcov) <- names(theta_curr)
  }

  return(list(est = theta_curr, se = se, vcov = vcov, p_vl = p_vl_curr,
              iterations = iter, se_max_iterations = se_max_iterations,
              se_diagnostics = se_diagnostics))
}
