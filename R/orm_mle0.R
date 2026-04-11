#' MLE0 Estimator under two-phase sampling designs for ordinal outcomes
#'
#' A special case of SMLE that treats all covariates as missing for Phase 1-only subjects
#' and models the covariate distribution as a discrete distribution over the unique
#' covariate combinations observed in Phase 2. Does not require a B-spline basis.
#'
#' @param formula Formula for the ordinal outcome Y.
#' @param data Data frame. NA allowed only in \code{x_name} for S=0 subjects.
#' @param x_name Name of the expensive covariate X (used to identify Phase 2 subjects via missingness).
#' @param family Character value specifying the distribution family, which is one of the following: "logit", "probit"
#' @param theta_init Optional initial vector for theta (beta, cutpoints).
#' @param se_calc Logical; whether to compute standard errors. Defaults to TRUE.
#' @param se_method Character value specifying the method for SE estimation: "forward" uses second-order forward-difference Hessian (default) and "numDeriv" uses Richardson extrapolation via \code{numDeriv::hessian()}.
#' @param verbose Logical; if TRUE, prints convergence info after estimation. Defaults to FALSE.
#' @param max_iter Maximum number of EM iterations. Defaults to 500.
#' @param tol Convergence tolerance for the optimizer. Defaults to 1e-6.
#' @param se_max_iter Maximum number of EM iterations for finding optimal q_j in SE estimation. Defaults to 1000.
#' @param se_tol Convergence tolerance for the optimizer in SE estimation. Defaults to 1e-8.
#'
#' @return A list containing:
#'   \item{est}{Estimated parameter vector (beta, cutpoints).}
#'   \item{se}{Standard errors (if \code{se_calc = TRUE}).}
#'   \item{vcov}{Variance-covariance matrix (if \code{se_calc = TRUE}).}
#'   \item{q_j}{Estimated discrete distribution over support points.}
#'   \item{iterations}{Number of EM iterations.}
#'   \item{se_max_iterations}{Maximum inner EM iterations during SE estimation.}
#' @export
orm_mle0 <- function(formula, data, x_name, family = "probit",
                     theta_init = NULL, se_calc = TRUE, se_method = "forward",
                     verbose = TRUE, max_iter = 500, tol = 1e-6,
                     se_max_iter = 1000, se_tol = 1e-8) {
  # -- 0. Validate ----
  se_method <- match.arg(se_method, c("forward", "numDeriv"))
  if (!(family %in% c("probit", "logistic"))) {
    stop("The 'family' must be \"probit\" or \"logistic\".")
  }
  cols2check <- intersect(setdiff(all.vars(formula), x_name), colnames(data))
  na_counts <- colSums(is.na(data[, cols2check, drop = FALSE]))
  if (any(na_counts > 0)) {
    missing_cols <- names(na_counts[na_counts > 0])
    stop(paste("Missing values detected in the following columns:",
               paste(missing_cols, collapse = ", "),
               "- Only the column specified in 'x_name' can have missingness."))
  }

  # Identify selection indicator (S=1 if x_name is NOT NA, S=0 if x_name is NA)
  s_indicator <- as.numeric(!is.na(data[[x_name]]))
  n1_idx <- which(s_indicator == 1); n0_idx <- which(s_indicator == 0)
  n <- nrow(data)

  # -- 1. Construct Data Matrices ----
  mf <- stats::model.frame(formula, data, na.action = na.pass)
  yvec <- as.numeric(stats::model.response(mf))
  Xmat_tmp <- stats::model.matrix(formula, mf)
  Xmat <- Xmat_tmp[, -1, drop = FALSE]  # Remove intercept

  yvec_s1 <- yvec[n1_idx]; yvec_s0 <- yvec[n0_idx]
  Xmat_s1 <- Xmat[n1_idx, , drop = FALSE]

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

  # Support points: unique rows of the full model matrix from Phase 2
  row_key <- apply(round(Xmat_s1, 12), 1, paste, collapse = "\t")
  unique_key <- unique(row_key)
  support_mat <- Xmat_s1[match(unique_key, row_key), , drop = FALSE]
  s1_support_idx <- match(row_key, unique_key)
  m <- nrow(support_mat)

  # q_j: frequency of each support point in Phase 2
  q_j_curr <- as.vector(tabulate(s1_support_idx, nbins = m)) / length(n1_idx)

  # -- 3. EM Loop ----
  converged <- FALSE
  for (iter in 1:max_iter) {
    theta_old <- theta_curr
    q_j_old <- q_j_curr

    ## -- E-STEP: compute psi_ij for S=0 ----
    prob_y_s0 <- compute_py_given_support(theta_curr, yvec_s0, support_mat, family)
    psi_ij <- compute_psi_mle0(prob_y_s0, q_j_curr)

    ## -- M-STEP ----
    # Update theta
    theta_curr <- update_theta_mle0(theta_curr, psi_ij, yvec_s1, yvec_s0, Xmat_s1, support_mat, family)
    # Update q_j
    q_j_curr <- update_q_mle0(psi_ij, s1_support_idx, n)

    # Check convergence
    if (max(abs(theta_curr - theta_old)) < tol) {
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

  # -- 4. Standard Error Estimation ----
  se <- NULL
  vcov <- NULL
  se_max_iterations <- NULL
  if (se_calc) {
    se_results <- estimate_se_mle0(
      theta = theta_curr, q_j = q_j_curr, yvec_s1, yvec_s0, Xmat_s1,
      support_mat, s1_support_idx, family, se_max_iter, se_tol,
      method = se_method, verbose = verbose
    )
    se <- se_results$se
    vcov <- se_results$vcov
    se_max_iterations <- se_results$se_max_iterations
    names(se) <- rownames(vcov) <- colnames(vcov) <- names(theta_curr)
  }

  return(list(est = theta_curr, se = se, vcov = vcov, q_j = q_j_curr, iterations = iter, se_max_iterations = se_max_iterations))
}
