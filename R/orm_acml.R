#' Ascertainment Corrected Maximum Likelihood (ACML)under two-phase sampling designs for ordinal outcomes
#'
#' @param formula Formula for the ordinal outcome Y.
#' @param data Data frame containing all variables in \code{formula}. No NA allowed.
#' @param pi_values Numeric vector of sampling probabilities for ODS.
#' @param family Character value specifying the distribution family, which is one of the following: "logit", "probit"
#' @param theta_init Optional initial vector for theta (beta, cutpoints).
#'
#' @return A list with containing estimates.
#' @export
orm_acml <- function(formula, data, pi_values, family = "probit", theta_init = NULL) {
  # -- 0. Validate ----
  if (!(family %in% c("probit", "logistic"))) {
    stop("The 'family' must be \"probit\" or \"logistic\".")
  }
  na_counts <- colSums(is.na(data[, all.vars(formula), drop = FALSE]))
  if (any(na_counts > 0)) {
    missing_cols <- names(na_counts[na_counts > 0])
    stop(paste("Missing values detected in the following columns:",
               paste(missing_cols, collapse = ", ")))
  }

  # -- 1. Construct Data Matrices ----
  mf <- stats::model.frame(formula, data)
  yvec <- stats::model.response(mf)
  Xmat <- stats::model.matrix(formula, mf)[, -1, drop = FALSE]  # Remove intercept
  K <- length(pi_values)
  p <- ncol(Xmat)

  # -- 2. Initialize Parameters ----
  use_defaults <- TRUE
  if (!is.null(theta_init)) {
    theta_len <- p + K - 1
    if (length(theta_init) == theta_len) {
      use_defaults <- FALSE
    } else {
      warning(paste("theta_init must have length", theta_len))
    }
  }
  if (use_defaults) {
    fit_init <- MASS::polr(formula, data, method = family)
    theta_init <- as.vector(c(fit_init$coefficients, fit_init$zeta))
  }

  # -- 3. Optimization ----
  fit <- stats::optim(
    par = theta_init, fn = acml_probit_nll, gr = acml_probit_grad,
    yvec = yvec, Xmat = Xmat, pi_perY = pi_values, family,
    method = "BFGS", hessian = TRUE
  )

  if(fit$convergence != 0) {
    warning("Optimization did not reach convergence (code ", fit$convergence, ").")
  }

  # -- 4. Standard Error Estimation ----
  vcov_mat <- solve(fit$hessian)
  se <- sqrt(diag(vcov_mat))

  return(list(est = fit$par, se = se, vcov = vcov_mat, convergence = fit$convergence))
}
