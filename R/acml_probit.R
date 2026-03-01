#' ACML under probit
#'
#' @param formula A model formula specifying the probit model.
#' @param data A data frame containing the variables in the formula.
#' @param pi_values Numeric vector of sampling probabilities for ODS.
#' @param theta_init Optional numeric vector of initial values for theta.
#'
#' @return A list with parameter estimates.
#' @export
acml_probit <- function(formula, data, pi_values, theta_init = NULL) {
  # -- 0. Validate ----
  na_counts <- colSums(is.na(data[, all.vars(formula), drop = FALSE]))
  if (any(na_counts > 0)) {
    missing_cols <- names(na_counts[na_counts > 0])
    stop(paste("Missing values detected in the following columns:",
               paste(missing_cols, collapse = ", ")))
  }

  # -- 1. Construct Data Matrices ----
  mf <- model.frame(formula, data)
  yvec <- model.response(mf)
  Xmat <- model.matrix(formula, mf)[, -1, drop = FALSE]  # remove intercept (column 1)

  # -- 2. Setup Dimensions ----
  K <- length(pi_values)
  p <- ncol(Xmat)

  # -- 3. Initialize Parameters ----
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
    fit_init <- MASS::polr(formula, data, method = "probit")
    theta_init <- as.vector(c(fit_init$coefficients, fit_init$zeta))
  }
  beta_names <- colnames(Xmat)
  if (is.null(beta_names)) beta_names <- paste0("beta", 1:p)
  alpha_names <- paste0("(Intercept:", 1:(K - 1), ")")
  names(theta_init) <- c(beta_names, alpha_names)

  # -- 4. Optimization ----
  fit <- optim(
    par = theta_init,
    fn = acml_probit_nll,
    gr = acml_probit_grad,
    Xmat = Xmat,
    yvec = yvec,
    pi_perY = pi_values,
    method = "BFGS",
    hessian = TRUE
  )

  # -- 5. Process Results ----
  if(fit$convergence != 0) {
    warning("Optimization did not reach convergence (code ", fit$convergence, ").")
  }

  # Calculate standard errors from the inverse Hessian
  vcov_mat <- solve(fit$hessian)
  se <- sqrt(diag(vcov_mat))

  res <- data.frame(
    Term = names(theta_init),
    Estimate = fit$par,
    Std.Err = se,
    Z = fit$par / se,
    P = 2 * (1 - pnorm(abs(fit$par / se)))
  )

  return(list(summary = res, fit = fit, vcov = vcov_mat))
}
