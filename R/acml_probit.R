#' ACML under probit
#'
#' @param formula A model formula specifying the probit model.
#' @param data A data frame containing the variables in the formula.
#' @param pi_values Numeric vector of sampling probabilities for ODS.
#'
#' @return A list with parameter estimates.
#' @export
acml_probit <- function(formula, data, pi_values) {
  # -- 0. Validate ----
  na_counts <- colSums(is.na(data[, all.vars(formula), drop = FALSE]))
  if (any(na_counts > 0)) {
    missing_cols <- names(na_counts[na_counts > 0])
    stop(paste("Missing values detected in the following columns:",
               paste(missing_cols, collapse = ", ")))
  }

  # -- 1. Construct Data Matrices ----
  mf <- model.frame(formula, data)
  y_vec <- model.response(mf)
  X_mat <- model.matrix(formula, mf)[, -1, drop = FALSE]  # remove intercept (column 1)

  # -- 2. Setup Dimensions ----
  K <- length(pi_values)
  p <- ncol(X_mat)

  # -- 3. Initialize Parameters ----
  # Initialize p betas at 0 and K-1 alphas spread from -1 to 1
  start_params <- c(numeric(p), seq(-1, 1, length.out = K - 1))
  beta_names <- colnames(X_mat)
  if (is.null(beta_names)) beta_names <- paste0("beta", 1:p)
  alpha_names <- paste0("(Intercept:", 1:(K - 1), ")")
  names(start_params) <- c(beta_names, alpha_names)

  # -- 4. Optimization ----
  fit <- optim(
    par = start_params,
    fn = acml_probit_nll,
    gr = acml_probit_grad,
    X_mat = X_mat,
    y_vec = y_vec,
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
    Term = names(start_params),
    Estimate = fit$par,
    Std.Err = se,
    Z = fit$par / se,
    P = 2 * (1 - pnorm(abs(fit$par / se)))
  )

  return(list(summary = res, fit = fit, vcov = vcov_mat))
}
