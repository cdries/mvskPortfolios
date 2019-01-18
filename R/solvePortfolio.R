
solvePortfolio <- function(p, w0, m1, M2, M3, M4, lb, ub, lin_eq, lin_eqC,
                           nlin_eq, lin_ieq, lin_ieqC, nlin_ieq, options) {

  ### optimization options
  if (!("maxeval" %in% names(options))) options$maxeval = 10000
  if (!("check_derivatives" %in% names(options))) options$check_derivatives = FALSE
  if (!("print_level" %in% names(options))) options$print_level = 0
  options$algorithm <- "NLOPT_LD_SLSQP"


  ### set up constraints
  # box constraints
  if (is.null(lb)) lb <- rep(0, p)
  if (is.null(ub)) ub <- rep(1, p)

  # equality constraints
  g_eq <- function(w) {

    cts <- jac <- NULL

    # linear constraints
    if (!is.null(lin_eq)) {
      cts <- lin_eq %*% w - lin_eqC
      jac <- lin_eq
    }

    # non-linear constraints
    if (!is.null(nlin_eq)) {
      nlin_eq_res <- nlin_eq(w)
      cts <- c(cts, nlin_eq_res$constraints)
      jac <- rbind(jac, nlin_eq_res$jacobian)
    }

    return (list("constraints" = cts, "jacobian" = jac))
  }

  # inequality constraints
  g_ineq <- function(w) {

    cts <- jac <- NULL

    # linear constraints
    if (!is.null(lin_ieq)) {
      cts <- lin_ieq %*% w - lin_ieqC
      jac <- lin_ieq
    }

    # non-linear constraints
    if (!is.null(nlin_ieq)) {
      nlin_ieq_res <- nlin_ieq(w)
      cts <- c(cts, nlin_ieq_res$constraints)
      jac <- rbind(jac, nlin_ieq_res$jacobian)
    }

    return (list("constraints" = cts, "jacobian" = jac))
  }


  ### Select objective function
  if (w0 == "maxDiv") {
    fn <- function(w) fDR(w, M2)
  } else if (w0 == "ERC") {
    # TODO
  } else {
    stop("Choose a valid portfolio objective")
  }


  ### optimize portfolio
  w00 <- rep(1 / p, p)
  sol <- nloptr::nloptr(x0 = w00, eval_f = fn, lb = lb, ub = ub,
                        eval_g_eq = g_eq, eval_g_ineq = g_ineq, opts = options)

  return (list("wopt" = sol$solution, "val" = sol$objective, "name" = w0))
}

