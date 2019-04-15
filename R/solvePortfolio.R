
solvePortfolio <- function(p, w0, m1, M2, M3, M4, lb, ub, lin_eq, lin_eqC,
                           nlin_eq, lin_ieq, lin_ieqC, nlin_ieq, options, param) {

  ### input
  # p         : dimension of the portfolio
  # w0        : name of portfolio to optimize
  # m1        : vector with expected returns
  # M2        : covariance matrix
  # M3        : coskewness matrix
  # M4        : cokurtosis matrix
  # lb        : vector with lower bounds for the weights
  # ub        : vector with upper bounds for the weights
  # lin_eq    : equality constraints: eq w = eqC (should be matrix!), see details
  # lin_eqC   : equality constraints: eq w = eqC, see details
  # nlin_eq   : function with non-linear equality constraint (returns objective value and jacobian)
  # lin_ieq   : inequality constraints: ieq w leq ieqC (should be matrix!)
  # lin_ieqC  : inequality constraints: ieq w leq ieqC
  # nlin_ieq  : function with non-linear inequality constraints (returns objective value and jacobian)
  # options   : optimization options
  # param     : additional parameters in a named list for the objective function (such as gamma for EU)
  #
  ### output
  # vector with optimal weights

  ### optimization options
  if (!("maxeval" %in% names(options))) options$maxeval = 10000
  if (!("check_derivatives" %in% names(options))) options$check_derivatives = FALSE
  if (!("print_level" %in% names(options))) options$print_level = 0
  options$algorithm <- "NLOPT_LD_SLSQP"


  ### set up constraints
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

    list(constraints = cts, jacobian = jac)
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

    list(constraints = cts, jacobian = jac)
  }


  ### Select objective function
  fn <- get_href(w0, m1, M2, M3, M4, param)


  ### optimize portfolio
  w00 <- rep(1 / p, p)
  sol <- nloptr::nloptr(x0 = w00, eval_f = fn, lb = lb, ub = ub,
                        eval_g_eq = g_eq, eval_g_ineq = g_ineq, opts = options)

  sol$solution
}

