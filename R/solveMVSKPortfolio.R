
solveMVSKPortfolio <- function(p, w0, g, m1, M2, M3, M4, indmom, lb, ub, lin_eq, lin_eqC,
                               nlin_eq, lin_ieq, lin_ieqC, nlin_ieq, options,
                               href, kappa, relative, param, mompref) {

  ### input
  # p         : dimension of the portfolio
  # w0        : numeric vector of length w0 with reference portfolio weights
  # g         : either function g(delta) or vector for direction of improvement
  # m1        : vector with expected returns
  # M2        : covariance matrix
  # M3        : coskewness matrix
  # M4        : cokurtosis matrix
  # indmom    : boolean vector of length 4 indicating which moments to use
  # lb        : vector with lower bounds for the weights
  # ub        : vector with upper bounds for the weights
  # lin_eq    : equality constraints: eq w = eqC (should be matrix!), see details
  # lin_eqC   : equality constraints: eq w = eqC, see details
  # nlin_eq   : function with non-linear equality constraint (returns objective value and jacobian)
  # lin_ieq   : inequality constraints: ieq w leq ieqC (should be matrix!)
  # lin_ieqC  : inequality constraints: ieq w leq ieqC
  # nlin_ieq  : function with non-linear inequality constraints (returns objective value and jacobian)
  # options   : optimization options
  # href      : either function or one of the preset names; href(w) leq href(w0) + kappa
  # kappa     : margin on href
  # relative  : determines if margin is relative or absolute;
  #             in case of relative: href(w) leq (1 + sign(href(w0) * kappa) href(w0)
  # param     : list with additional parameters for the href function
  # mompref   : moment preferences (+1 for higher, -1 for lower)
  #
  ### output
  # objective : minus the diversification ratio
  # gradient  : the gradient of the objective with respect to w

  ### optimization options
  if (!("maxeval" %in% names(options))) options$maxeval = 10000
  if (!("check_derivatives" %in% names(options))) options$check_derivatives = FALSE
  if (!("print_level" %in% names(options))) options$print_level = 0
  options$algorithm <- "NLOPT_LD_SLSQP"


  ### set up constraints
  # equality constraints
  g_eq <- function(x) {

    w <- x[-1]
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

    if (!is.null(cts)) jac <- cbind(0, jac)

    list(constraints = cts, jacobian = jac)
  }

  # inequality constraints
  if (is.null(href)) {
    objw0 <- NULL
    fn <- function(w) NULL
  } else {
    if (is.function(href)) {
      fn <- href
    } else {
      if (href == "TEvol" && !("wref" %in% names(param))) param <- list(wref = w0)
      fn <- get_href(href, m1, M2, M3, M4, param)
    }
    fn0 <- fn(w0)$objective
    if (relative) objw0 <- (1 + sign(fn0) * kappa) * fn0 else objw0 <- (fn0 + kappa)
  }
  mw0 <- getmom(indmom, w0, m1, M2, M3, M4)
  if (is.null(mompref)) sgm <- c(-1, 1, -1, 1)[indmom] else sgm <- -mompref

  g_ineq <- function(x) {

    delta <- x[1]
    w <- x[-1]
    objRtemp <- fn(w)

    # MVSK constraints
    if (is.numeric(g)) {
      gf <- delta * g
      gfgr <- g
    } else {
      gtmp <- g(delta)
      gf <- gtmp$objective
      gfgr <- gtmp$gradient
    }
    mw <- getmom(indmom, w, m1, M2, M3, M4)
    obj <- sgm * (mw - mw0) + gf
    obj <- c(obj, objRtemp$objective - objw0)

    momgr <- getmomgr(indmom, w, m1, M2, M3, M4)
    momgr <- rbind(cbind(gfgr, momgr * sgm), c(0, objRtemp$gradient))

    # linear and non-linear constraints
    cts <- jac <- NULL
    if (!is.null(lin_ieq)) {
      cts <- lin_ieq %*% w - lin_ieqC
      jac <- lin_ieq
    }
    if (!is.null(nlin_ieq)) {
      nlin_ieq_res <- nlin_ieq(w)
      cts <- c(cts, nlin_ieq_res$constraints)
      jac <- rbind(jac, nlin_ieq_res$jacobian)
    }
    if (!is.null(jac)) jac <- cbind(0, jac)

    # combine constraints
    cts <- c(obj, cts)
    jac <- rbind(momgr, jac)[1:length(cts),]

    list(constraints = cts, jacobian = jac)
  }


  ### objective function
  fobj <- function(x) {
    obj <- -x[1]
    gr <- rep(0, length(x))
    gr[1] <- -1
    list(objective = obj, gradient = gr)
  }


  ### optimize portfolio
  sol <- nloptr::nloptr(x0 = c(0, w0), eval_f = fobj, lb = c(-1, lb), ub = c(1, ub),
                        eval_g_eq = g_eq, eval_g_ineq = g_ineq, opts = options)
  wopt <- sol$solution[-1]
  delta <- sol$solution[1]
  moms <- getmom(indmom, wopt, m1, M2, M3, M4)
  constr <- sol$eval_g_ineq(sol$solution)$constraints

  list(w = wopt, delta = delta, moms = moms, constr = constr)
}
