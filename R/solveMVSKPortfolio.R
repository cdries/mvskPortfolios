
solveMVSKPortfolio <- function(p, initport, kappa, g, m1, M2, M3, M4, indmom, lb, ub, lin_eq,
                               lin_eqC, nlin_eq, lin_ieq, lin_ieqC, nlin_ieq, options) {

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

    return (list("constraints" = cts, "jacobian" = jac))
  }

  # inequality constraints
  objw0 <- -(1 - kappa) * initport$val
  mw0 <- getmom(indmom, initport$wopt, m1, M2, M3, M4)
  sgm <- c(-1, 1, -1, 1)[indmom]

  if (initport$name == "maxDiv") {
    fn <- function(w) fDR(w, M2)
  } else if (initport$name == "ERC") {
    # TODO
  } else {
    stop("Choose a valid portfolio objective")
  }

  g_ineq <- function(x) {

    delta <- x[1]
    w <- x[-1]
    objRtemp <- fn(w)

    # MVSK constraints
    mw <- getmom(indmom, w, m1, M2, M3, M4)
    obj <- sgm * (mw - mw0) + delta * g
    obj <- c(obj, objw0 + objRtemp$objective)

    momgr <- getmomgr(indmom, w, m1, M2, M3, M4)
    momgr <- rbind(cbind(g, momgr * sgm), c(0, objRtemp$gradient))

    # linear and non-linear constraints
    cts <- jac <- NULL
    if (!is.null(lin_ieq)) {
      cts <- lin_ieq %*% w - lin_ieqC
      jac <- lin_ieq
    }
    if (!is.null(nlin_ieq)) {
      nlin_ieq_res <- nlin_eq(w)
      cts <- c(cts, nlin_ieq_res$constraints)
      jac <- rbind(jac, nlin_ieq_res$jacobian)
    }
    if (!is.null(jac)) jac <- cbind(0, jac)

    # combine constraints
    cts <- c(obj, cts)
    jac <- rbind(momgr, jac)

    return (list("constraints" = cts, "jacobian" = jac))
  }

  # objective function
  fobj <- function(x) {
    obj <- -x[1]
    gr <- rep(0, length(x))
    gr[1] <- -1
    return (list("objective" = obj, "gradient" = gr))
  }


  ### optimize portfolio
  sol <- nloptr::nloptr(x0 = c(0, initport$wopt), eval_f = fobj, lb = c(-1, lb), ub = c(1, ub),
                        eval_g_eq = g_eq, eval_g_ineq = g_ineq, opts = options)
  wopt <- sol$solution[-1]
  delta <- sol$solution[1]
  moms <- getmom(indmom, wopt, m1, M2, M3, M4)

  return (list("wopt" = wopt, "delta" = delta, "moms" = moms, "name" = initport$name))
}
