
solveMVSKPortfolio <- function(p, w0, kappa, g, m1, M2, M3, M4, indmom, lb, ub, lin_eq, lin_eqC, nlin_eq,
                               lin_ieq, lin_ieqC, nlin_ieq, options, relative, fnPerf = NULL, mpref = NULL) {

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
  if (is.null(kappa)) {
    objw0 <- NULL
    fn <- function(w) NULL
  } else {
    if (relative) objw0 <- -(1 - kappa) * fnPerf$val else objw0 <- -(fnPerf$val + kappa)
    if (fnPerf$name == "maxDiv") {
      fn <- function(w) fDR(w, M2)
    } else if (fnPerf$name == "ERC") {
      fn <- function(w) fERC(w, M2)
    } else if (fnPerf$name == "HI") {
      fn <- function(w) fHI(w)
    } else {
      stop("Choose a valid portfolio objective")
    }
  }
  mw0 <- getmom(indmom, w0, m1, M2, M3, M4)
  if (is.null(mpref)) sgm <- c(-1, 1, -1, 1)[indmom] else sgm <- -mpref

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
    obj <- c(obj, objw0 + objRtemp$objective)

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
  sol <- nloptr::nloptr(x0 = c(0, w0), eval_f = fobj, lb = c(-1, lb), ub = c(1, ub),
                        eval_g_eq = g_eq, eval_g_ineq = g_ineq, opts = options)
  wopt <- sol$solution[-1]
  delta <- sol$solution[1]
  moms <- getmom(indmom, wopt, m1, M2, M3, M4)
  constr <- sol$eval_g_ineq(sol$solution)$constraints

  return (list("wopt" = wopt, "delta" = delta, "moms" = moms, "constr" = constr))
}
