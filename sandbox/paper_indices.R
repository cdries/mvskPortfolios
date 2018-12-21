rm(list = ls())

# this file implements a backtest on the indices data

library(mvskPortfolios)
library(PerformanceAnalytics)


# load data ---------------------------------------------------------------

data(indices)
indices <- indices[-c(1, 2),]           # first two weeks contain NA
x <- as.matrix(indices)


# settings and initialisation ---------------------------------------------

p <- ncol(x)
lb <- rep(0, p)
ub <- rep(0.2, p)

n <- 10 * 52
kappa <- c(0, 0.01, 0.05, 0.1, 1)
w_array <- array(NA, dim = c(nrow(x), p, 1 + length(kappa)))

w0 <- "maxDiv"


# run through time --------------------------------------------------------

for (tt in n:(nrow(x) - 1)) {

  if (tt %% 10 == 0) cat("---------", tt, "out of", nrow(x) - 1, "---------\n")

  # subset data
  z <- x[(tt - n + 1):tt,]

  # estimate moments
  M2 <- cov(z)
  M3 <- M3.MM(z, as.mat = FALSE)
  M4 <- M4.MM(z, as.mat = FALSE)

  # equally-weighted portfolio
  w_array[tt + 1,, 1] <- rep(1 / p, p)

  # efficient portfolios
  res <- mvskPortfolio(M2 = M2, M3 = M3, M4 = M4, w0 = w0, kappa = kappa, lb = lb, ub = ub)
  w_array[tt + 1,, -1] <- t(res$summ$w)
}


# compute return series ---------------------------------------------------

ind <- !is.na(w_array[, 1, 1])
z <- x[ind, ]
w_oos <- w_array[ind,,]
dates <- index(indices)[ind]

ret_mat <- apply(w_oos, 3, function(a) rowSums(a * z))


# examine results ---------------------------------------------------------

hist_geoommu <- (apply(ret_mat, 2, function(x) prod(1 + x)^(52 / length(x))) - 1)
hist_sd <- apply(ret_mat, 2, stats::sd)
hist_m3 <- apply(ret_mat, 2, function(x) mean((x - mean(x))^3))
hist_m4 <- apply(ret_mat, 2, function(x) mean((x - mean(x))^4))
hist_skew <- apply(ret_mat, 2, function(x) mean((x - mean(x))^3)) / hist_sd^3
hist_exkurt <- apply(ret_mat, 2, function(x) mean((x - mean(x))^4)) / hist_sd^4 - 3

tb <- rbind(hist_geoommu * 100, hist_sd * sqrt(52), hist_m3 * 1e6, hist_m4 * 1e7,
            hist_skew, hist_exkurt)
colnames(tb) <- c("EW", paste0("VSK - kappa ", kappa))
rownames(tb) <- c("geometric mean (%)", "standard deviation (ann.)", "m3 (1e-6)", "m4 (1e-7)",
                  "standardized skewness", "excess kurtosis")
print(tb)

cols <- c("black", "red", "darkblue", "blue", "cyan", "darkgreen")
strats <- colnames(tb)

# cumulative performance
plot(dates, cumprod(1 + ret_mat[, 1]), type = 'l', lwd = 1.2, col = "black",
     main = "cumulative performance", xlab = "time", ylab = "cumulative return")
for (ii in 2:ncol(ret_mat)) {
  lines(dates, cumprod(1 + ret_mat[, ii]), lwd = 1.2, col = cols[ii])
}
legend("topleft", strats, col = cols, lwd = rep(2, length(cols)))

# drawdown curves
drd <- apply(ret_mat, 2, function(x) (cumprod(1 + x) - cummax(cumprod(1 + x))) /
               cummax(cumprod(1 + x)) * 100)
plot(dates, drd[, 1], type = 'l', lwd = 1.2, col = "black", xlab = "time", ylab = "drawdown",
     main = "drawdown curves (% from previous high)")
for (ii in 2:ncol(drd)) {
  lines(dates, drd[, ii], lwd = 1.2, col = cols[ii])
}
legend("bottomright", strats, col = cols, lwd = rep(2, length(cols)))

