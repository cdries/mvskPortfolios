
rm(list = ls())

# load data ---------------------------------------------------------------

data(indices)
indices <- indices[-c(1, 2),]           # first two weeks contain NA
x <- as.matrix(indices)


# summary statistics ------------------------------------------------------

summ <- matrix(NA, nrow = ncol(x), ncol = 4)
rownames(summ) <- colnames(indices)
colnames(summ) <- c("ann. geometric mean (%)", "ann. standard deviation", "standardized skewness", "excess kurtosis")
summ[, 1] <- apply(x, 2, function(a) (prod(1 + a))^(52 / length(a)) - 1) * 100
summ[, 2] <- apply(x, 2, sd) * sqrt(52)
summ[, 3] <- apply(x, 2, function(a) mean((a - mean(a))^3)) / apply(x, 2, sd)^3
summ[, 4] <- apply(x, 2, function(a) mean((a - mean(a))^4)) / apply(x, 2, sd)^4 - 3

xtable::xtable(summ, digits = 2)
