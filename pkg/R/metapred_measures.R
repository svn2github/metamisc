### Performance / error functions
# Error function: Mean Squared Error
mse <- brier <- function(p, y, ...) mean((p - y)^2)


# Error function: Variance of prediction error
vare <- function(p, y, ...) var(p - y)

absmean <- function(perf.measures, ...) {
  pm <- unlist(perf.measures)
  abs(mean(pm))
}

# Measure 1: Coefficient of variation of prediction error.
coefVarPred <- function(p, y, data = NULL, abs = TRUE, ...)
  coefVar(x = p - y, abs = abs, ...) 

### Heterogeneity, generalizability, pooled performance functions

# Measure 1: Coefficient of variation (=scaled sd)
# In general sense, abs needs not be TRUE, but for metapred it should,
# such that higher values are worse performance.
coefVar <- function(x, abs = TRUE, ...) {
  cv <- sd(x)/mean(x)
  if (isTRUE(abs)) abs(cv) else cv
}

coefVarMean <- function(x, abs = TRUE, ...) 
  coefVar(x, abs = abs) + if (abs) abs(mean(x)) else mean(x)


pooledvar <- function(perf.measures, N, ...) {
  pm <- unlist(perf.measures)
  pm
  ## TODO: Extract sample size for each cluster and apply corresponding to the right perf.measures
  
}