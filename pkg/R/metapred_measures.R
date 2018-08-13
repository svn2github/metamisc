##############################                Performance / error functions                  ###############################
# ### By convention, all performance measures:
# # Arguments:
# p     numeric vector of predicted probabilities
# y     numeric/integer vector of observed outcome, of same length as p.
# ...   for future compatibility
# 
# # Return:
# numeric of length 1.
##############################                

# Error function: Mean Squared Error
mse <- brier <- function(p, y, ...) mean((p - y)^2)

rmse <- function(p, y, ...)
  sqrt(mse(p = p, y = y, ...))




# Error function: Variance of prediction error
var.e <- function(p, y, ...) var(p - y)



# Measure 1: Coefficient of variation of prediction error.
coef.var.pred <- function(p, y, abs = TRUE, ...)
  coef.var(x = p - y, abs = abs, ...) 

############################## Heterogeneity, generalizability, pooled performance functions ###############################
# ### By convention, all generalizability measures:
# # Required arguments:
# x     numeric vector of performance in different strata
# ...   for compatibility.
#
# # Possible arguments, that are always passed through successfully:
# N     scalar integer, total sample size
# n     named numeric/integer vector, sample size per strata. Has same length as x.
# data  Full data set. Use only if other arguments are not enough.
#
# # Possible arguments that are passed through only for certain measures:
# se, TBI
# 
# # Return:
# numeric of length 1.
##############################   

# Measure 0: mean
abs.mean <- function(x, ...)
  abs(mean(unlist(x)))


# Measure 1: Coefficient of variation (=scaled sd)
# In general sense, abs needs not be TRUE, but for metapred it should,
# such that higher values are worse performance.
coef.var <- function(x, abs = TRUE, ...) {
  cv <- sd(x)/mean(x)
  if (isTRUE(abs)) abs(cv) else cv
}

coef.var.mean <- function(x, abs = TRUE, ...) 
  coef.var(x, abs = abs) + if (abs) abs(mean(x)) else mean(x)

# Measure 2 (?): GINI coefficient
# No code or import necessary
# GiniMd(x, na.rm = FALSE)


weighted.abs.mean <- function(x, n, ...)
  abs.mean(x * sqrt(n - 1)) / sum(sqrt(n - 1))


pooled.var <- function(x, n, ...) {
  pm <- unlist(x)
  pm
  ## TODO: Extract sample size for each cluster and apply corresponding to the right performance measures
  ## TODO: use rubins rules.
}

rubins.rules <- function(x, ...)
  x + var(x) * (1 + 1/length(x))

# squared.diff #a penalty equal to the mean squared differences between coefficients. Alternatively