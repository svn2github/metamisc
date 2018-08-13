##############################                Performance / error functions                  ###############################
# ### By convention, all performance measures:
# # Arguments:
# p     numeric vector of predicted probabilities
# y     numeric/integer vector of observed outcome
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
vare <- function(p, y, ...) var(p - y)



# Measure 1: Coefficient of variation of prediction error.
coefVarPred <- function(p, y, data = NULL, abs = TRUE, ...)
  coefVar(x = p - y, abs = abs, ...) 

############################## Heterogeneity, generalizability, pooled performance functions ###############################
# ### By convention, all generalizability measures:
# # Required arguments:
# x     numeric vector of performance in different strata
# ...   for compatibility.
#
# # Possible arguments:
# N     integer, total sample size
# n     numeric/integer vector, sample size per strata. Must be same length as x
# data  Full data set. Use only if other arguments are not enough.
# 
# # Return:
# numeric of length 1.
##############################   

# Measure 0: mean
absmean <- function(x, ...)
  abs(mean(unlist(x)))


# Measure 1: Coefficient of variation (=scaled sd)
# In general sense, abs needs not be TRUE, but for metapred it should,
# such that higher values are worse performance.
coefVar <- function(x, abs = TRUE, ...) {
  cv <- sd(x)/mean(x)
  if (isTRUE(abs)) abs(cv) else cv
}

coefVarMean <- function(x, abs = TRUE, ...) 
  coefVar(x, abs = abs) + if (abs) abs(mean(x)) else mean(x)


pooledvar <- function(x, n, ...) {
  pm <- unlist(x)
  pm
  ## TODO: Extract sample size for each cluster and apply corresponding to the right perf.measures
  
}