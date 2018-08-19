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

#' @importFrom pROC auc
auc <- function(p, y, ...) {
  if (is.matrix(p))
    p <- p[, 1]
  if (is.matrix(y))
    y <- y[, 1]
  pROC::auc(response = y, predictor = p)
}

calibration.intercept <- cal.int <- function(p, y, estFUN, family, ...)
  pred.recal(p = p, y = y, estFUN = estFUN, family = family, which = "intercept")

# Slope.only is a trick to make this functin work for metapred.
# Slope.only should otherwise always be false! Also: this messes up the variances,
# making meta-analysis impossible!
# multiplicative slope!
calibration.slope <- cal.slope <- function(p, y, estFUN, family, slope.only = TRUE, ...) {
  # refit <- pred.recal(p = p, y = y, estFUN = estFUN, family = family, which = "slope")
  # if (slope.only) {
  #   refit[[1]] <- refit[[1]][[2]]
  # }
  # refit
  
  refit <- pred.recal(p = p, y = y, estFUN = estFUN, family = family, which = "slope")
  if (slope.only) {
    refit[[1]] <- refit[[1]][2]
    refit$variances <- variances(refit)[2]
  }
  refit
}

# additive slope!
calibration.add.slope <- cal.add.slope <- function(p, y, estFUN, family, slope.only = TRUE, ...)  {
  
  refit <- pred.recal(p = p, y = y, estFUN = estFUN, family = family, which = "add.slope")
  if (slope.only) {
    refit[[1]] <- refit[[1]][2]
    refit$variances <- variances(refit)[2]
  }
  refit
}
  






############################## Heterogeneity, generalizability, pooled performance functions ###############################
# ### By convention, all generalizability measures:
# # Required arguments:
# x       list of class "listofperf", list of performance in different strata. 
#           Note that it has its own unlist method. Practically all (except plot) should call unlist first!
# ...     for compatibility.
#
# # Possible arguments, that are always passed through successfully:
# N       scalar integer, total sample size
# n       named numeric/integer vector, sample size per strata. Has same length as x. (OR CHANGE  TO LIST AS WELL??)
# data    Full data set. Use only if other arguments are not enough.
# coef    data.frame containing coefficients of stratified models. Rows are strata, columns are coefs.
# coef.se data.frame containing se of coefficients of stratified models. Rows are strata, columns are coefs.
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
  x <- unlist(x) ### back to unlist? now that a method has been implemented?!
  cv <- sd(x)/mean(x)
  if (isTRUE(abs)) abs(cv) else cv
}

coef.var.mean <- function(x, abs = TRUE, ...)  {
  x <- unlist(x)
  coef.var(x, abs = abs) + if (abs) abs(mean(x)) else mean(x)
}
  

# Measure 2 (?): GINI coefficient
# No code or import necessary
# GiniMd(x, na.rm = FALSE)


weighted.abs.mean <- function(x, n, ...) 
  abs.mean(x <- unlist(x) * sqrt(n - 1)) / sum(sqrt(n - 1))


pooled.var <- function(x, n, ...) {
  x <- unlist(x)
  
  ## TODO: Extract sample size for each cluster and apply corresponding to the right performance measures
  ## TODO: use rubins rules.
}

rubins.rules <- function(x, n, ...) {
  x <- unlist(x)
  x + var(x) * (1 + 1/n)
}

# squared.diff #a penalty equal to the mean squared differences 
squared.diff <- function(x, ...) {
  x <- unlist(x)
  mse(x, mean(x))
}

# New, for auc and intercept, and slope tbi
#  Forest plot of list of performance measures. Currently only works for auc from the pROC package.
#' @importFrom metafor rma.uni
#' @export
plot.listofperf <- function(x, ...) { # xlab tbi from perfFUN
  xlab <- paste(list(...)$perfFUN.name, "in validation strata")
  
  if (is.null(names(x))) # The # is to show users that the numbers are not their own. (no longer necessary)
    names(x) <- paste("#", seq_along(x), sep = "") 
  
  z <- ci.listofperf(object = x, ...)
    
  if (inherits(x[[1]], "auc")) { # To be replaced by child function.
    vm <- valmeta(measure = "cstat", cstat = z$theta, cstat.95CI = z[, c("theta.ci.lb", "theta.ci.ub")])
    est <- vm$est
    ci.lb <- vm$ci.lb
    ci.ub <- vm$ci.ub
    pi.lb <- vm$pi.lb
    pi.ub <- vm$pi.ub
  } else if (inherits(x[[1]], "lm")) {
    mf <- predict(metafor::rma.uni(yi = sapply(x, coef), vi = sapply(x, variances))) # NOTE: DOES IT USE t or normal dist???
    est <- mf$pred
    ci.lb <- mf$ci.lb
    ci.ub <- mf$ci.ub
    pi.lb <- mf$cr.lb
    pi.ub <- mf$cr.ub
  }

  
  fp <- metamisc::forest(theta       = z$theta,
                         theta.ci.lb = z$theta.ci.lb,
                         theta.ci.ub = z$theta.ci.ub,
                         theta.slab  = names(x),
                         theta.summary       = est,
                         theta.summary.ci.lb = ci.lb,
                         theta.summary.ci.ub = ci.ub,
                         theta.summary.pi.lb = pi.lb,
                         theta.summary.pi.ub = pi.ub,
                         xlab  = xlab,
                         ...)
  plot(fp)
  NaN
}

# Old, for auc only
# #  Forest plot of list of performance measures. Currently only works for auc from the pROC package.
# #' @importFrom pROC ci
# #' @export
# plot.listofperf <- function(x, ...) { # xlab tbi from perfFUN
#   xlab <- paste(list(...)$perfFUN.name, "in validation strata")
#   
#   if (is.null(names(x))) # The # is to show users that the numbers are not their own. (no longer necessary)
#     names(x) <- paste("#", seq_along(x), sep = "") 
#   
#   if (in)
# 
#   z <- lapply(x, pROC::ci)  
#   theta       <- sapply(z, `[[`, 2)
#   theta.ci.lb <- sapply(z, `[[`, 1)
#   theta.ci.ub <- sapply(z, `[[`, 3)
#   
#   vm <- valmeta(measure = "cstat", cstat = theta, cstat.95CI = data.frame(theta.ci.lb, theta.ci.ub))
#   
#   fp <- metamisc::forest(theta       = theta,
#                          theta.ci.lb = theta.ci.lb,
#                          theta.ci.ub = theta.ci.ub,
#                          theta.slab  = names(x),
#                          theta.summary = vm$est,
#                          theta.summary.ci.lb = vm$ci.lb,
#                          theta.summary.ci.ub = vm$ci.ub,
#                          theta.summary.pi.lb = vm$pi.lb,
#                          theta.summary.pi.ub = vm$pi.ub,
#                          sort  = FALSE, 
#                          xlab  = xlab,
#                          ...)
#   plot(fp)
#   NaN
# }

#' ##############################               OLD, x is still vector here                     ###############################





#' ##############################                Performance / error functions                  ###############################
#' # ### By convention, all performance measures:
#' # # Arguments:
#' # p     numeric vector of predicted probabilities
#' # y     numeric/integer vector of observed outcome, of same length as p.
#' # ...   for future compatibility
#' # 
#' # # Return:
#' # numeric of length 1.
#' ##############################                
#' 
#' # Error function: Mean Squared Error
#' mse <- brier <- function(p, y, ...) mean((p - y)^2)
#' 
#' rmse <- function(p, y, ...)
#'   sqrt(mse(p = p, y = y, ...))
#' 
#' 
#' 
#' 
#' # Error function: Variance of prediction error
#' var.e <- function(p, y, ...) var(p - y)
#' 
#' 
#' 
#' # Measure 1: Coefficient of variation of prediction error.
#' coef.var.pred <- function(p, y, abs = TRUE, ...)
#'   coef.var(x = p - y, abs = abs, ...) 
#' 
# #' #' @importFrom pROC auc
#' auc <- function(p, y, ...) {
#'   if (is.matrix(p))
# #'     p <- p[, 1]
#'   if (is.matrix(y))
# #'     y <- y[, 1]
#'   pROC::auc(response = y, predictor = p)
#' }
#' 
#' 
#' 
#' 
#' 
#' ############################## Heterogeneity, generalizability, pooled performance functions ###############################
#' # ### By convention, all generalizability measures:
#' # # Required arguments:
#' # x       numeric vector of performance in different strata
#' # ...     for compatibility.
#' #
#' # # Possible arguments, that are always passed through successfully:
#' # N       scalar integer, total sample size
#' # n       named numeric/integer vector, sample size per strata. Has same length as x.
#' # data    Full data set. Use only if other arguments are not enough.
#' # coef    data.frame containing coefficients of stratified models. Rows are strata, columns are coefs.
#' # coef.se data.frame containing se of coefficients of stratified models. Rows are strata, columns are coefs.
#' #
#' # # Possible arguments that are passed through only for certain measures:
#' # se, TBI
#' # 
#' # # Return:
#' # numeric of length 1.
#' ##############################   
#' 
#' # Measure 0: mean
#' abs.mean <- function(x, ...)
#'   abs(mean(unlist(x)))
#' 
#' 
#' # Measure 1: Coefficient of variation (=scaled sd)
#' # In general sense, abs needs not be TRUE, but for metapred it should,
#' # such that higher values are worse performance.
#' coef.var <- function(x, abs = TRUE, ...) {
#'   cv <- sd(x)/mean(x)
#'   if (isTRUE(abs)) abs(cv) else cv
#' }
#' 
#' coef.var.mean <- function(x, abs = TRUE, ...) 
#'   coef.var(x, abs = abs) + if (abs) abs(mean(x)) else mean(x)
#' 
#' # Measure 2 (?): GINI coefficient
#' # No code or import necessary
#' # GiniMd(x, na.rm = FALSE)
#' 
#' 
#' weighted.abs.mean <- function(x, n, ...)
#'   abs.mean(x * sqrt(n - 1)) / sum(sqrt(n - 1))
#' 
#' 
#' pooled.var <- function(x, n, ...) {
#'   pm <- unlist(x)
#'   pm
#'   ## TODO: Extract sample size for each cluster and apply corresponding to the right performance measures
#'   ## TODO: use rubins rules.
#' }
#' 
#' rubins.rules <- function(x, n, ...)
#'   x + var(x) * (1 + 1/n)
#' 
#' plotauc <- function(x, ...) {
#'   print(x)
#'   print(ci(x))
#'   NULL
#' }
#' 
#' 
#' # squared.diff #a penalty equal to the mean squared differences 
#' squared.diff <- function(x, ...)
#'   mse(x, mean(x))
#' 
# #' #' @importFrom pROC ci
#' forest.listofperf <- function(x, title, ...) {
#'   z <- lapply(x, pROC::ci) 
# #'   metamisc::forest(theta       = sapply(z, `[[`, 2),
# #'                    theta.ci.lb = sapply(z, `[[`, 1),
# #'                    theta.ci.ub = sapply(z, `[[`, 3),
# #'                    theta.slab  = names(x),
#'                    title = title,
#'                    sort  = FALSE,
#'                    xlab  = "AUC")
#' }
#' 
