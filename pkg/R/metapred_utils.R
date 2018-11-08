# Internal function for centering
# x numeric vector to be centered in data sets.
# center.in indices of data sets.
center <- function(x, center.in) {
  if (length(center.in) != length(x))
    stop("length(center.in) should match length(x).")
  for (trial in sort(unique(center.in))) {
    selection.id <- center.in == trial
    selection <- x[selection.id]
    x[selection.id] <- selection - mean(selection, na.rm = T)
  }
  x
}

# Old. Works but choice of variables is limited and does not make a lot of sense.
# # Centers data within studies / clusters
# # Note that the center indicator is also centered. This is not intended, but can be worked around.
# # data data.frame. data set.
# # center.in numeric vector corresponding to cluster indices.
# # center.1st logical. Should the 1st variable (response) be centered
# # center.rest logical. Should the other variables (predictors) be centered?
# centerData <- function(data, center.in, center.1st = FALSE, center.rest = FALSE) {
#   if (!is.data.frame(data) && !is.matrix(data))
#     stop("data should be a data.frame or matrix.")
#   if (length(center.in) != nrow(data))
#     stop("length(center.in) should match nrow(data).")
#   if (isTRUE(center.1st))
#     data[ , 1] <- center(data[ , 1], center.in)
#   if (isTRUE(center.rest) && ncol(data) > 2)
#     for (col in 2:ncol(data))
#       data[ , col] <- center(data[ , col], center.in)
#     data
# }

# Centers data within studies / clusters
# # Note that the center indicator is also centered. This is not intended, but can be worked around. Sadly, this means
# that a categorical cluster variable will also be tried to be centered, causing an error.
# data data.frame. data set.
# center.i numeric vector corresponding to cluster indicators.
# center.which numeric or integer Which variables should be centered? Defaults to all except for the first column (assumed to be
# the outcome)
# centerData.old <- function(data, center.i, center.which = NULL) {
#   if (!is.data.frame(data) && !is.matrix(data))
#     stop("data should be a data.frame or matrix.")
#   if (length(center.i) != nrow(data))
#     stop("length(center.i) should match nrow(data).")
#   if (is.null(center.which))
#     center.which <- seq_len(ncol(data))[-1]
#   if (isTRUE(any(!!center.which))) {
#     if (all(center.which <= ncol(data)) && all(center.which > 0))
#       for (col in center.which)
#         data[ , col] <- center(data[ , col], center.i) 
#       else stop("center.which should indicate (numeric values of) columns of data set that are to be centered. Use 0 or FALSE for none.")  
#   }
#   return(data)
# }


# Deprecated
# Centers data within studies / clusters
# data data.frame. data set.
# cluster.var name of variable in data, vector corresponding to cluster indicators. 
# Overrides cluster.vec If both are NULL, all are in the same cluster. 
# cluster.vec vector of cluster indicators. Overridden by cluster.var. If both are NULL, all are in the same cluster.
# center.which numeric or integer Which variables should be centered? Defaults to all numeric except for the first 
# column (assumed to be the outcome)
# centerData <- function(data, cluster.var = NULL, cluster.vec = NULL, center.which = NULL) {
#   if (!is.data.frame(data) && !is.matrix(data))
#     stop("data should be a data.frame or matrix.")
#   
#   if (!is.null(cluster.var)) {
#     cluster.vec <- data[ , cluster.var]
#     center.col <- which(colnames(data) == cluster.var)
#   } else {
#     center.col <- 0
#     if (is.null(cluster.vec)) cluster.vec <- rep(1, nrow(data))
#   }
# 
#   if (length(cluster.vec) != nrow(data))
#     stop("length(cluster.vec) should match nrow(data).")
#   
#   # if (is.null(center.which)) # old, gives error when categorical variables are supplied
#   #   center.which <- seq_len(ncol(data))[-1]
#   
#   if (is.null(center.which)) { # new, experimental
#     center.which <- unlist(lapply(data, is.numeric))
#     center.which[1] <- FALSE
#     center.which <- which(center.which)
#   }
#   
#   if (isTRUE(any(!!center.which))) {
#     if (all(center.which <= ncol(data)) && all(center.which > 0)) {
#       for (col in center.which)
#         if (col != center.col) 
#           data[ , col] <- center(x = data[ , col], center.in = cluster.vec) 
#     } else stop("center.which should indicate (numeric values of) columns of data set that are to be centered. Use 0 or FALSE for none.")  
#   }
#   return(data)
# }

# Center covariates within clusters
# data data.frame
# y.name character, name of outcome variable
# cluster.name character, name of cluster variable.
centerCovs <- function(data, y.name, cluster.name) {
  to.center <- which((!(colnames(data) == cluster.name | colnames(data) == y.name) ) & apply(data, 2, is.numeric))
  cluster.vec <- data[ , cluster.name]
  
  for (col in to.center)
    data[ , col] <- center(data[ , col], center.in = cluster.vec)
  
  data
}


# Deprecated
# coerces data set to data.list
# data data.frame. data set.
# strata.i numeric. stratum indicators.
# Returns data.list
# asDataList <- function(data, strata.i) {
#   data.list <- list()
#   strata <- sort(unique(strata.i))
#   for (i in 1:length(strata))
#     data.list[[i]] <- data[strata.i == strata[i], ]
#   names(data.list) <- strata
#   data.list
# }

# Deprecated
# gets only the relevant data from a data.list
# data.list list of data sets
# ccs numeric. covariate column selection
# cl numeric. indices of clusters to be selected.
# returns data.list.
# getDataList <- function(data.list, ccs, cl) 
#   lapply(data.list[cl], getData, predictors = ccs)


# These functions are used for making various names. Any change should be made here, such that
# these functions can also be used to retrieve objects from a list.
# st.u vector. unique indices or names of clusters.
# f numeric. fold index.
# type character. type of cv (just a chosen name)
# data.list list of data sets.
# data data.frame.
# covariate.columns indices of covariate columns.
# All return character.
# getFoldName <- function(st.u, f = NULL, type = NULL)
#   paste(getclName(st.u = st.u), getcvName(f = f, type = type), sep = if (length(f) > 0 || length(type) > 0) ". " else "")

# For l10, fixed and successive, the validation folds are unique. Bootstrap reuses them, and must add an iteration number.
getFoldName <- function(st.u, f = NULL, type = NULL) {
  if (isTRUE(type == "l1o") || isTRUE(type == "leaveOneOut") || isTRUE(type == "fixed") || isTRUE(type == "successive"))
    paste(getclName(st.u = st.u))
  else 
    paste(getclName(st.u = st.u), getcvName(f = f, type = type), sep = if (length(f) > 0 || length(type) > 0) ". " else "")
}
  

# getclName <- function(st.u)
#   paste("cl", toString(st.u), sep = " ")

getclName <- function(st.u)
  paste(toString(st.u), sep = " ")

getcvName <- function(f, type = NULL)
  paste(type, f, sep = " ")

# Deprecated
# getCovariateNames <- function(data.list, covariate.columns) {
#   if (!length(covariate.columns))  return(NULL)
#   if (!is.null(colnames(data.list[[1]]))) return(colnames(data.list[[1]])[covariate.columns])
#   warning("Covariate names could not be found.")
#   return(paste("X", covariate.columns - 1, sep = ""))
# }
# Deprecated
# getModelName <- function(data.list, covariate.columns)
#   if (is.null(name <- getCovariateNames(data.list, covariate.columns) ))
#     return("intercept only") else return(toString(name))

getStepName <- function(x)
  paste("s", x, sep = "")

# Deprecated
# getPredictorNames <- function(f, data)
#   colnames(stats::model.frame.default(formula = f, data = data))[-1]

# Deprecated
# This one is mostly important for getting the intercept-only formula:
# data data.frame
# predictors indices of predictors
# returns formula
# getFormula <- function(data, predictors = NULL) {
#   if (length(predictors) == 0 || identical(predictors, 0)) {
#     f <- stats::formula(data)
#     f[3] <- 1 # [3] is the right hand side.
#   } else f <- stats::formula(data[ , c(1, predictors)])
#   stats::as.formula(f)
# }

# Deprecated
# Gets only the requested data, as well as the response variable.
# data data.frame or matrix
# pred.indices indexes of preditor columns
# Returns data.frame.
# getData <- function(data, pred.indices) {
#   if (!is.data.frame(data) && !is.matrix(data))
#     stop("data must be a data.frame or matrix.")
#   sel <- c(1, pred.indices)
#   d <- data.frame(data[ , sel])
#   names(d) <- names(data)[sel]
#   d
# }

# Deprecated
# gets formula of for data set. Replacement of formula.data.frame, which returns bogus for
# 1 column data.frames.
# data data.frame
# Returns formula.
# getFullFormula <- function(data) {
#   data <- as.data.frame(data)
#   # if (!is.data.frame(data) && !is.matrix(data))
#   #   stop("data must be a data.frame or matrix.")
#   
#   if (identical(ncol(data), 1L))
#   {
#     xnames <- names(data)[-1]
#     if (identical(length(xnames), 0L))
#       xnames <- "1"
#     left <- paste(names(data)[1], "~")
#     return(paste(left, xnames))
#     
#   } else return(stats::formula(data))
# }

getCoefs  <- function(fit, ...) {
  if (inherits(fit, "multinom"))
    return(coefMultinom(fit, ...))
  else return(coef(fit))
}
# Needs some work. Should also return some coefficient names.
coefMultinom <- function(fit, ...)
  as.vector(t(coef(fit)))

# Perhaps unnecessary:
getVars   <- function(fit, ...) diag(vcov(fit))
getCoVars <- function(fit, ...) vcov(fit)
getSE     <- function(fit, ...) sqrt(getVars(fit))

# Deprecated
# Coerces l to one string
# l list or vector of strings.
# returns 1 string, character.
# oneStr <- function(l, sep = "") {
#   if (length(l) > 0) out <- l[[1]] else return("")
#   
#   if (length(l) > 1) {
#     for (i in 2:length(l)) out <- paste(out, l[[i]], sep = sep)
#   }
#   return(out)
# }

# Deprecated
# Necessary for making a formula from a model.frame.
# f formula
# returns: formula without backticks.
# removeFormulaBackticks <- function(f)
# {
#   g <- gsub("`", "", f)
#   h <- oneStr(list(g[[2]], " ~ ", g[[3]]), sep = "")
#   stats::update.formula(f, h)
# }

# Deprecated
# f formula
# p.name character. name of predictor
# Returns: formula
# removePredictor <- function(f, p.name)
#   stats::update.formula(f, paste(". ~ . -", p.name) )

# Deprecated
# f formula
# p.name list or vector of character names of predictors
# Returns: formula
# removePredictors <- function(f, p.names)
#   removePredictor(f, paste(oneStr(unlist(p.names), sep = " - "), sep = " - "))

### The following functions are for generating the folst.u for the cross-validation in metapred
# st.u Numeric or character vector. Unique names of the strata / clusters.
# k Numeric. Differs per function
#   bootstrap: number of bootstraps
#   fixed: indices of data sets for validation.
#   successive (still a hidden function): Number of validation sets.
#   leaveOneOut = l1o: iecv: internal-external cross-validation of data sets/clusters.
leaveOneOut <- l1o <- function(st.u, ...) {
  st.u <- sort(st.u)
  if (length(st.u) < 2)
    stop("iecv not possible for fewer than 2 strata.")
  indexes <- seq_len(length(st.u))
  out <- list(dev = list(), dev.i = list(), val = as.list(st.u[indexes]), val.i = as.list(indexes))
  
  for (i in indexes) {
    out$dev[[i]]   <- st.u[-indexes[i]]
    out$dev.i[[i]] <- indexes[-i]
  }
  
  out
}

bs <- bootstrap <- function(st.u, k = NULL, ...) {
  st.u <- sort(st.u)
  if (is.null(k))
    k <- 200
  if (length(st.u) < 2)
    stop("Bootstrapping data sets is impossible for < 2 data sets.")
  dev <- dev.i <- val <- val.i <- list()
  
  i <- 1
  while (length(dev) < k) {
    indexes <- sample(1:length(st.u), replace = T)
    if (length(unique(indexes)) >= length(st.u))
      next
    dev[[i]] <- st.u[indexes]
    val[[i]] <- st.u[-indexes]
    dev.i[[i]] <- indexes
    val.i[[i]] <- seq_len(length(st.u))[-indexes]
    i <- i + 1
  }
  list(dev = dev, dev.i = dev.i , val = val, val.i = val.i)
}

fixed <- function(st.u, k = NULL, ...) {
  st.u <- sort(st.u)
  if (length(st.u) < 2)
    stop("Selecting a validation set is impossible for < 2 data sets.")
  if (is.null(k))
    k <- length(st.u)
  indexes <- seq_len(length(st.u))
  list(dev = list(st.u[-k]), dev.i = list(indexes[-k]), val = list(st.u[k]), val.i = list(indexes[k]))
}

ws <- within.sample <- function(st.u, ...) { # Necessary for testing recalibration functions.
  st.u <- sort(st.u)
  if (length(st.u) < 2)
    stop("Selecting a validation set is impossible for < 2 data sets.")
  indexes <- seq_len(length(st.u))
  list(dev = as.list(st.u[indexes]), dev.i = as.list(indexes), val = as.list(st.u[indexes]), val.i = as.list(indexes))
}

successive <- function(st.u, k = NULL, ...) {
  st.u <- sort(st.u)
  if (is.null(k)) k <- 1
  k <- as.integer(k)
  if (k < 1) stop("k must be >= 1")
  if (length(st.u) < (k + 1)) stop("Cross-validation requires k + 1 data sets, where k is the number of test sets.")
  out <- list(dev = list(), dev.i = list(), val = list(), val.i = list())
  
  for (i in seq_len(length(st.u) - k) ) {
    sel <- 1:i
    out$dev.i[[i]] <- sel
    out$dev[[i]] <- st.u[sel]
  }
  for (i in seq_len(length(out$dev)) ) {
    sel <- (i + 1):(i + k)
    out$val.i[[i]] <- sel
    out$val[[i]] <- st.u[sel]
  }
  out
}

# Remove observations (rows) that have at least one missing value.
# Necessary for metapred(), to ensure that the same observations are used
# df data.frame
# Returns: data.frame
remove.na.obs <- function(df) 
  df[apply(df, 1, function(df) sum(is.na(df))) == 0, ]

# Gets the predict method.
# fit Model fit object.
# two.stage logical. Is the model a two-stage model?
# predFUN Optional function, which is immediately returned
# ... For compatibility only.
getPredictMethod <- function(fit, two.stage = TRUE, predFUN = NULL, ...) {
  # A user written function may be supplied:
  if (!is.null(predFUN)) {
    if (is.function(predFUN)) {
      return(predFUN)
    } else return(get(as.character(predFUN), mode = "function"))
  }

  # If two-stage, the fit is used only to extract the link function.
  # If one-stage, fit's prediction method may be used.
  if (two.stage) {# Preferably mp.cv.dev should not be here. But it currently has to.
    if (inherits(fit, c("glm", "lm", "mp.cv.dev"))) 
      return(predictGLM)
    else stop("No prediction method has been implemented for this model type yet for two-stage
              meta-analysis. You may supply one with the predFUN argument.")
  } else return(predict)
}

# Prediction function for two-stage metapred GLM objects
# object glm model fit object
# newdata newdata to predict for, of class "data.frame"
# b vector of coef. Overrides coef of object
# f formula used for selecting relevant variables from newdata. Overrides object
# ... For compatibility only.
# Returns vector of predicted values.
predictGLM <- function(object, newdata, b = NULL, f = NULL, type = "response", ...) {
  if (is.null(b)) b <- coef(object)
  if (is.null(f)) f <- formula(object)
  X <- model.matrix(stats::as.formula(f), data = newdata)
  X <<- X
  b <<- b
  
  lp <- X %*% b
  
  if (identical(type, "response")) {
    if (is.null(fam <- family(object)))
      return(lp)
    else
      return(fam$linkinv(lp))
  } else if (identical(type, "link"))
    return(lp)
  # if (is.null(object$family)) lp
  # else object$family$linkinv(lp)
}

# Univariate Random Effects Meta-Analysis
# coefficients data.frame or matrix, containing coef
# variances data.frame or matrix, containing variances
# method Method for meta-analysis.
# ... Optional arguments for rma().
#' @importFrom metafor rma
urma <- function(coefficients, variances, method = "DL", ...)
{
  if (!(is.data.frame(coefficients) || is.matrix(coefficients)) || !(is.data.frame(variances) || is.matrix(variances)) )
    stop("coefficients and variances must both be a data.frame or matrix.")
  if (!identical(dim(coefficients), dim(variances)))
    stop("coefficients and variances must have the same dimensions.")
  
  meta.b <- meta.se <- rep(NA, ncol(coefficients))
  for (col in 1:ncol(coefficients)) {
    r <- metafor::rma(coefficients[ , col] , variances[ , col], method = method, ...)
    meta.b[col]  <- r$beta
    meta.se[col] <- r$se
  }
  
  meta.v <- meta.se^2
  
  names(meta.b) <- names(meta.v) <- names(meta.se) <- colnames(coefficients)
  list(coefficients = meta.b, variances = meta.v, se = meta.se)
}


which.abs.min <- function(x) 
  which.min(abs(x))

which.abs.max <- function(x)
  which.max(abs(x))

which.1 <- function(x)
  which.abs.min(x - 1)

#' @author Valentijn de Jong
#' @method unlist   listofperf
#' @export
unlist.listofperf <- function(x, ...) 
  sapply(x, `[[`, 1)


# Safe way of getting a function. 
# For some reason, this function can find functions that match.fun cannot.
# x function or character name thereof
# Returns function or error.
get.function <- function(x, ...) {
  if (is.function(x))
    return(x)
  else
    return(get(as.character(x), mode = "function"))
}