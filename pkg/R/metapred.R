# TODO  1
# @Valentijn: perhaps better to combine param perfFUN, genFUN and selFUN into one parameter with 3 or 4 distinct options
# @Valentijn: I suggest to omit intercept recalibration. The intercept issue can be addressed directly by specifying distinct 
# error functions that do or do not account for mis-calibration in intercept term
# Add summary.metapred
# categorical variables
#   + categorical variables as one or as multiple.
# add function: areTRUE. st.i == cl may otherwise lead to problems, if st.i has NA

# TODO  2
# add more output for metapred
# tol
# centering: now only works for numeric variables
# recal.int + interaction gives error: 
# mp <- metapred(data = DVTipd.reordered, 
#                            strata = "cluster",
#                            formula = dvt ~ ddimdich * histdvt, 
#                            family = binomial, 
#                            recal.int = T)
# metapred.summaries docs (see glm.summaries {stats}), including subset
# add stepnumbers to print.metapred or print.fit
# Write test for new centerCovs function.
# add fitted.metapred()

# TODO 3
# variances for intercept recalibration.
# One-stage MA: predStep: is  type = response correct?
# perf method
# penalties
# performance measurement, current is Brier for binomial
# add more options to penalty step
# , cl.name in modelStep
# Remove Reduce() from perfStep, and make perfStep() compatible with multiple data sets (for cvFUN = fixed, or cvFUN = bs)

# Changed:
# Added is.metapred()
# Added more robust tests for stepwise.
# Automatically remove observations with missing values.
# @Valentijn: I would center covariates by default, this often helps to improve generalizability and also speeds estimation
#   Valentijn: Done.
# @Valentijn: You can change default of meta.method to DL. This is a lot faster and has limited implications on estimated means.
#   Valentijn: Done.

###### Outline
### Top: high-level functions
#
# metapred
# 
# mp.fit
# 
# mp.step + mp.global
# 
# mp.cv
# 
# mp.cv.val
# mp.cv.dev	+ mp.cv.recal
# 
# mp.cv.meta.fit + mp.recal.meta.fit
# 
# mp.meta.fit
#
# mp.stratified.fit
#
# mp.stratum.fit
#
### bottom : low level functions
###### 


#' Generalized Stepwise Regression for prediction models in clustered data
#'
#' Generalized stepwise regression for obtaining a prediction model with adequate performance across data sets. Requires
#' data from individuals in multiple studies.
#' 
#' @author Valentijn de Jong
#' 
#' @references Debray TPA, Moons KGM, Ahmed I, Koffijberg H, Riley RD. A framework for developing, implementing, 
#' and evaluating clinical prediction models in an individual participant data meta-analysis. 
#' \emph{Stat Med}. 2013;32(18):3158-80. 
#'
#' @param data data.frame containing the data.
#' @param strata Name of the strata (e.g. studies or clusters) variable, as character. Used for two-stage MA only.
#' @param formula \code{formula} of the first model to be evaluated. \code{metapred} will start at \code{formula} and update it
#' using terms of \code{scope}. Defaults to full main effects model, where the first column in \code{data} is assumed to be
#' the outcome.
#' \code{metapred} assumes the first column in the data set is the outcome, and all remaining columns
#' (except \code{strata}) are predictors. See \link[stats]{formula} for formulas in general.
#' @param estFUN Function for estimating the model in the first stage. Currently "lm" and "glm" are supported.
#' @param scope \code{formula}. The difference between \code{formula} and \code{scope} defines the range of models examined in the 
#' stepwise search. Defaults to NULL, which leads to the intercept-only model. If \code{scope} is nested in \code{formula}, 
#' this implies backwards selection will be applied (default). If \code{scope} is nested in \code{formula}, this implies forward 
#' selection will be applied. If equal, no stepwise selection is applied. 
#' @param retest Logical. Should added (removed) terms be retested for removal (addition)? \code{TRUE} implies bi-directional 
#' stepwise search.
#' @param max.steps Integer. Maximum number of steps (additions or removals of terms) to take. Defaults to 1000, which is
#' essentially as many as it takes. 0 implies no stepwise process.
#' @param center logical. Should numeric covariates be centered?
#' @param recal.int Logical. Should the intercept be recalibrated?
#' @param cvFUN Cross-validation method, on the study (i.e. cluster or stratum) level. "
#' l1o" for leave-one-out cross-validation (default). "bootstrap" for bootstrap. Or "fixed", for one or more data sets
#' which are only used for validation. A user written function may be supplied as well.
#' @param cv.k Parameter for cvFUN. For \code{cvFUN="bootstrap"}, this is the number of bootstraps. For \code{cvFUN="fixed"}, 
#' this is a vector of the indices of the (sorted) data sets. Not used for \code{cvFUN="l1o"}.
#' @param metaFUN Function for computing the meta-analytic coefficient estimates in two-stage MA. 
#' Default: \link[metafor]{rma.uni}
#' from the metafor package is used. Default settings are univariate random effects, estimated with "DL". Method can be
#' passed trough the \code{meta.method} argument.
#' @param meta.method Name of method for meta-analysis. Default is "DL". For more options see \link[metafor]{rma.uni}.
#' @param predFUN Function for predicting new values. Defaults to the appropriate link functions for two-stage MA where
#' \code{glm()} or \code{lm()} is used in the first stage. For one-stage models \code{predict()} is used.
#' @param perfFUN Function for computing the performance of the prediction models. Default: mean squared error 
#' (\code{perfFUN="mse"}).
#' Other options are \code{"vare"}.
#' @param genFUN Function computing generalizability measure using the performance measures. Default: (absolute) mean  
#' (\code{genFUN="absmean"}). Choose \code{squareddiff} for a penalty equal to the mean squared differences between 
#' coefficients. Alternatively, choose \code{pooledvar} for a weighted average of variance terms.
#' @param selFUN Function for selecting the best method. Default: lowest value for \code{genFUN}. Should be set to
#' "which.max" if high values for \code{genFUN} indicate a good model.
#' @param ... To pass arguments to estFUN (e.g. family = "binomial"), or other methods.
#'
#' @return \code{metapred} A list of class \code{metapred}, containing the final coefficients in \code{coef}, and the stepwise
#' tree of estimates of the coefficients \code{(coef)}, performance measures \code{(perf)}, generalizability measures
#' \code{(gen)} in \code{stepwise}, and more.
#'
#' @examples 
#' data(DVTipd)
#' DVTipd$cluster <- letters[1:4] # Add a fictional clustering to the data.
#' #f <- dvt ~ sex + vein + malign
#' f <- dvt ~ histdvt + ddimdich
#' # Leaving scope empty implies backwards selection
#' mp <- metamisc::metapred(DVTipd, strata = "cluster", formula = f, family = binomial)
#' 
#'\dontrun{
#' # Some additional examples:
#' metapred(DVTipd, strata = "cluster", formula = dvt ~ 1, scope = f, family = binomial) # Forwards
#' metapred(DVTipd, strata = "cluster", formula = f, scope = f, family = binomial) # no selection
#' metapred(DVTipd, strata = "cluster", formula = f, max.steps = 0, family = binomial) # no selection
#' metapred(DVTipd, strata = "cluster", formula = f, recal.int = TRUE, family = binomial)
#' metapred(DVTipd, strata = "cluster", formula = f, meta.method = "REML", family = binomial)
#'}
#' # By default, metapred assumes the first column is the outcome.
#' DVTipd.reordered <- DVTipd[c("dvt", "ddimdich", "histdvt", "cluster")]
#' mp <- metamisc::metapred(DVTipd.reordered, strata = "cluster", family = binomial)
#' fitted <- predict(mp, newdata = DVTipd.reordered)
#' 
#' 
#' @import stats
#'
#' @importFrom stats formula var
#'
#' @export
metapred <- function(data, strata, formula = NULL, estFUN = "glm", scope = NULL, retest = FALSE, max.steps = 1000, 
                     center = TRUE, recal.int = FALSE, cvFUN = NULL, cv.k = NULL,  # tol = 0,
                     metaFUN = NULL, meta.method = NULL, predFUN = NULL, perfFUN = NULL, genFUN = NULL, selFUN = "which.min",
                     ...) {
  call   <- match.call()
  data   <- remove.na.obs(as.data.frame(data))
  
  if (is.null(formula)) formula <- stats::formula(data[ , -which(colnames(data) == strata)])  
  if (is.null(scope)) scope <- f2iof(formula)
  updates <- getFormulaDiffAsChar(formula, scope)
  
  strata.i <- data[, strata]
  strata.u <- sort(unique(strata.i))
  # data     <- centerData(data, cluster.vec = strata.i, center.which = center) # change for 1 stage? # Deprecated
  if (center)
    data <- centerCovs(data = data, y.name = f2o(formula), cluster.name = strata)
  
  if (is.null(cvFUN))   cvFUN   <- "l1o"
  if (is.null(metaFUN)) metaFUN <- "urma"
  if (is.null(perfFUN)) perfFUN <- "mse"
  if (is.null(genFUN))  genFUN  <- "absmean"
  if (is.null(meta.method)) meta.method <- "DL"
  # Change to "-" when perfFUN <- R2 or some other measure for which greater = better.
  
  estFUN.name <- estFUN
  estFUN  <- match.fun(estFUN)
  cvFUN   <- get(cvFUN)
  perfFUN <- get(perfFUN)
  genFUN  <- get(genFUN)
  selFUN  <- get(selFUN)
  metaFUN <- get(metaFUN)
  
  folds <- cvFUN(strata.u, k = cv.k)
  if (!isTRUE(length(folds$dev) > 0) || !isTRUE(length(folds$dev[[1]]) > 0))
    stop("At least 1 cluster must be used for development.")
  
  fit <- mp.fit(formula = formula, data = data, to.change = updates, st.i = strata.i, st.u = strata.u, folds = folds,
                recal.int = recal.int, retest = retest, max.steps = max.steps, tol = 0, 
                estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, predFUN = predFUN, perfFUN = perfFUN,
                genFUN = genFUN, selFUN = selFUN, ...)
  
  predFUN <- getPredictMethod(fit$stepwise$s0$cv$unchanged, two.stage = TRUE, predFUN = predFUN)
  
  out <- c(fit, list(call = call, strata = strata, data = data,
                     options = list(cv.k = cv.k, meta.method = meta.method, recal.int = recal.int,
                                    center = center, max.steps = max.steps, retest = retest), # add: tol
                     FUN = list(cvFUN = cvFUN, predFUN = predFUN, perfFUN = perfFUN, metaFUN = metaFUN, genFUN = genFUN, # perfFUN = perfFUN
                                selFUN = selFUN, estFUN = estFUN, estFUN.name = estFUN.name))) 
  class(out) <- c("metapred")
  return(out)
}

# For prediction of newdata. Not used internally.
# Object metapred object
# newdata data.frame, defaults to NULL, which gives the fitted data
# strata character, name of strata variable in newdata. Defaults to name in fitted object.
# type character. Type of prediction. This is intended to override the default of glm and lm.
# ("response" or "link" possible; "terms" not implemented)
# recal.int logical. Recalibrate the intercept before prediction? Defaults to same as development for metapred object,
# center locigal, Center covariates, before prediction? Defaults to same as development of metapred object,
# ... For compatibility only.
#' @author Valentijn de Jong
#' @importFrom stats predict
#' @method predict   metapred
#' @export
predict.metapred <- function(object, newdata = NULL, strata = NULL, type = "response", 
                             recal.int = NULL, center = NULL, ...) {
  if (is.null(newdata))
    newdata <- object[["data"]]
  if (is.null(strata))
    strata <- object$strata
  if (is.null(recal.int))
    recal.int <- object$options$recal.int
  if (is.null(center))
    center <- object$options$center
  
  # newdata <- centerData(newdata, cluster.var = strata, center.which = center) # Deprecated
  newdata <- centerCovs(data = newdata, y.name = f2o(formula(object)), cluster.name = strata)
  
  if (isTRUE(recal.int))
    object <- recalibrate(object = object, newdata = newdata)
  
  object$FUN$predFUN(object = object, newdata = newdata, type = type, ...)
}

# #' Extract the regression coefficients
# #' 
# #' The \code{coef} function extracts the estimated coefficients of the final model of objects of class \code{"metapred"}.
# #' @param object A fitted \code{metapred} object
# #' @param \ldots Optional arguments (currently not supported).
# #' 
# #' @method coef metapred
# #' @export
#' @author Valentijn de Jong
#' @method coef   metapred
#' @importFrom stats coef
#' @export
coef.metapred <- function(object, stratified = FALSE, ...) {
  if (stratified) 
    coef(object$global.model$stratified.fit)
  else
    coef(object$global.model)
}

# Note returned formula is for the final model!
#' @author Valentijn de Jong
#' @importFrom stats formula
#' @method formula   metapred
#' @export
formula.metapred <- function(x, ...)
  x$global.model$formula

# Returns family or NULL
#' @author Valentijn de Jong
#' @importFrom stats family
#' @method family   metapred
#' @export
family.metapred <- function(object, ...) {
  if (!is.null((f <- object$stepwise$s0$cv[[1]]$family)))
    f
  else 
    NULL
}

#' @author Valentijn de Jong
#' @method print   metapred
#' @export
print.metapred <- function(x, ...) {
  cat("Call: ")                       
  print(x$call) # cat cannot handle a call
  cat("\n") 
  print.mp.fit(x)
}

#' @author Valentijn de Jong
#' @method summary   metapred
#' @export
summary.metapred <- function(object, ...) {
  cat("Call: ")                       
  print(object$call) # cat cannot handle a call
  cat("\n") 
  summary.mp.fit(object)
}

# Test whether object is meta.pred
#' @author Valentijn de Jong
#' @importFrom methods is
#' @method is   metapred
#' @export
is.metapred <- function(object)
  inherits(object, "metapred")

# Implementation of the subset method
#' @author Valentijn de Jong
#' @method subset   metapred
#' @export
subset.metapred <- function(x, select = "best", ...) {
  if (identical(select, "best")) 
    return(mp.step.get.best(x[["stepwise"]][[x[["best.step"]]]]))
  if (identical(select, "global"))
    return(x[["global.model"]])
  stop("select must equal 'best' or 'global'.")
}

### Performance / error functions
# Error function: Mean Squared Error
mse <- brier <- function(p, y, data = NULL, ...) mean((p - y)^2)


# Error function: Variance of prediction error
vare <- function(p, y, data = NULL, ...) var(p - y)

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

# Gets the predict method.
# fit Model fit object.
# two.stage logical. Is the model a two-stage model?
# predFUN Optional function, which is immediately returned
# ... For compatibility only.
getPredictMethod <- function(fit, two.stage = TRUE, predFUN = NULL, ...) {
  # A user written function may be supplied:
  if (!is.null(predFUN))
    if (is.function(predFUN))
      return(get(predFUN))
  else stop("predFUN should be a function.")
  
  # If two-stage, the fit is used only to extract the link function.
  # If one-stage, fit's prediction method may be used.
  if (two.stage) {# Preferably mp.cv.dev should not be here. But it has to.
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
# b data.frame or matrix, containing coef
# v data.frame or matrix, containing variances
# method Method for meta-analysis.
# ... Optional arguments for rma().
#' @importFrom metafor rma
urma <- function(b, v, method = "DL", ...)
{
  if (!(is.data.frame(b) || is.matrix(b)) || !(is.data.frame(v) || is.matrix(v)) )
    stop("b and v must both be a data.frame or matrix.")
  if (!identical(dim(b), dim(v)))
    stop("b and v must have the same dimensions.")
  
  meta.b <- meta.se <- rep(NA, ncol(b))
  for (col in 1:ncol(b)) {
    r <- metafor::rma(b[ , col] , v[ , col], method = method, ...)
    meta.b[col]  <- r$beta
    meta.se[col] <- r$se
  }
  
  meta.v <- meta.se^2
  
  names(meta.b) <- names(meta.v) <- names(meta.se) <- colnames(b)
  list(b = meta.b, v = meta.v, se = meta.se)
}

# Perform fit a model using cross-validation, for metapred
# formula formula to start with
# data data.frame, containing dev and val data
# to.change predictor terms to add or remove
# st.i numeric vector, cluster indicators
# st.u unique values in st.i
# folds list, folds as in utils
# recal.int logical, recalibrate intercept?

# retest logical, should removed (added) predictors be added (removed) again?
# max.steps numeric, maximum number of steps (predictor additions/removals) to be taken
# tol numeric, tolerance, minimum change in generalizability to accept a different model # To be implemented

# estFUN function, used for obtaining predict method
# predFUN function, user supplied predict method, overrides estFUN's predict()
# perfFUN function, function for computing performance, defaults to mse = mean square error
# genFUN function, function for computing generalizability, defaults to absolute mean.
# ... options for predictMethod, perfFUN, and genFUN.
# Returns object of class mp.cv.val, which is a validated mp.cv.dev
mp.fit <- function(formula, data, to.change, st.i, st.u, folds, recal.int = FALSE, 
                   retest = FALSE, max.steps = 3, tol = 0,
                   estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                   perfFUN = mse, genFUN = absmean, selFUN = which.min, ...) {
  out <- steps <- list()
  
  ## Step 0
  # As to.change = c("") yields unchanged formula.
  step.count <- 0
  steps[[getStepName(step.count)]] <- mp.step(formula = formula, data = data, to.change = c(""), 
                                     st.i = st.i, st.u = st.u, folds = folds, recal.int = recal.int, 
                                     retest = FALSE,
                                     estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                                     perfFUN = mse, genFUN = absmean, selFUN = selFUN, ...)
  steps[[getStepName(step.count)]][["step.count"]] <- step.count
  out[["best.step"]] <- getStepName(step.count)
  current.model <- mp.step.get.best(steps[[1]])
  current.model[["to.change"]] <- to.change
  gen.diff <- Inf # TBI

  
  if (!identical(length(to.change), 0L))
    repeat {
      ## Loop management
      if (isTRUE(gen.diff <= tol)) { # TBI
        out[["stop.reason"]] <- "improvement <= tolerance."
        break
      }
      if (isTRUE(length(current.model[["to.change"]]) <= 0)) {
        out[["stop.reason"]] <- "all changes were implemented."
        break
      }
      if (isTRUE(step.count >= max.steps)) {
        out[["stop.reason"]] <- "max.steps was reached."
        break
      }
      step.count <- step.count + 1
      
      ## Perform a step
      steps[[getStepName(step.count)]] <- mp.step(formula = current.model[["formula"]], data = data,
                                                  to.change = current.model[["to.change"]],
                                                  st.i = st.i, st.u = st.u, folds = folds, recal.int = recal.int,
                                                  retest = retest,
                                                  estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL,
                                                  perfFUN = mse, genFUN = absmean, selFUN = selFUN, ...)
      steps[[getStepName(step.count)]][["step.count"]] <- step.count
      ## Model selection
      # This step
      best.new.model <- mp.step.get.best(steps[[getStepName(step.count)]])
      
      # Overall
      if (mp.which.best.change(list(current.model, best.new.model), selFUN = selFUN) <= 1) { # TBI: 0 for gen.diff < tol
        out[["stop.reason"]] <- "no improvement was possible."
        out[["best.step"]] <- getStepName(step.count - 1)
        break
      } else {
        current.model <- best.new.model
        out[["best.step"]] <- getStepName(step.count)
      }
    }
  
  # Return a global model and the stepwise list
  out[["global.model"]] <- mp.global(current.model, metaFUN = metaFUN, meta.method = meta.method)
  out[["stepwise"]] <- steps
  out[["step.count"]] <- step.count
  class(out) <- "mp.fit"
  
  out
}

#' @author Valentijn de Jong
#' @method print   mp.fit
#' @export
print.mp.fit <- function(x, ...) {
  for (i in seq_along(x[["stepwise"]])) {
    if (i == 1) {
      cat("Started with model:\n")
      print(x$stepwise[[i]][["start.formula"]])
    } else if (i > 2) {
      cat("\n")
      cat("Continued with model:\n")
      print(x$stepwise[[i]][["start.formula"]])
    }
    print.mp.step(x$stepwise[[i]], show.f = FALSE)
  }
  
  cat("\n")
  cat("Cross-validation stopped after", x[["step.count"]], "steps, as", x[["stop.reason"]])
  cat(" Final model:\n")
  print(x[["global.model"]])
}

#' @author Valentijn de Jong
#' @method summary   mp.fit
#' @export
summary.mp.fit <- function(object, ...) {
  for (s in seq_along(object[["stepwise"]]))
    summary.mp.step(object[["stepwise"]][[s]])
  cat("\n")
  cat("Cross-validation stopped after", object[["step.count"]], "steps, as", object[["stop.reason"]])
  cat(" Final model:\n")
  print(object[["global.model"]])
}

# "Method" for selecting best mp.cv from a mp.cv object
# step mp.step
# Returns mp.cv object, with best value for genFUN.
# mp.step.best.change <- function(step, ...)
#   step[["cv"]][[which(sapply(step[["cv"]], `[[`, "best.change"))]]

# Select best model
mp.which.best.change <- function(cvs, selFUN = which.min, ...)
  selFUN(sapply(cvs, `[[`, "gen"))

mp.step.get.best <- function(step, selFUN = which.min, ...)
  step[["cv"]][[mp.which.best.change(step[["cv"]], selFUN = selFUN, ...)]]

# Perform one step in the fitting process of mp.fit
# formula formula to start with
# data data.frame, containing dev and val data
# to.change predictor terms to add or remove
# st.i numeric vector, cluster indicators
# st.u unique values in st.i
# folds list, folds as in utils
# recal.int logical, recalibrate intercept?
# two.stage logical, is it a two stage model? For future use.
# estFUN function, used for obtaining predict method
# predFUN function, user supplied predict method, overrides estFUN
# perfFUN function, function for computing performance, defaults to mse = mean square error
# genFUN function, function for computing generalizability, defaults to absolute mean.
# ... options for predictMethod, perfFUN, and genFUN.
# Returns object of class mp.cv.val, which is a validated mp.cv.dev
mp.step <- function(formula, data, to.change, st.i, st.u, folds, recal.int = FALSE, 
                    two.stage = TRUE, retest = FALSE,
                    estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                    perfFUN = mse, genFUN = absmean, selFUN = which.min, ...) {
  cv <- out <- list()
  out[["start.formula"]] <- formula
  
  for (fc in seq_along(to.change) )
  {
    change <- to.change[fc]
    
    # Produce formula for changes and no changes:
    if (identical(to.change, "")) {
      name <- "unchanged"
      new.formula <- formula
    } else {
      name <- as.character(change)
      new.formula <- updateFormula(formula, change) 
    }
    
    # Run
    cv[[name]] <- mp.cv(formula = new.formula, data = data, st.i = st.i, st.u = st.u,
                        folds = folds, recal.int = recal.int,
                        estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method,
                        predFUN = predFUN, perfFUN = perfFUN, genFUN = genFUN, ...)
    # Save changes
    cv[[name]][["to.change"]]   <- if (retest) to.change else to.change[-fc]
    cv[[name]][["changed"]]     <- change
  }
  
  out[["best.change"]] <- mp.which.best.change(cv, selFUN = selFUN)
  out[["cv"]] <- cv
  
  class(out) <- "mp.step"
  out
}

#' @author Valentijn de Jong
#' @method summary   mp.step
#' @export
summary.mp.step <- function(object, ...) {
  cv <- object[["cv"]]
  
  # Model comparison is handled by print.mp.step()
  if (!is.null(object[["step.count"]]))
    cat("\n", "Step ", object[["step.count"]], ". ", sep = "")
  if (identical(cv[[1]][["changed"]], "")) {
    # cat("\n")
    cat("Tested 1 model in this step. ") 
  } else {
    # cat("\n")
    cat("Tested", length(cv), "models in this step. ")
    print.mp.step(object)
    cat("\n")
  }
  
  # Model coefficients are handled here, using print.mp.cv()
  for (i in seq_along(cv)) {
    ch <- cv[[i]][["changed"]]
    if (identical(ch, "")) ch <- "nothing"
    cat("Changing", ch, "yields:\n")
    print.mp.cv(cv[[i]])
    if (isTRUE(i == object[["best.change"]][[1]]))
      cat("This model has best generalizability in this step.")
    cat("\n")
  }
}

#' @author Valentijn de Jong
#' @method print   mp.step
#' @export
print.mp.step <- function(x, show.f = TRUE, ...) {
  if (show.f) {
    cat("Starting with model:\n")
    print(x[["start.formula"]])
  }
  cat("\n")
  Change.to.model <- names(x[["cv"]])
  Generalizability <- sapply(x[["cv"]], `[[`, "gen")
  Mean.performance <- sapply(x[["cv"]], `[[`, "mean.perf")
  print(data.frame(Change.to.model, Generalizability, Mean.performance), row.names = FALSE)
}

# Turn a cross-validated model into a full or 'global' model
# mp.cv.dev object of class mp.cv.dev
# metaFUN function for meta-analysis
# meta.method option for metaFUN
# ... optional arguments for metaFUN
mp.global <- function(cv.dev, metaFUN = urma, meta.method = "DL") {
  out <- c(cv.dev, mp.meta.fit(cv.dev[["stratified.fit"]], metaFUN, meta.method) )
  class(out) <- c("mp.global", class(out))
  out
}

#' @author Valentijn de Jong
#' @method print   mp.global
#' @export
print.mp.global <- function(x, ...) {
  cat("Meta-analytic model of prediction models estimated in", x$n.clusters, "strata. Coefficients:" , "\n")
  print(x[["coefficients"]])
}

#' @author Valentijn de Jong
#' @method summary   mp.global
#' @export
summary.mp.global <- function(object, ...) {
  print.mp.global(object, ...)
  cat("\n")
  print.mp.cv(object, ...)
}

# Make a new meta-model and cross-validate it, for metapred
# formula formula
# data data.frame
# st.i vector, cluster indicators. char or numeric
# st.u unique values in st.i
# recal.int logical, recalibrate intercept?
# estFUN function, for estimating, e.g. glm
# metaFUN function, for producing meta model
# meta.method character, option for metaFUN
# predFUN function, user supplied predict method, overrides estFUN
# perfFUN function, function for computing performance, defaults to mse = mean square error
# genFUN function, function for computing generalizability, defaults to absolute mean.
# ... options for predictMethod, perfFUN, and genFUN.
# Returns object of classes mp.cv, mp.cv.val, mp.cv.dev, which is a list of meta-analytic models developed on dev folds,
# and a validated on val folds 
mp.cv <- function(formula, data, st.i, st.u, folds, recal.int = FALSE, two.stage = TRUE,
                  estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                  perfFUN = mse, genFUN = absmean, ...) {
  
  out <- mp.cv.dev(formula = formula, data = data, st.i = st.i, st.u = st.u, folds = folds, two.stage = two.stage,
                   estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, ...)
  
  out <- mp.cv.val(cv.dev = out, data = data, st.i = st.i, folds = folds, recal.int = recal.int, two.stage = two.stage,
                   estFUN = glm, predFUN = predFUN, perfFUN = mse, genFUN = absmean, ...)
  
  class(out) <- c("mp.cv", class(out))
  out
}

#' @author Valentijn de Jong
#' @method print   mp.cv
#' @export
print.mp.cv <- function(x, ...) {
  print.mp.cv.dev(x, ...)
  cat("\n")
  print.mp.cv.val(x, ...)
}

# Validate a mp.cv.dev
# cv.dev, mp.cv.dev
# data data.frame, containing old and new data
# st.i numeric vector, cluster indicators
# folds list, folds as in utils
# recal.int logical, recalibrate intercept?
# two.stage logical, is it a two stage model? For future use.
# estFUN function, used for obtaining predict method
# predFUN function, user supplied predict method, overrides estFUN
# perfFUN function, function for calculating performance, defaults to mse = mean square error
# genFUN function, function for calculating generalizability, defaults to absolute mean.
# add.perfFUN list of functions, additional performance functions. Evaluated but not used for selection.
# add.genFUn list of functions, additional generalizability functions. Evaluated but not used for selection.
# ... options for predictMethod, perfFUN, and genFUN.
# Returns object of class mp.cv.val, which is a validated mp.cv.dev
mp.cv.val <- function(cv.dev, data, st.i, folds, recal.int = FALSE, two.stage = TRUE,
                      estFUN = glm, predFUN = NULL, perfFUN = mse, add.perfFUN = list(), genFUN = absmean, ...) { # add: add.genFUN
  # Recalibrate?
  if (isTRUE(recal.int))
    cv.dev <- mp.cv.recal(cv.dev = cv.dev, newdata = data, estFUN = estFUN, folds = folds)
  
  # Predict outcome 
  predictMethod <- getPredictMethod(fit = cv.dev, two.stage = two.stage, predFUN = predFUN, ...)
  p <- list()
  for (i in seq_len(cv.dev[["n.cv"]]))
    p[[i]] <- predictMethod(object = cv.dev, newdata = data[folds[["val"]][[i]] == st.i, ], type = "response",
                            b = coef(cv.dev[["cv"]][[i]]),
                            f = formula(cv.dev), two.stage = two.stage, ...)
  
  # Necessary for performance computation
  outcome <- f2o(formula(cv.dev))
  
  cv.dev[["perf"]] <- data.frame(matrix(ncol = 2, nrow = cv.dev[["n.cv"]]))
  row.names(cv.dev[["perf"]]) <- names(cv.dev[["cv"]])
  colnames(cv.dev[["perf"]])  <- c("val.strata", "perf")
  
  # Compute performance for predictor selection
  # perFUN receives args that might be useful, but not necessary for defaults
  for (i in seq_len(cv.dev[["n.cv"]])) 
    cv.dev[["perf"]][i, ] <- c(getclName(folds[["val"]][[i]]), 
                               perfFUN(p[[i]], data[folds[["val"]][[i]] == st.i, outcome],
                                       data = data, fit = cv.dev[["cv"]][[i]], ...)) 
  cv.dev[["perf"]][, "perf"] <- as.numeric(cv.dev[["perf"]][, "perf"])
  
  # Compute additional performance measures.
  # ...
  
  # And finally, the generalizability (and mean performance)
  cv.dev[["gen"]]       <-  genFUN(cv.dev[["perf"]][, "perf"], ...)
  cv.dev[["mean.perf"]] <- absmean(cv.dev[["perf"]][, "perf"], ...)
  
  class(cv.dev) <- c("mp.cv.val", class(cv.dev))
  cv.dev
}

#' @author Valentijn de Jong
#' @method print   mp.cv.val
#' @export
print.mp.cv.val <- function(x, ...) {
  if (!is.null(x[["perf"]])) {
    cat("Cross-validation at stratum level yields the following performance: \n")
    print(x[["perf"]])
  }
  if (!is.null(x[["gen"]])) {
    cat("\n")
    cat("Generalizability value:",  x[["gen"]], "\n")
  }
}

# Make a new meta-model to be cross-validated for metapred
# formula formula
# data data.frame
# estFUN function for estimation, e.g glm
# st.i vector, cluster indicators. char or numeric
# st.u unique values in st.i
# folds list, folds as in utils
# two.stage logical, is it a two stage method?
# estFUN function, for estimating, e.g. glm
# metaFUN function, for producing meta model
# meta.method character, option for metaFUN
# Returns mp.cv.dev, which is a list of meta-analytic models developed on dev folds
mp.cv.dev <- function(formula, data, st.i, st.u, folds, two.stage = TRUE, 
                      estFUN = glm, metaFUN = urma, meta.method = "DL", ...) {
  out <- list(...) # possibly contains family
  if (!is.null(out$family)) {
    if (is.character(out$family)) # Ensure that family is the output of family(). Taken directly from glm.
      out$family <- get(out$family, mode = "function", envir = parent.frame())
    if (is.function(out$family)) 
      out$family <- out$family()
    if (is.null(out$family$family)) {
      print(family)
      stop("'family' not recognized")
    }
  }
  out[["formula"]] <- formula
  
  ## Cluster part 
  out[["stratified.fit"]] <- mp.stratified.fit(formula = formula, data = data, st.i = st.i, st.u = st.u, estFUN = estFUN, ...)  
  out[["n.clusters"]] <- length(out[["stratified.fit"]])
  
  ## cv meta-part
  out[["cv"]] <- mp.cv.meta.fit(stratified.fit = out[["stratified.fit"]], folds = folds, metaFUN = metaFUN, meta.method = meta.method)
  
  out[["n.cv"]] <- length(out$cv)
  class(out) <- "mp.cv.dev"
  out
}

# For some reason I have to manually force it to find the right S3 methods.
#' @author Valentijn de Jong
#' @method print   mp.cv.dev
#' @export
print.mp.cv.dev <- function(x, ...) {
  print.mp.stratified.fit(x[["stratified.fit"]]) 
  cat("\n")
  print.mp.cv.meta.fit(x[["cv"]])
}

#' @author Valentijn de Jong
#' @method family   mp.cv.dev
#' @export
family.mp.cv.dev <- function(object, ...) 
  object$family

# Recalibrate a mp.meta.fit
# mp.meta.fit mp.meta.fit
# newdata # ne data.frame, containing only the validation sample (unlike this example!)
# formula, original formula
# estFUN estimation function, e.g. glm
# Returns same object with updated coefficients.
mp.recal.meta.fit <- function(meta.fit, formula, newdata, estFUN, ...) {
  meta.fit[["formula"]] <- formula
  meta.fit[["orig coef"]] <- coef(meta.fit)
  meta.fit[["coefficients"]] <- recalibrate(meta.fit, newdata = newdata, estFUN = estFUN, ...)[["coefficients"]]
  meta.fit[["formula"]] <- NULL ### necessary ??
  meta.fit
}

# Recalibrate a mp.cv.dev
# mp.cv.dev mp.cv.dev
# newdata data.frame
# formula, original formula
# estFUN estimation function, e.g. glm
# folds folds list, as in utils.
# Returns same object with updated coefficients.

mp.cv.recal <- function(cv.dev, newdata, folds, estFUN) {
  for (i in seq_along(cv.dev[["cv"]]))
    cv.dev[["cv"]][[i]] <- mp.recal.meta.fit(meta.fit = cv.dev[["cv"]][[i]], 
                                             formula = cv.dev[["formula"]],
                                             newdata = newdata[folds[["val.i"]][[i]] == i, ], 
                                             estFUN = estFUN, family = if (!is.null(cv.dev$family)) cv.dev$family else NULL)
  cv.dev
}

# Make new cv.models for mp.meta.fit for metapred
# This function selects the right fold for development, and calls the fitting function.
# stratified.fit list of mp.stratum.fit objects
# folds list of fold divisions, as given by l1o, or bootstrap in the utils.
# metaFUN function for estimating meta-analytic models, e.g. urma (this file)
# meta.method options for metaFUN
# Returns a list of objects of class mp.cv.meta.fit, which are meta-analytic prediction models

mp.cv.meta.fit <- function(stratified.fit, folds, metaFUN = urma, meta.method = "DL") {
  out <- list()
  
  for (fo in seq_along(folds[["dev.i"]]))
    out[[getclName(folds[["dev"]][[fo]])]] <- mp.meta.fit(stratified.fit = stratified.fit[folds[["dev.i"]][[fo]]], 
                                                          metaFUN = metaFUN, meta.method = meta.method)
  class(out) <- "mp.cv.meta.fit"
  out
}

#' @author Valentijn de Jong
#' @method coef   mp.cv.meta.fit
#' @export
coef.mp.cv.meta.fit <- function(object, ...) 
  t(as.data.frame(lapply(object, `[[`, "coefficients"), check.names = FALSE))

#' @author Valentijn de Jong
#' @method print   mp.cv.meta.fit
#' @export
print.mp.cv.meta.fit <- function(x, ...) {
  cat("Meta-analytic models, estimated in", length(x), "fold combinations. Coefficients: \n")
  print(coef(x))
  
  if (!is.null(x[[1]][["orig coef"]])) {
    cat("\n")
    cat("Original coefficients before recalibration: \n")
    print(t(as.data.frame(lapply(x, `[[`, "orig coef"), check.names = FALSE) ))
  }
}

# Make new meta model (i.e. model fitted on multiple clusters) for ?? for metapred
# stratified.fit mp.stratified.fit
# metaFUN function for estimating meta-analytic models, e.g. urma (this file)
# meta.method options for metaFUN
# Returns object of class mp.cv.model, which is a meta-analytic prediction model
mp.meta.fit <- function(stratified.fit, metaFUN = urma, meta.method = "DL") {
  out <- list()
  
  b <- coef.mp.stratified.fit(stratified.fit) # again i have to point it to the right method.
  v <- t(as.data.frame(lapply(stratified.fit, `[[`, "v", drop = FALSE)))
  
  meta <- metaFUN(b = b, v = v, method = meta.method) 
  
  out[["coefficients"]] <- meta$b
  out[["v"]] <- meta$v
  
  class(out) <- "mp.meta.fit"
  out
}

#' @author Valentijn de Jong
#' @method print   mp.meta.fit
#' @export
print.mp.meta.fit <- function(x, ...) {
  cat("Coefficients: ", "\n")
  print(coef(x))
}

# Make new stratified models for mp.meta.fit for metapred
# formula formula
# data data.frame
# st.i vector, strata indicators. char or numeric
# st.u unique values in st.i
# estFUN function for estimation, e.g glm
# ... other, e.g. family for glm.
# Returns mp.stratified.fit
mp.stratified.fit <- function(formula, data, st.i, st.u, estFUN, ...) {
  out <- list()
  
  for (st.this in st.u)
    out[[getclName(st.this)]] <- mp.stratum.fit(estFUN(formula = formula, data = data[st.i == st.this, ], ...) )
  
  class(out) <- "mp.stratified.fit"
  out
}

#' @author Valentijn de Jong
#' @method coef   mp.stratified.fit
#' @export
coef.mp.stratified.fit <- function(object, ...)
  t(as.data.frame(lapply(object, `[[`, "coefficients", drop = FALSE)))

#' @author Valentijn de Jong
#' @method print   mp.stratified.fit
#' @export
print.mp.stratified.fit <- function(x, ...) {
  cat("Prediction models estimated in", length(x), "strata. Coefficients:\n")
  print(coef(x))
}

# Make a new mp.stratum.fit for mp.stratified.fit for metapred.
# fit model fit object, e.g. glm object
# Returns mp.stratum.fit
mp.stratum.fit <- function(fit) {
  out <- list()
  
  out[["coefficients"]] <- getCoefs(fit)
  out[["v"]]            <- getVars(fit)
  out[["covar"]]        <- getCoVars(fit)
  out[["vcov"]]         <- vcov(fit)
  
  class(out) <- "mp.stratum.fit"
  out
}

#' @author Valentijn de Jong
#' @method print   mp.stratum.fit
#' @export
print.mp.stratum.fit <- function(x, ...) {
  cat("Coefficients: ", "\n")
  print(coef(x))
}