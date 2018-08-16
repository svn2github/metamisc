# TODO 1: Thomas' ideas
# @Valentijn: perhaps better to combine param perfFUN, genFUN and selFUN into one parameter with 3 or 4 distinct options
# @Valentijn: I suggest to omit intercept recalibration. The intercept issue can be addressed directly by specifying distinct 
# error functions that do or do not account for mis-calibration in intercept term

# TODO 2: metapred objects
# add more output for metapred
# tol, aka tolerance for stopping.
# metapred.summaries docs (see glm.summaries {stats}), including subset
# add stepnumbers to print.metapred or print.fit
# Write tests for fitted.whatever
# add docs for fitted.whatever.
# check what happens for recalibrated mp.cv in fitted.mp.cv
# add rownames to data, to prevent error/mismatch in fitted.metapred/fitted.mp.cv # Maybe already fixed by as.data.frame()
# add the unchanged model to the next step's model list. Model comparison will be easier. Subsetting will be easier: this allows
#   a 'chain' of "changed" predictors to be added, thereby making the output clearer. 
# change "changed" for global models as well.

# TODO 3: General
# add function: areTRUE. st.i == cl may otherwise lead to problems, if st.i has NA
# variances for intercept recalibration.
# perf method
# penalties
# performance measurement, current is Brier for binomial
# add more options to penalty step
# Write test for new centerCovs function.

# TODO 4: Maybe later
# One-stage MA: predStep: is  type = response correct?
# centering: now only centers numeric variables. Should it also work for dummies of categorical variables?

# TODO 5: DONE:
# recal.int + interaction gives error. FIXED.
# add fitted.metapred(). DONE.
# categorical variables Done!
# Add summary.metapred. DONE.
# Added is.metapred()
# Added more robust tests for stepwise.
# Automatically remove observations with missing values.
# @Valentijn: I would center covariates by default, this often helps to improve generalizability and also speeds estimation
#   Valentijn: Done.
# @Valentijn: You can change default of meta.method to DL. This is a lot faster and has limited implications on estimated means.
#   Valentijn: Done.
# , cl.name in modelStep
# Remove Reduce() from perfStep, and make perfStep() compatible with multiple data sets (for cvFUN = fixed, or cvFUN = bs)

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
#' @param data data.frame containing the data. Note that \code{metapred} removes observations with missing data \emph{listwise},
#' to ensure that the same data is used in each model.
#' @param strata Name of the strata (e.g. studies or clusters) variable, as character.
#' @param formula \code{formula} of the first model to be evaluated. \code{metapred} will start at \code{formula} and update it
#' using terms of \code{scope}. Defaults to full main effects model, where the first column in \code{data} is assumed to be
#' the outcome and all remaining columns (except \code{strata}) predictors. See \link[stats]{formula} for formulas in general.
#' @param estFUN Function for estimating the model in the first stage. Currently "lm" and "glm" are supported.
#' @param scope \code{formula}. The difference between \code{formula} and \code{scope} defines the range of models examined in the 
#' stepwise search. Defaults to NULL, which leads to the intercept-only model. If \code{scope} is nested in \code{formula}, 
#' this implies backwards selection will be applied (default). If \code{scope} is nested in \code{formula}, this implies forward 
#' selection will be applied. If equal, no stepwise selection is applied. 
#' @param retest Logical. Should added (removed) terms be retested for removal (addition)? \code{TRUE} implies bi-directional 
#' stepwise search.
#' @param max.steps Integer. Maximum number of steps (additions or removals of terms) to take. Defaults to 1000, which is
#' essentially as many as it takes. 0 implies no stepwise selection.
#' @param center logical. Should numeric predictors be centered?
#' @param recal.int Logical. Should the intercept be recalibrated in each validation?
#' @param cvFUN Cross-validation method, on the study (i.e. cluster or stratum) level. "
#' l1o" for leave-one-out cross-validation (default). "bootstrap" for bootstrap. Or "fixed", for one or more data sets
#' which are only used for validation. A user written function may be supplied as well.
#' @param cv.k Parameter for cvFUN. For \code{cvFUN="bootstrap"}, this is the number of bootstraps. For \code{cvFUN="fixed"}, 
#' this is a vector of the indices of the (sorted) data sets. Not used for \code{cvFUN="l1o"}.
#' @param metaFUN Function for computing the meta-analytic coefficient estimates in two-stage MA. 
#' By default, \link[metafor]{rma.uni}, from the metafor package is used. Default settings are univariate random effects, 
#' estimated with "DL". Method can be passed trough the \code{meta.method} argument.
#' @param meta.method Name of method for meta-analysis. Default is "DL". For more options see \link[metafor]{rma.uni}.
#' @param predFUN Function for predicting new values. Defaults to the predicted probability of the outcome, using the link 
#' function of \code{glm()} or \code{lm()}.
#' @param perfFUN Function for computing the performance of the prediction models. Default: mean squared error 
#' (\code{perfFUN="mse"}).Other options are \code{"vare"}.
#' @param genFUN Function computing generalizability measure using the performance measures. Default: (absolute) mean  
#' (\code{genFUN="abs.mean"}). Choose \code{squared.diff} for a penalty equal to the mean squared differences between 
#' coefficients. Alternatively, choose \code{pooled.var} for a weighted average of variance terms.
#' @param selFUN Function for selecting the best method. Default: lowest value for \code{genFUN}. Should be set to
#' "which.max" if high values for \code{genFUN} indicate a good model.
#' @param ... To pass arguments to estFUN (e.g. family = "binomial"), or to other FUNctions.
#'
#' @return A list of class \code{metapred}, containing the final model in \code{global.model}, and the stepwise
#' tree of estimates of the coefficients, performance measures, generalizability measures in \code{stepwise}.
#' 
#' @details Use \link{subset} to obtain an individual prediction model from a \code{metapred} object.
#' 
#'  Note that \code{formula.changes} is currently unordered; it does not represent the order of changes in the stepwise 
#'  procedure.
#'  
#'  \code{metapred} is still under development, use with care.
#' 
#' @examples 
#' data(DVTipd)
#' DVTipd$cluster <- letters[1:4] # Add a fictional clustering to the data.
#' #f <- dvt ~ sex + vein + malign
#' f <- dvt ~ histdvt + ddimdich
#' # Leaving scope empty implies backwards selection
#' mp <- metapred(DVTipd, strata = "cluster", formula = f, family = binomial)
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
#' mp <- metapred(DVTipd.reordered, strata = "cluster", family = binomial)
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
                     metaFUN = NULL, meta.method = NULL, predFUN = NULL, perfFUN = NULL, genFUN = NULL,
                     selFUN = "which.min",
                     ...) {
  call <- match.call()
  dots <- list(...)
  data <- remove.na.obs(as.data.frame(data))
  
  if (is.null(formula)) formula <- stats::formula(data[ , -which(colnames(data) == strata)])  
  if (is.null(scope)) scope <- f2iof(formula)
  updates <- getFormulaDiffAsChar(formula, scope)
  
  strata.i <- as.vector(data[, strata])
  strata.u <- sort(unique(strata.i))
  if (center)
    data <- centerCovs(data = data, y.name = f2o(formula), cluster.name = strata)
  
  if (is.null(cvFUN))   cvFUN   <- "l1o"
  if (is.null(metaFUN)) metaFUN <- "urma"
  if (is.null(perfFUN)) perfFUN <- "mse"
  if (is.null(genFUN))  genFUN  <- "abs.mean"
  if (is.null(meta.method)) meta.method <- "DL"
  # Change to "-" when perfFUN <- R2 or some other measure for which greater = better.
  
  estFUN.name <- estFUN
  estFUN  <- match.fun(estFUN)
  cvFUN   <- get(cvFUN)
  perfFUN <- get(perfFUN)
  genFUN  <- get(genFUN)
  selFUN  <- get(selFUN)
  metaFUN <- get(metaFUN)
  
  genFUN.add <- dots[["genFUN.add"]] 
  dots[["genFUN.add"]] <- NULL
  
  folds <- cvFUN(strata.u, k = cv.k)
  if (!isTRUE(length(folds[["dev"]]) > 0) || !isTRUE(length(folds[["dev"]][[1]]) > 0))
    stop("At least 1 cluster must be used for development.")
  
  
  
  # fit <- mp.fit(formula = formula, data = data, remaining.changes = updates, st.i = strata.i, st.u = strata.u, folds = folds,
  #               recal.int = recal.int, retest = retest, max.steps = max.steps, tol = 0, 
  #               estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, predFUN = predFUN, perfFUN = perfFUN,
  #               genFUN = genFUN, genFUN.add = genFUN.add, selFUN = selFUN, ... = unlist(dots)) # genFUN.add = genFUN.add,
  mp.dots <<- dots
  mp.genFUN.add <<- genFUN.add
  mp.args <- c(list(formula = formula, data = data, remaining.changes = updates, st.i = strata.i, st.u = strata.u, folds = folds,
                    recal.int = recal.int, retest = retest, max.steps = max.steps, tol = 0, 
                    estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, predFUN = predFUN, perfFUN = perfFUN,
                    genFUN = genFUN, genFUN.add = genFUN.add, selFUN = selFUN), dots)
  
  mp.args <<- mp.args
  fit <- do.call(mp.fit, args = mp.args ) 
  
  predFUN <- getPredictMethod(fit$stepwise$s0$cv$unchanged, two.stage = TRUE, predFUN = predFUN)
  formula.final <- fit$global.model$formula
  
  out <- c(fit, list(call = call, strata = strata, data = data, folds = folds, # add nobs and strata.nobs
                     formula.start = formula, scope = scope, formula = formula.final,
                     formula.changes = getFormulaDiffAsChar(formula.final, formula), 
                     # NOTE: formula.changes is currently unordered!
                     options = list(cv.k = cv.k, meta.method = meta.method, recal.int = recal.int,
                                    center = center, max.steps = max.steps, retest = retest), # add: tol
                     FUN = list(cvFUN = cvFUN, predFUN = predFUN, perfFUN = perfFUN, metaFUN = metaFUN, genFUN = genFUN, 
                                selFUN = selFUN, estFUN = estFUN, estFUN.name = estFUN.name, genFUN.add = genFUN.add)))
  class(out) <- c("metapred")
  return(out)
}

# #' The function \code{summary} can be used to obtain an extensive summary of the stepwise process and the final model. 
# #' \code{subset(x, select = "best.cv")} can be used to obtain a summary of the model with best generalizability in the
# #' cross-validation, whereas \code{subset(x, select = "global")} gives a summary of the global model (i.e. fitted on all 
# #' strata).

# For prediction of newdata. 
# Object metapred object
# newdata data.frame, defaults to NULL, which implies the fitted data
# strata character, name of strata variable in newdata. Defaults to name in fitted object.
# type character. Type of prediction. This is intended to override the default of glm and lm.
# ("response" or "link" possible; "terms" not implemented)
# recal.int logical. Recalibrate the intercept before prediction? Defaults to same as development for metapred object,
# center logical, Center covariates, before prediction? Defaults to same as development of metapred object,
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
  if (is.null(center))
    center <- object$options$center
  if (center)
    newdata <- centerCovs(data = newdata, y.name = f2o(formula(object)), cluster.name = strata)
  if (is.null(recal.int))
    recal.int <- object$options$recal.int
  if (isTRUE(recal.int))
    object <- recalibrate(object = object, newdata = newdata)
  
  object$FUN$predFUN(object = object, newdata = newdata, type = type, ...)[, 1] # [, 1] such that a vector is returned.
}

#' Extract Model Fitted Values
#' 
#' Extract the fitted values of a \code{metapred} object. By default returns fitted values of the model in the 
#' cross-validation procedure.
#' 
#' Function still under development, use with caution.
#' 
#' Only returns type = "response".
#' 
#' @author Valentijn de Jong
#' @importFrom stats fitted
#' @method fitted   metapred
#' @param object object of class metapred 
#' @param select character. Select fitted values from "cv" (default) or from "global" model.
#' @param step character or numeric. Name or number of step to select if \code{select} = "cv". Defaults to best step.
#' @param model character or numeric. Name or number of model to select if \code{select} = "cv". Defaults to
#' best model.
#' @param as.stratified logical. \code{select} = "cv" determines whether returned predictions are stratified in a list 
#' (\code{TRUE}, default) or in their original order (\code{FALSE}).
#' @param type character. Type of fitted value.
#' @param ... For compatibility only.
#' @export
fitted.metapred <- function(object, select = "cv", step = NULL, model = NULL, 
                            as.stratified = TRUE, type = "response", ...) {
  if (isTRUE(select == "cv")) {
    ftd <- fitted(subset.metapred(x = object, select = select, step = step, model = model, type = type, ...))
    if (as.stratified)
      return(ftd)
    else {
      ftd.v <- Reduce(rbind, ftd) #as vector
      return(ftd.v[match(rownames(ftd.v), rownames(object[["data"]])) ] ) # return to original ordering.
    }
  }

  if (isTRUE(select == "global"))
    return(predict.metapred(object = object, newdata = NULL, type = type, ...))
  stop("select must equal 'cv' or 'global'.")
}

# #' Only returns type = "response".
# #' @author Valentijn de Jong
# #' @importFrom stats residuals
# #' @method residuals   metapred
# #' @export
# residuals.metapred <- function(object, select = "cv", step = NULL, model = NULL, as.stratified = TRUE, ...) {
#   y <- object$data[ , f2o(formula(object))]
#   ftd <- fitted.metapred(object = object, select = select, step = step, 
#                        model = model, as.stratified = as.stratified, ...)
#   if (as.stratified) {
#     ftd <- fitted(mp, as.stratified = TRUE)
#     y <- mp$data[ , f2o(formula(mp))]
#     
#     st.i <- mp[["data"]][[mp[["strata"]]]]
#     # now somehow sort them into a list with similar dim and dimnames as ftd, using st.i
#     
#     stop("To be implemented")
#   } else {
#     y - ftd
#   }
# }

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
# #' @author Valentijn de Jong
# #' @method subset   metapred
# #' @export
# old.subset.metapred <- function(x, select = "best.cv", ...) {
#   if (identical(select, "best.cv")) 
#     return(mp.step.get.best(x[["stepwise"]][[x[["best.step"]]]]))
#   if (identical(select, "global"))
#     return(x[["global.model"]])
#   stop("select must equal 'best.cv' or 'global'.")
# }

#' Subsetting metapred fits
#' 
#' Return a model from the cross-validation procedure or the final 'global' model. Caution: This function is 
#' still under development.
#' 
#' @author Valentijn de Jong
#' 
#' @param x metapred object
#' @param select Which type of model to select: "cv" (default) or "global"
#' @param step  Which step should be selected? Defaults to the best step. 
#' numeric is converted to name of the step: 0 for an unchanged model, 1 for the first change... 
#' @param model Which model change should be selected? NULL (default, best change) or character name of variable
#' or (integer) index of model change. 
#' @param ... For compatibility only.
#' 
#' @examples 
#' data(DVTipd)
#' DVTipd$cluster <- letters[1:4] # Add a fictional clustering to the data.
#' mp <- metapred(DVTipd, strata = "cluster", formula = dvt ~ histdvt + ddimdich, family = binomial)
#' subset(mp) # best cross-validated model
#' subset(mp, select = "global") # Final model fitted on all strata.
#' subset(mp, step = 1) # The best model of step 1
#' subset(mp, step = 1, model = "histdvt") # The model in which histdvt was removed, in step 1.
#' 
#' @return An object of class \code{mp.cv} for select = "cv" and an object of class \code{mp.global} for select = "global". 
#' In both cases, additional data is added to the resulting object, thereby making it suitable for further methods.
#' @export
subset.metapred <- function(x, select = "cv", step = NULL, model = NULL, ...) {
  # Select
  if (identical(select, "global"))
    out <- x[["global.model"]]
  else {
    if (!identical(select, "cv"))
      stop("select must equal 'cv' or 'global'.")
    
    # Step
    if (is.null(step))
      step <- x$best.step
    if (is.numeric(step))
      step <- getStepName(step)
    
    # Model
    if (is.null(model))
      model <- x[["stepwise"]][[step]][["best.change"]]
    
    out <- x[["stepwise"]][[step]][["cv"]][[model]]
  }
  
  # This is the real reason for why this function exists. To add this stuff to a mp.cv or mp.global object.
  # Normally those do not have this data, for memory/performance reasons.
  out$data    <- x$data
  out$strata  <- x$strata
  out$folds   <- x$folds
  out$FUN     <- x$FUN
  out$options <- x$options

  out
}

# Perform fit a model using cross-validation, for metapred
# formula formula to start with
# data data.frame, containing dev and val data
# remaining.changes predictor terms to add or remove
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
mp.fit <- function(formula, data, remaining.changes, st.i, st.u, folds, recal.int = FALSE, 
                   retest = FALSE, max.steps = 3, tol = 0,
                   estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                   perfFUN = mse, genFUN = abs.mean, genFUN.add = list(), selFUN = which.min, ...) {
  out <- steps <- list()
  
  ## Step 0
  # As remaining.changes = c("") yields unchanged formula.
  step.count <- 0
  steps[[getStepName(step.count)]] <- mp.step(formula = formula, data = data, remaining.changes = c(""), 
                                     st.i = st.i, st.u = st.u, folds = folds, recal.int = recal.int, 
                                     retest = FALSE,
                                     estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                                     perfFUN = mse, genFUN = genFUN, genFUN.add = genFUN.add, 
                                     selFUN = selFUN, ...)
  steps[[getStepName(step.count)]][["step.count"]] <- step.count
  out[["best.step"]] <- getStepName(step.count)
  current.model <- mp.step.get.best(steps[[1]])
  current.model[["remaining.changes"]] <- remaining.changes
  gen.diff <- Inf # TBI

  
  if (!identical(length(remaining.changes), 0L))
    repeat {
      ## Loop management
      if (isTRUE(gen.diff <= tol)) { # TBI
        out[["stop.reason"]] <- "improvement <= tolerance."
        break
      }
      if (isTRUE(length(current.model[["remaining.changes"]]) <= 0)) {
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
                                                  remaining.changes = current.model[["remaining.changes"]],
                                                  st.i = st.i, st.u = st.u, folds = folds, recal.int = recal.int,
                                                  retest = retest,
                                                  estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL,
                                                  perfFUN = mse, genFUN = genFUN, genFUN.add = genFUN.add,
                                                  selFUN = selFUN, ...)
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
  out[["nobs.strata"]] <- sapply(out[["global.model"]][["stratified.fit"]], nobs)
  out[["nobs"]] <- sum(out[["nobs.strata"]])
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

mp.step.get.change <- function(step, ...)
  mp.step.get.best(step)[["changed"]]

# Perform one step in the fitting process of mp.fit
# formula formula to start with
# data data.frame, containing dev and val data
# remaining.changes predictor terms to add or remove
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
mp.step <- function(formula, data, remaining.changes, st.i, st.u, folds, recal.int = FALSE, 
                    two.stage = TRUE, retest = FALSE,
                    estFUN = glm, metaFUN = urma, meta.method = "DL", predFUN = NULL, 
                    perfFUN = mse, genFUN = abs.mean, genFUN.add = list(), 
                    selFUN = which.min, ...) {
  cv <- out <- list()
  out[["start.formula"]] <- formula
  
  for (fc in seq_along(remaining.changes) )
  {
    change <- remaining.changes[fc]
    
    # Produce formula for changes and no changes:
    if (identical(remaining.changes, "")) {
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
                        predFUN = predFUN, perfFUN = perfFUN, genFUN = genFUN,
                        genFUN.add = genFUN.add, ...)
    # Save changes
    cv[[name]][["remaining.changes"]] <- if (retest) remaining.changes else remaining.changes[-fc]
    cv[[name]][["changed"]] <- change
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
  out <- c(cv.dev, mp.meta.fit(stratified.fit = cv.dev[["stratified.fit"]], 
                               metaFUN = metaFUN, meta.method = meta.method) )
  class(out) <- c("mp.global", class(out))
  out
}

# #' @author Valentijn de Jong
# #' @method family   mp.global
# #' @export
# family.mp.global <- function(object, ...)
#   object$family

# Fitted? Predict?

#' @author Valentijn de Jong
#' @method predict   mp.global
#' @export
predict.mp.global <- function(object, newdata = NULL, strata = NULL, type = "response", 
                             recal.int = NULL, center = NULL, ...) 
  predict.metapred(object = object, newdata = newdata, strata = strata, type = type, 
                    recal.int = recal.int, center = center, ...)

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
                  perfFUN = mse, genFUN = abs.mean, genFUN.add = list(), ...) {
  out <- mp.cv.dev(formula = formula, data = data, st.i = st.i, st.u = st.u, folds = folds, two.stage = two.stage,
                   estFUN = estFUN, metaFUN = metaFUN, meta.method = meta.method, ...)
  
  out <- mp.cv.val(cv.dev = out, data = data, st.i = st.i, folds = folds, recal.int = recal.int, two.stage = two.stage,
                   estFUN = glm, predFUN = predFUN, perfFUN = mse, genFUN = genFUN, genFUN.add = genFUN.add, ...)
  
  class(out) <- c("mp.cv", class(out))
  out
}

# Obtain fitted values of a (\code{subset}ted from metapred) mp.cv object
# object mp.cv
# returns list of fitted values, with length equal to number of strata. 
#' @author Valentijn de Jong
#' @method fitted   mp.cv
#' @export
fitted.mp.cv <- function(object, ...) {
  data <- object[["data"]]
  if (is.null(data))
    stop("Use subset(metapred()) to obtain the cv.")
  strata <- object$strata
  center <- object$options$center
  if (center)
    data <- centerCovs(data = data, y.name = f2o(formula(object)), cluster.name = strata)

  fitted.mp.cv.dev(object, data = data, folds = object$folds, st.i = data[[object$strata]], ...)
}

# fitted.mp.cv.val <- function(object, newdata, folds, st.i, recal.int = FALSE, predFUN = NULL, ...) 
#   fitted.mp.cv.dev(object = object, newdata = newdata, folds = folds, st.i = st.i, predFUN = predFUN)

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
                      estFUN = glm, predFUN = NULL, perfFUN = mse, add.perfFUN = list(), # add.perfFUN = list()
                      genFUN = abs.mean, genFUN.add = list(), ...) {
  dots <- list(...)
  
  # Recalibrate?
  if (isTRUE(recal.int))
    cv.dev <- mp.cv.recal(cv.dev = cv.dev, newdata = data, estFUN = estFUN, folds = folds)
  
  # Predict outcome # Moved to fitted.mp.cv.dev()
  # predictMethod <- getPredictMethod(fit = cv.dev, two.stage = two.stage, predFUN = predFUN, ...)
  # p <- list()
  # for (i in seq_len(cv.dev[["n.cv"]]))
  #   p[[i]] <- predictMethod(object = cv.dev, newdata = data[folds[["val"]][[i]] == st.i, ], type = "response",
  #                           b = coef(cv.dev[["cv"]][[i]]),
  #                           f = formula(cv.dev), two.stage = two.stage, ...)
  p <- fitted.mp.cv.dev(object = cv.dev, data = data, folds = folds, st.i = st.i, predFUN = predFUN) # TBI
  
  # Necessary for performance computation
  outcome <- f2o(formula(cv.dev))
  
  cv.dev[["perf"]] <- data.frame(matrix(ncol = 2, nrow = cv.dev[["n.cv"]]))
  cv.dev[["nobs.val"]] <- sapply(p, length)
  row.names(cv.dev[["perf"]]) <- names(cv.dev[["cv"]])
  colnames(cv.dev[["perf"]])  <- c("val.strata", "perf")
  
  # Compute performance for predictor selection
  # perFUN receives args that might be useful, but not necessary for defaults
  for (i in seq_len(cv.dev[["n.cv"]])) 
    cv.dev[["perf"]][i, ] <- c(getclName(folds[["val"]][[i]]), 
                               perfFUN(p[[i]], data[folds[["val"]][[i]] == st.i, outcome],
                                       data = data, fit = cv.dev[["cv"]][[i]], ...)) 
  cv.dev[["perf"]][, "perf"] <- as.numeric(cv.dev[["perf"]][, "perf"])
  
  # # Compute additional performance measures.
  # TBI
  
  # And finally, the generalizability (and mean performance)
  cv.dev[["gen"]]       <- genFUN(x   = cv.dev[["perf"]][, "perf"], data = data, N = nrow(data), n = cv.dev[["nobs.val"]], ...) 
  cv.dev[["mean.perf"]] <- abs.mean(x = cv.dev[["perf"]][, "perf"], ...)
  
  # Optionally, additional generalizability measures
  if (length(genFUN.add <- dots[["genFUN.add"]]))
    if (!is.list(genFUN.add))
      genFUN.add <- list(genFUN.add)
  
  gen.add <- list()
    for (fun.id in seq_along(genFUN.add))
      gen.add[[names(genFUN.add)[[fun.id]]]] <- genFUN.add[[fun.id]](x = cv.dev[["perf"]][, "perf"], data = data, N = nrow(data), n = cv.dev[["nobs.val"]], ...)
  
  cv.dev[["gen.add"]] <- gen.add
    
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
  out[["cv"]] <- mp.cv.meta.fit(stratified.fit = out[["stratified.fit"]], folds = folds, 
                                metaFUN = metaFUN, meta.method = meta.method)
  
  out[["n.cv"]] <- length(out$cv)
  # out[["nobs.cv"]]
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

# Obtain fitted values of a mp.cv.dev object
# object mp.cv.dev object.
# Returns list of predicted values, with length equal to number of strata.
#' @author Valentijn de Jong
#' @method fitted   mp.cv.dev
#' @export
fitted.mp.cv.dev <- function(object, data, folds, st.i, predFUN = NULL, ...) {
  predictMethod <- getPredictMethod(fit = object, two.stage = TRUE, predFUN = predFUN, ...)
  p <- list()
  for (i in seq_len(object[["n.cv"]]))
    p[[i]] <- predictMethod(object = object, 
                            newdata = data[folds[["val"]][[i]] == st.i, ], 
                            type = "response",
                            b = coef(object[["cv"]][[i]]),
                            f = formula(object), 
                            two.stage = TRUE, 
                            ...)
  # If some other fold method than l1o is applied, this may otherwise cause an error.
  # better method TBI.
  if (all(unlist(lapply(folds$val, length)) == 1))
    names(p) <- unlist(folds$val)
  p
}

# #' @author Valentijn de Jong
# #' @method family   mp.cv.dev
# #' @export
# family.mp.cv.dev <- function(object, ...)
#   object$family

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
                                             estFUN = estFUN, 
                                             family = if (!is.null(cv.dev$family)) cv.dev$family else NULL)
  cv.dev
}

# Make new cv.models for mp.meta.fit for metapred
# This function selects the right fold for development, and calls the fitting function.
# stratified.fit list of mp.stratum.fit objects
# folds list of fold divisions, as given by l1o, or bootstrap in the utils.
# metaFUN function for estimating meta-analytic models, e.g. urma (this file)
# meta.method options for metaFUN
# Returns an object of class mp.cv.meta.fit, which is a list of meta-analytic prediction models
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
  variances <- t(as.data.frame(lapply(stratified.fit, `[[`, "variances", drop = FALSE)))
  
  meta <- metaFUN(coefficients = b, variances = variances, method = meta.method) 
  
  out[["coefficients"]] <- meta$coefficients
  out[["variances"]]    <- meta$variances
  out[["nobs.strata"]]  <- sapply(stratified.fit, nobs)
  out[["nobs"]] <- sum(out[["nobs.strata"]])
  
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
  out[["variances"]]    <- getVars(fit)
  out[["covar"]]        <- getCoVars(fit)
  out[["vcov"]]         <- vcov(fit)
  out[["nobs"]]         <- nobs(fit, use.fallback = TRUE)
  
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

# Too many args that might not be available for this method. I have abandoned it.
# predict.mp.stratum.fit <- function(object, newdata = NULL, family = NULL, 
# formula = NULL, method = NULL, type = "response") {
#   if (!is.null(family))
#     object$family <- family
#   if (is.null(formula))
#     stop("predict.mp.stratum.fit needs a formula")
#   else
#     object$formula <- formula
#   if (is.null(method))
#     stop("predict.mp.stratum.fit needs a prediction method")
# 
#   method(object = object, newdata = newdata, type = type)
# }

# As I would have had to implement it a million times:
family.default <- function(object, ...) 
  object$family

#' Standard errors and variances
#' 
#' Obtain standard errors or variances of a model fit
#' 
#' @aliases variances
#' 
#' @author Valentijn de Jong
#' 
#' @usage se(object, ...)
#' variances(object, ...)
#' 
#' @param object A model fit object
#' @param ... other arguments
#' 
#' @return For \code{se} the standard errors of \code{object}, and for 
#' \code{variances} the variances.
#' @export
se <- function(object, ...)
  UseMethod("se", object)

#' @export
se.default <- function(object, ...)
  sqrt(variances(object, ...))

#' @export
variances <- function(object, ...)
  UseMethod("variances", object)

#' @export
variances.default <- function(object, ...) {
  if (!is.null(v <- object$variances))
    v
  else 
    diag(vcov(object, ...))
}

#' @export
variances.metapred <- function(object, ...)
  variances(object[["global.model"]])

#' @export
variances.mp.cv <- function(object, select = "cv", ...) 
  variances(object[[select, exact = FALSE]])

#' @export
variances.mp.stratified.fit <- function(object, ...)
  t(as.data.frame(lapply(object, variances), check.names = FALSE))

#' @export
variances.mp.cv.meta.fit <- function(object, ...)
  variances.mp.stratified.fit(object, ...)
