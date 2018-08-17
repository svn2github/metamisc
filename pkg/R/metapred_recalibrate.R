# computeRecal is an internal function for recalibrate() and metapred() to recalibrate the 
# intercept and/or coefficients.
# recalibrate() is to be used by predict.metapred() and is exported.
# object Model fit object, of class glm, lm or metapred
# newdata Data to use for recalibration.
# wholefit should whole fit be returned?

# To be replaced with one or two more basic functions: 1 for whole fit, 1 for returning coefficients only.

computeRecal <- function(object, newdata, b = NULL, f = ~ 1, estFUN = NULL, f.orig = NULL, wholefit = FALSE,  ...) {
  dots <- list(...)
  if (is.null(b)) b <- coef(object)
  if (is.null(estFUN)) {
    if (inherits(object, "metapred"))
      estFUN <- object$FUN$estFUN
    else estFUN <- as.character(class(object)[[1]])
  }
  estFUN <- match.fun(estFUN)

  # Make offset (linear predictor)
  X <- model.matrix(if (is.null(f.orig)) formula(object) else as.formula(f.orig), data = newdata)
  lp <- X %*% b

  # offset must be in newdata.
  osdata   <- cbind(newdata, lp)
  f <- update.formula(formula(object), formula(f))
  if (is.null(object$family))
  { 
    if (!is.null(dots[["family"]])) {
      refit <- estFUN(f, data = osdata, offset = lp, family = dots[["family"]])
    } else 
      refit <- estFUN(f, data = osdata, offset = lp)
  } else 
    refit <- estFUN(f, data = osdata, offset = lp, family = object$family)
  
  if (isTRUE(wholefit))
    return(refit)
  else
    return(coef(refit))
}

# p vector of predicted probs or outcomes
# y outcome vector
# estFUN estimation function
# family family
# which is "intercept" or "slope".
# ... For compatibility only
pred.recal <- function(p, y, estFUN, family = gaussian, which = "intercept", ...) {
  if (is.character(family))  # Taken directly from glm()
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  estFUN <- match.fun(estFUN)
  lp <- family$linkfun(p)
  alpha <- rep(mean(lp), length(lp))
  data <- data.frame(y = y, lp = lp, alpha = alpha)
  
  if (identical(which, "intercept"))
    return(estFUN(formula = y ~ 1,  data = data, family = family, offset = lp))
  if (identical(which, "slope"))
    return(estFUN(formula = y ~ lp, data = data, family = family, offset = alpha))
  if (identical(which, "add.slope"))
    return(estFUN(formula = y ~ lp - 1, data = data, family = family, offset = lp))
}

# gl.recal <- pred.recal(predict(gl, type = "response"), gl$data$X1, "glm", binomial)
# int.recal <- (gl.recal)[[1]][[1]]
# all.equal(int.recal, 0) # intercept recalibration works.

# gl.recal <- pred.recal(predict(gl, type = "response"), gl$data$X1, "glm", binomial, which = "slope")
# slo.recal <- gl.recal[[1]][[2]]
# all.equal(slo.recal, 1) # Perfect

# gl.recal <- metamisc:::pred.recal(predict(gl, type = "response"), gl$data$X1, "glm", binomial, which = "add.slope")
# slo.recal <- gl.recal[[1]][[2]]
# all.equal(slo.recal, 1) #



## Split this into two functions.
# shrink <- function(object, newdata, method = "chisq", b = NULL, estFUN = NULL, ...) {
#   call <- match.call()
#   if (is.null(b)) b <- coef(object)
#   if (is.null(object$orig.coef))
#     object$orig.coef <- list()
#   if (!is.list(object$orig.coef))
#     stop("object is incompatible")
#   
#   if (identical(method, "chisq")) {
#     cs <- object$null.deviance - object$deviance
#     df <- object$df.null - object$df.residual
#     object$shrinkage.factor <- (cs - df) / cs
#   } 
#   else {
#   f <- formula( ~ lp)
#   
#   object$orig.coef[[length(object$orig.coef) + 1]] <- coef(object)
#   br <- computeRecal(object = object, newdata = newdata, f = f, estFUN = estFUN, ...)
#   
#   object$shrinkage.factor <- br[2] + 1
#   b[1] <- br[1]
#   b[-1] <- b[-1] * object$shrinkage.factor
#   object$coefficients <- b
#   }
# 
#   
#   if (is.call(object$call))
#   {
#     object$original.call <- object$call
#     object$call <- call
#   }
#   object
# }





# computeRecal(g, d3, estFUN = glm)

# \code{recalibrate} assumes coefficients are stored in \code{object$coefficients}
# and that \code{estFUN} accepts an \code{offset} argument.

#' Recalibrate a Prediction Model
#'
#' \code{recalibrate} is used to recalibrate a prediction model of classes \code{metapred, glm} or  \code{lm}.
#'
#' @param object A model fit object to be recalibrated, of class \code{metapred, glm} or \code{lm}, and more.
#' @param newdata data.frame containing new data set for updating.
#' @param f formula. Which coefficients of the model should be updated? Default: intercept only. Left-hand side may
#' be left out. See \link[stats]{formula} for details.
#' @param estFUN Function for model estimation. If left \code{NULL}, the function is automatically retrieved
#' for \code{metapred} objects. For other objects, the function with name corresponding to the first class
#' of the object is taken. E.g. \code{glm()} for \code{glm} objects.
#' @param ... Optional arguments to pass to \code{estFUN}.
#'
#' @details Currently only the coefficients are updated and the variances and other aspects are left untouched. 
#' For updating the entire model and all its statistics, see \link[stats]{update}.
#'
#' @return Recalibrated model fit object, of the same class as \code{object}. Generally, updated coefficients can
#' be retrieved with \code{coef()}.
#' 
#' @examples
#' data(DVTipd)
#' DVTipd$cluster <- 1:4 # Add a fictional clustering to the data set.
#' # Suppose we estimated the model in three studies: 
#' DVTipd123 <- DVTipd[DVTipd$cluster <= 3, ]
#' mp <- metamisc:::metapred(DVTipd123, strata = "cluster", f = dvt ~ vein + malign, 
#' family = binomial)
#' # and now want to recalibrate it for the fourth:
#' DVTipd4 <- DVTipd[DVTipd$cluster == 4, ]
#' metamisc:::recalibrate(mp, newdata = DVTipd4)
#'
#' @export
recalibrate <- function(object, newdata, f = ~ 1, estFUN = NULL, ...) {
  call <- match.call()
  if (is.null(object$orig.coef))
    object$orig.coef <- list()
  if (!is.list(object$orig.coef))
    stop("object is incompatible with recalibrate.")
  f <- as.formula(f)

  object$orig.coef[[length(object$orig.coef) + 1]] <- coef(object)
  br <- computeRecal(object = object, newdata = newdata, f = f, estFUN = estFUN, ...)
  i <- match(names(br), names(object$coefficients))
  object$coefficients[i] <- object$coefficients[i] + br

  if (is.call(object$call))
  {
    object$original.call <- object$call
    object$call <- call
  }
  object
}

########### Deprecated ###############


# computeInt is an internal function for recalibrate() and metapred() to recalibrate the intercept.
# recalibrate() is to be used by predict.metapred() only.

# computeInt <- function(object, newdata, b = NULL, estFUN = NULL,  ...) {
#   if (is.null(b)) b <- coef(object)
#
#   # Make offset (linear predictor)
#   pred.f  <- formula(object)
#   pred.mf <- as.matrix(stats::model.frame(formula = pred.f, data = newdata))
#   pred.mf[ , 1] <- 1
#   os <- pred.mf %*% b
#
#   # Estimate intercept recalibration, using offset.
#   osdata   <- cbind(newdata, os)
#   recal.f  <- getFormula(data = osdata, predictors = NULL)
#   b[1] + coef(estFUN(recal.f, data = osdata, offset = os, ...))
# }

# recalibrate2 <- function(object, newdata, intercept = TRUE, estFUN = NULL, ...) {
#   call <- match.call()
#   coefficients <- FALSE # for future development.
#   if (isTRUE(intercept) || isTRUE(coefficients))
#   {
#     if (is.null(object$orig.coef))
#       object$orig.coef <- list()
#     if (!is.list(object$orig.coef))
#       stop("object is incompatible with recalibrate.")
#     if (is.null(estFUN)) {
#       if (inherits(object, "metapred"))
#         estFUN <- object$FUN$estFUN
#       else estFUN <- as.character(class(object)[[1]])
#     }
#     estFUN <- match.fun(estFUN)
#   }
#
#   if (isTRUE(intercept))
#   {
#     int <- computeInt(object = object, newdata = newdata, estFUN = estFUN, ...)
#     object$orig.coef[[length(object$orig.coef) + 1]] <- coef(object)
#     object$coefficients[1] <- int
#   }
#
#   if (isTRUE(coefficients))
#     stop("coefficientsficient recalibration is not implemented yet.")
#
#   if (isTRUE(intercept) || isTRUE(coefficients))
#     if (is.call(object$call))
#     {
#       object$original.call <- object$call
#       object$call <- call
#     }
#   object
# }
