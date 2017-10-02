# computeRecal is an internal function for recalibrate() and metapred() to recalibrate the intercept.
# recalibrate() is to be used by predict.metapred() and is exported.
# object Model fit object, of class glm, lm or metapred
# newdata Data to use for recalibration.


computeRecal <- function(object, newdata, b = NULL, f = ~ 1, estFUN = NULL,  ...) {
  if (is.null(b)) b <- coef(object)
  if (is.null(estFUN)) {
    if (inherits(object, "metapred"))
      estFUN <- object$FUN$estFUN
    else estFUN <- as.character(class(object)[[1]])
  }
  estFUN <- match.fun(estFUN)

  # Make offset (linear predictor)
  pred.f  <- formula(object)
  pred.mf <- as.matrix(stats::model.frame(formula = pred.f, data = newdata))
  pred.mf[ , 1] <- 1
  os <- pred.mf %*% b

  # offset must be in newdata.
  osdata   <- cbind(newdata, os)
  f <- update.formula(formula(object), f)
  if (is.null(object$family))
    br <- coef(estFUN(f, data = osdata, offset = os))
  else
    br <- coef(estFUN(f, data = osdata, offset = os, family = object$family))
  i <- match(names(br), names(b))
  b[i] <- b[i] + br
  b
}

# computeRecal(g, d3, estFUN = glm)

# \code{recalibrate} assumes coefficients are stored in \code{object$coefficients}
# and that \code{estFUN} accepts an \code{offset} argument.

#' Recalibrate a Prediction Model
#'
#' \code{recalibrate} is used to recalibrate a prediction model of classes \code{metapred, glm} or  \code{lm}.
#'
#' @param object A model fit object to be recalibrated, of class \code{metapred, glm} or \code{lm}, and more.
#' @param newdata data.frame containing new data set for updating.
#' @param f formula. What part of the model should be updated? Default: intercept only. Left-hand side may
#' be left out.
#' @param estFUN Function for model estimation. If left \code{NULL}, the function is automatically retrieved
#' for \code{metapred} objects. For other objects, the function with name corresponding to the first class
#' of the object is taken. E.g. \code{glm()} for \code{glm} objects.
#' @param ... Optional arguments to pass to \code{estFUN}.
#'
#' @details Currently only the coefficients are updated and leaves  the variances
#' and other aspects untouched. For updating the entire model and all its statistics, see \link[stats]{update}.
#'
#' @return Updated model fit object, of the same class as \code{object}. Generally, updated coefficients can
#' be retrieved with \code{coef()}.
#'
#' @export
recalibrate <- function(object, newdata, f = ~ 1, estFUN = NULL, ...) {
  call <- match.call()
  if (is.null(object$original.coefficients))
    object$original.coefficients <- list()
  if (!is.list(object$original.coefficients))
    stop("object is incompatible with recalibrate.")
  f <- as.formula(f)

  object$original.coefficients[[length(object$original.coefficients) + 1]] <- coef(object)
  object$coefficients <- computeRecal(object = object, newdata = newdata, f = f, estFUN = estFUN, ...)

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
#     if (is.null(object$original.coefficients))
#       object$original.coefficients <- list()
#     if (!is.list(object$original.coefficients))
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
#     object$original.coefficients[[length(object$original.coefficients) + 1]] <- coef(object)
#     object$coefficients[1] <- int
#   }
#
#   if (isTRUE(coefficients))
#     stop("coefficient recalibration is not implemented yet.")
#
#   if (isTRUE(intercept) || isTRUE(coefficients))
#     if (is.call(object$call))
#     {
#       object$original.call <- object$call
#       object$call <- call
#     }
#   object
# }
