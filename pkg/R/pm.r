pm <- function(x, formula, family=binomial()) {
  UseMethod("pm")
}


pm.default <- function(x, formula, family=binomial()) {
  if ("glm" %in% class(x)) {
    out <- x
    class(out) <- c(class(out),"pm")
    return(out)
  }

  if (class(formula) != "formula") stop("Invalid formula object!")
  if (class(family) != "family") stop("Invalid family object!")
  if (family$link=="logit" & !("(Intercept)" %in% names(x))) warning("The model does not have an intercept term!")
  out <- list()
  out$coefficients <- x
  out$formula <- formula
  out$family <- family
  class(out) <- "pm"
  return(out)
}

