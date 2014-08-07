pmmodel <- function(coefficients, formula, family=binomial()) {
  UseMethod("pmmodel")
}


pmmodel.default <- function(coefficients, formula, family=binomial()) {
  if (class(formula) != "formula") stop("Invalid formula object!")
  if (class(family) != "family") stop("Invalid family object!")
  if (! ("(Intercept)" %in% names(coefficients))) warning("The model does not have an intercept term!")
  out <- list()
  out$coefficients <- coefficients
  out$formula <- formula
  out$family <- family
  class(out) <- "pmmodel"
  return(out)
}