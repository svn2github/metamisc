pm <- function(x, formula, family=binomial()) {
  UseMethod("pm")
}


pm.default <- function(x, formula, family=binomial()) {
  if (missing(x)) stop("The variable 'x' is missing!")
  if (missing(formula)) stop("No formula provided!")
  if (missing(family)) stop("No family provided!")
  if (class(x) != "numeric") stop("Invalid set of coefficients!")
  if (class(formula) != "formula") stop("Invalid formula object!")
  if (class(family) != "family") stop("Invalid family object!")
  
  if (family$link=="logit" & !("(Intercept)" %in% names(x))) warning("The model does not have an intercept term!")
  out <- list()
  out$coefficients <- x
  out$formula <- formula
  out$family <- family
  out$k.upd <- 0
  class(out) <- "pm"
  return(out)
}

update.pm <- function(x, data, method="is", ...) {
  
  if (missing(data)) stop("No updating data provided!")
  if (missing(method)) stop("No updating method provided!")
  
  calc.lp <- function(coefficients, ds, fmla) {
    outcome = rownames(attr(terms(fmla),"factors"))[1]
    dfrTmp = model.frame(ds)
    x <- as.matrix(model.matrix(fmla, data=dfrTmp))
    out <- as.data.frame(array(NA, dim=c(dim(x)[1], 3)))
    colnames(out) <- c("lp", "yhat", "y")
    out$y <- ds[,match(outcome,colnames(ds))]
    
    names.beta <- if(class(coefficients)=="matrix") colnames(coefficients) else names(coefficients)
    beta = as.numeric(coefficients)
    beta = as.matrix(beta)
    beta[which(is.na(beta))] = 0 #Replace NAs by zero
    
    out$lp = x%*%beta[match(colnames(x), names.beta)]
    out$yhat <- inv.logit(out$lp)
    return(out)
  }
  
  predictions <- calc.lp (x$coefficients, data, x$formula)
  
  if (method=="i") {
    m.ll <- glm("y~1", offset=lp, data=predictions, family=x$family)
    x$k.upd <- 1
    if ("(Intercept)" %in% names(x$coefficients)) {
      x$coefficients["(Intercept)"] <- x$coefficients["(Intercept)"] + coefficients(m.ll)[1]
    } else {
      x$coefficients <- c(coefficients(m.ll)[1], x$coefficients) #Add an intercept term
      x$formula <- update(x$formula, ~ . + 1)
    }
  } else if (method=="is") {
    m.ll = glm(as.formula("y~lp"), data=predictions, family=x$family)
    x$coefficients <- x$coefficients*coefficients(m.ll)[2]
    x$k.upd <- 2
    if ("(Intercept)" %in% names(x$coefficients)) {
      x$coefficients["(Intercept)"] <- x$coefficients["(Intercept)"] + coefficients(m.ll)[1]
    } else {
      x$coefficients <- c(coefficients(m.ll)[1], x$coefficients) #Add an intercept term
      x$formula <- update(x$formula, ~ . + 1)
    }
  }
  
  return(x)
}

