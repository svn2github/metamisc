# formula to term labels
# utility function for obtaining terms (predictors, or interactions etc) of formula
# f formula
# returns term.labels (predictor names)
f2tl <- function(f) 
  attr(terms(f), "term.labels")

# Add or remove predictors
# f formula
# x predictors to be added or removed.
# NOTE: if one predictor is to be added, and one to be removed, an error is raised.
# Use two calls: updateFormula(updateFormula(f, xremove), xadd) or the other way around.
updateFormula <- function(f, x) {
  m <- x %in% f2tl(f)
  if (any(m)) if (!all(m)) stop("formula f must either contain all or no predictors x") # if any, all x must be contained in f 
  update.formula(f, paste("~ .", paste(if (any(m)) "-" else "+", as.character(x) ,
                                       sep = " ", collapse = " "), sep = " "  ) )
}

# Deprecated
# formulaAdd <- function(f, x) 
#   update.formula(f, paste("~ .", paste("+", as.character(x) , sep = " ", collapse = " "), sep = " "  ) )
# formulaRemove <- function(f, x)
#   update.formula(f, paste("~ .", paste("-", as.character(x) , sep = " ", collapse = " "), sep = " "  ) )


# Differences in vectors of characters
# x, y vector
# returns ordered vector of symmetric difference between x and y
setSymDiff <- function(x, y) 
  sort(unique(c(setdiff(x, y), setdiff(y, x))))

# Differences in formulas
# f, g formula
# returns ordered vector of symmetric difference between f and g 
getFormulaDiffAsChar <- function(f, g)
  setSymDiff(f2tl(f), f2tl(g))