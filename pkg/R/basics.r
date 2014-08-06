logit <- function(x) { if(is.numeric(x))  log(x/(1-x)) else stop("x is not numeric!") }

inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }

as.model <- function(coefficients, formula, family=binomial()) {
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

validate <- function(x, ds.ipd, time.calibration=NA) {
  
  cal.intercept <- function(y, lp, family) {
    y <- as.numeric(y)
    lp <- as.numeric(lp)
    predictions <- as.data.frame(cbind(y,lp))
    m.ll <- glm("y~1", offset=lp, data=predictions, family=family)
    results = array(NA, dim=2)
    names(results) = c("estimate", "se")
    results["estimate"] <- coefficients(m.ll)[1]
    results["se"] <- sqrt(vcov(m.ll)[1,1])
    return(results)
  }
  cal.slope <- function(y, lp, family) {
    predictions <- as.data.frame(cbind(y,lp))
    m.ll = glm(as.formula("y~lp"), data=predictions, family=family)
    results = array(NA, dim=2)
    names(results) = c("estimate", "se")
    results["estimate"] <- coefficients(m.ll)[2]
    results["se"] <- sqrt(vcov(m.ll)[2,2])
    return(results)
  }
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
  
  
  out <- list()
  if ("glm" %in% class(x)) {
    lp   <- as.numeric(predict(x, newdata=ds.ipd, type="link")) #calculate linear predictor
    yhat <- as.numeric(predict(x, newdata=ds.ipd, type="response"))
    outcome <- all.vars(formula(x))[1]
    family <- family(x)
    coefs <- coefficients(x)
    y <- ds.ipd[,outcome]
    predictions <- as.data.frame(cbind(lp, yhat, y))
    
    #calibration slope and intercept
    m.intercept <- cal.intercept(y, lp, family(x))
    m.slope <- cal.slope(y, lp, family(x))
  } else if ("pmmodel" %in% class(x)) {
    outcome <- all.vars(x$formula)[1]
    family <- x$family
    coefs <- x$coefficients
    predictions <- calc.lp(coefs, ds.ipd, x$formula)
    
    
    #calibration slope and intercept
    m.intercept <- cal.intercept(predictions$y, predictions$lp, x$family)
    m.slope <- cal.slope(predictions$y, predictions$lp, x$family)
  } else {
    stop("Invalid model class!!")
  }
  
  if (sum(c("pmmodel", "glm") %in% c("glm", "lm"))>0)
  {
    # Model discrimination
    roc.rule = roc(response=predictions$y, predictor=predictions$lp)
    
    num.observed <- sum(predictions$y)
    num.expected <- sum(predictions$yhat)
    p.observed   <- mean(predictions$y)
    
    # num events
    events.results <- c(num.observed, num.expected)
    names(events.results) <- c("num.observed", "num.expected")
    
    # O/E ratio and standard error
    OE <- num.observed/num.expected
    se.lnOE <- sqrt((1-p.observed)/num.observed)
    OEresults <- c(OE, exp(log(OE)+qnorm(0.025)*se.lnOE), exp(log(OE)+qnorm(0.975)*se.lnOE))
    names(OEresults) <- c("O:E", "2.5%CI", "97.5%CI")
    
    
    
    cal <- list(events=events.results, OE=OEresults, slope=m.slope, intercept=m.intercept)
    out <- list(coefficients=coefs, family=family, predictions=predictions, roc=roc.rule, cal=cal, ds.ipd=ds.ipd)
  }
    
  ## TODO: write code for time-to-event models 
  
  class(out) <- "validation"
  return(out)
}

print.validation <- function(x) {
  cat("Validation Data\n*************************************\n")
  cat(paste("Study size: ", dim(x$predictions)[1], " subjects (", 
            x$cal$events["num.observed"], " events)\n", sep=""))
  cat("\nDiscrimination\n*************************************\n")
  cat(paste("Area under the ROC curve: ", round(x$roc$auc,3), "\n", sep=""))
  ci.roc <- signif(ci(x$roc), digits = 3)
  cat(paste("95% confidence interval: ", ci.roc[1], "; ", ci.roc[3], "\n", sep=""))
  cat("\nCalibration\n*************************************\n")
  cat(paste("Observed versus expected: "), round(x$cal$OE["O:E"],3), "\n", sep="")
  cat(paste("Calibration-in-the-large: ", round(x$cal$intercept["estimate"],3), "\n"))
  cat(paste("Calibration slope: ", round(x$cal$slope["estimate"],3), "\n"))
}


#as.model <- function(coefficients, vars, model.class="glm", formula) {
#  model <- list()
#  if (model.class=="glm") {
#    model$coefficients <- coefficients
#    model$formula <- formula
#}  
#}


