logit <- function(x) { if(is.numeric(x))  log(x/(1-x)) else stop("x is not numeric!") }

inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }


validate <- function(x, ds.ipd, time.calibration=NA) {
  
  cal.intercept <- function(y, lp, family) {
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
  
  
  out <- list()
  if ("glm" %in% class(x)) {
    lp <- as.numeric(predict(x, newdata=ds.ipd, type="link")) #calculate linear predictor
    yhat <- as.numeric(predict(x, newdata=ds.ipd, type="response"))
    outcome <- all.vars(formula(x))[1]
    y <- ds.ipd[,outcome]
    
    # Model discrimination
    roc.rule = roc(response=y, predictor=lp)
    
    num.observed <- sum(y)
    num.expected <- sum(yhat)
    p.observed <- num.observed/length(y)
    
    # num events
    events.results <- c(num.observed, num.expected)
    names(events.results) <- c("num.observed", "num.expected")
    
    # O/E ratio and standard error
    OE <- num.observed/num.expected
    se.lnOE <- sqrt((1-p.observed)/num.observed)
    OEresults <- c(OE, exp(log(OE)+qnorm(0.025)*se.lnOE), exp(log(OE)+qnorm(0.975)*se.lnOE))
    names(OEresults) <- c("O:E", "2.5%CI", "97.5%CI")
    
    #calibration slope and intercept
    m.intercept <- cal.intercept(y, lp, x$family)
    m.slope <- cal.slope(y, lp, x$family)
    
    predictions = as.data.frame(cbind(lp,yhat,y))
    
    cal <- list(events=events.results, OE=OEresults, slope=m.slope, intercept=m.intercept)
    out <- list(coefficients=coefficients(x), family=x$family, predictions=predictions, roc=roc.rule, cal=cal, ds.ipd=ds.ipd)
  } else {
    stop("Invalid model class!!")
  }
  ## TODO: write code for time-to-event models 
  
  class(out) <- "validation"
  return(out)
}


#as.model <- function(coefficients, vars, model.class="glm", formula) {
#  model <- list()
#  if (model.class=="glm") {
#    model$coefficients <- coefficients
#    model$formula <- formula
#}  
#}


