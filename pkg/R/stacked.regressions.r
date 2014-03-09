stacked.regressions <- function (models=list(), names.models=NULL, outcome=NA, data=NULL) {
  require(pROC)
  
  if (class(models) != "list") stop("Models should be entered as a list")
  if (is.null(names.models)) names.models = c(paste("model",(1:length(models)),sep=""))
  if (length(models)==0) stop("No models entered!")
  if (length(names.models)!=length(models)) stop("The number of model names does not correspond to the amount of specified models!")
  if (length(models)==0) stop("No models entered!!")
  if (is.null(data)) stop("No validation data entered!")
  if (is.na(outcome)) stop("No outcome variable specified!")
  
  # Calculate variance inflation factor
  davis.vif <- function(fit) {
    v <- fit$vcov[-1,-1, drop=F] #Drop intercept
    nam <- names(fit$coefficients)[-1]
    
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
  }
  
  inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }
  
  #transpose a set of coefficients or model fit into a linear predictor 
  as.lp <- function (fit=NULL, #a named list of regression coefficients or fitted model of the class glm
                     outcome=NA, #the outcome variable
                     data=NULL #validation sample
  ) {      

    if ("lm" %in% class(fit)) {
      lp <- predict(fit, newdata=data, type="response")
      return(lp)
    } 
    if ("numeric" %in% class(fit)) {
      coefs <- fit[which(!is.na(fit))]
      if ("(Intercept)" %in% names(coefs)) {
        xcoefs <- coefs[-match("(Intercept)", names(coefs))] #model frame adds an intercept term by default
      } else {
        xcoefs <- c(0, coefs)
        names(xcoefs)[1] <- "-1"
      }
      
      fmla <- as.formula(paste(outcome,"~",paste(names(xcoefs), collapse="+")))
      dfrTmp <- model.frame(fmla, data)
      x <- as.matrix(model.matrix(fmla, data=dfrTmp))
      lp <- x%*%coefs[match(colnames(x),names(coefs))]   
      return(lp)
    } 
    
    stop ("Model type not supported:", class(fit))  
  }
  
  mle.regconstrained <- function(lp, y, ...) { # it is possible to add control parameters for optim
    model = function(theta, x, y, ...) {
      devs <- dbinom(y, 1, inv.logit(x%*%theta), log = TRUE)
      loglik <- sum(devs)
      return(-loglik)
    }
    v <- cbind(1, lp)
    theta.start <- c(0,rep(0,(dim(lp)[2])))
    names(theta.start) <- c("(Intercept)", colnames(lp))
    lower <- c(-Inf, rep(0, dim(lp)[2]))
    mle <- optim(par=theta.start, fn=model, y=y, x=as.matrix(v), hessian=T, lower=lower, method="L-BFGS-B", ...)
    lp.final <- v %*% mle$par
    out <- list(coefficients=mle$par, vcov=solve(mle$hessian), deviance=2*mle$value, lp=lp.final, fit=mle)
    return(out)
  }
  
  
  ## Generate the linear predictor of each model
  #s = lapply(models, as.lp, outcome, data)
  s <- vapply(models, as.lp, rep(0,dim(data)[1]), outcome, data) #faster
  y <- data[,outcome]
  if (length(models)==1) s = as.array(s) 
  colnames(s) <- names.models
  m.stacked <- mle.regconstrained(lp=s, y=y)
  vif <- davis.vif(m.stacked)
  
  ## Calculate apparent meta-model performance
  meta.roc    <- roc(response=y, predictor=m.stacked$lp)
  meta.roc.se <- sqrt(var.roc(meta.roc))
  cstat <- c(meta.roc$auc, meta.roc.se)
  names(cstat) = c("c-stat", "se(c-stat)")
  
  out <- list(weights=m.stacked$coefficients, vif=vif, deviance=m.stacked$deviance, performance=cstat, vcov=m.stacked$vcov)
  class(out) <- c("metamodel", class(out))
  
  return(out)
}










