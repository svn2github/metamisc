#Update SE(c.index) using Method 4 of Newcombe
restore.c.var.hanley <- function(cstat, N.subjects, N.events, restore.method="Newcombe.4", model="normal/logit") {
  n <- N.events #Number of events
  m <- N.subjects-N.events #Number of non-events
  
  if (missing(restore.method)) {
    restore.method <- "Newcombe.4"
  }
  
  if (restore.method=="Hanley" | restore.method=="Newcombe.2") {
    mstar <- m-1
    nstar <- n-1
  } else if (restore.method=="Newcombe.4") {
    mstar <- nstar <- N.subjects/2-1
  } else {
    stop ("Method not implemented yet!")
  }
  
  if (model=="normal/logit") {
    out <- (((1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n*cstat*(1-cstat)))
  } else if (model=="normal/identity") {
    out <- ((cstat*(1-cstat)*(1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n))
  } else {
    stop ("Meta-analysis model not implemented!")
  }
  
  return(out)
}

restore.c.var.se <- function(c.se, cstat, model="normal/logit") {
  if (model=="normal/identity") {
    return (c.se**2)
  }
  if (model=="normal/logit") {
    return((c.se/(cstat*(1-cstat)))**2)
  }
  stop("Invalid link function!")
}

restore.c.var.ci <- function(ci, level=0.95, model="normal/logit") {
  if (length(dim(ci))>=2) {
    upper <- ci[,2]
    lower <- ci[,1]
  } else {
    upper <- ci[2]
    lower <- ci[1]
  }
  
  if (model=="normal/identity") {
    return(((upper - lower)/(2*qnorm((1-level)/2)))**2)
  }
  if (model=="normal/logit") {
    return(((logit(upper) - logit(lower))/(2*qnorm((1-level)/2)))**2)
  }
  stop("Invalid link function!")
}



# See params of geom_smooth for more details
# Smoothed calibration plot: use  formula = obsy ~ splines::bs(predy, 3)
plotCalibration <- function(predy, obsy, modelname="Model", 
                            formula = obsy ~ predy, 
                            method="glm", se = se, 
                            level=0.95, 
                            fam=binomial, ...) {
  #require(ggplot2)
  #require(ggExtra) ## Do the extra outside the package
  
  # Use formula instead
  
  # Add density plot under gg plot
  #predy<-rnorm(300)
  #obsy<-rt(300,df=10)
  #family <- gaussian
  
  xy <- data.frame(predy, obsy)
  
  scatter <- ggplot(xy, aes(x=predy, y=obsy)) +
    labs(x = "Predicted", y="Observed") + 
    scale_x_continuous(limits=c(min(predy),max(obsy))) + 
    scale_y_continuous(limits=c(min(predy),max(obsy))) #+ 
  
  
  # Add reference line for perfect calibration
  scatter <- scatter + geom_abline(aes(slope=1, intercept=0, linetype="Perfect calibration"), size=1)
  scatter <- scatter + geom_smooth(aes(linetype=modelname), method = method, 
                                   se = se, level=level, method.args = list(family = fam), ...)

  
  #scatter <- scatter + ggMarginal(data = xy, x = "predy", y = "obsy", margins = "x", type="histogram", size=4)
  #scatter <- scatter + labs(x = "Predicted", y="Observed") 
  scatter
  
  
}
                                                                                                                               