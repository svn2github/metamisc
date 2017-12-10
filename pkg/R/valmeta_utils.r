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


restore.oe.OEN <- function(O, E, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  if (model == "normal/identity") {
    cc <- which(E==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    N[cc] <- N[cc]+correction
    out <- O/E
  } else if (model %in% c("normal/log", "poisson/log")) {
    cc <- which(E==0 | O==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    N[cc] <- N[cc]+correction
    out <- (log(O)-log(E))
  } else {
    out <- rep(NA, length(O)) 
  }
  
  # Apply extrapolation or omit study where t.val != t.ma
  if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
    out[which(t.val!=t.ma)] <- NA
  } else if (t.extrapolate) {
    out[which(t.val!=t.ma)] <- restore.oe.PoPe(Po=O/N, Pe=E/N, t.extrapolate=T, t.ma=t.ma, t.val=t.val, model=model)[which(t.val!=t.ma)]
  }
  
  return (out)
}

restore.oe.PoPe <- function (Po, Pe, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  if (missing(t.val)) {
    t.val <- rep(NA, length(Po))
  }
  
  if (model == "normal/identity") {
    out <- Po/Pe
    
    # Apply extrapolation using Poisson distribution
    if (t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      Po.new <- 1-exp(t.ma*log(1-Po)/t.val)
      Pe.new <- 1-exp(t.ma*log(1-Pe)/t.val)
      out[which(t.val!=t.ma)] <- (Po.new/Pe.new)[which(t.val!=t.ma)]
    } else if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  } else if (model %in% c("normal/log", "poisson/log")) {
    out <- log(Po)-log(Pe)
    
    if (t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      Po.new <- 1-exp(t.ma*log(1-Po)/t.val)
      Pe.new <- 1-exp(t.ma*log(1-Pe)/t.val)
      out[which(t.val!=t.ma)] <- log(Po.new/Pe.new)[which(t.val!=t.ma)]
    } else if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  } else {
    out <- rep(NA, length(Po)) 
  }
  
  return (out)
}

restore.oe.OPeN <- function(O, Pe, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.OEN(O=O, E=Pe*N, N=N, correction=correction, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}

restore.oe.EPoN <- function(E, Po, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.OEN(O=Po*N, E=E, N=N, correction=correction, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}

restore.oe.PoPeN <- function (Po, Pe, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.OEN(O=Po*N, E=Pe*N, N=N, correction=correction, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}

# Restore OE ratio from calibration-in-the-large
restore.oe.citl <- function(citl, O, Po, N, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  if (missing(t.val)) {
    t.val <- rep(NA, length(citl))
  }
  
  Po[is.na(Po)] <- (O/N)[is.na(Po)]
  
  # Apply extrapolation using Poisson distribution
  if (t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
    Po[which(t.val!=t.ma)] <- (1-exp(t.ma*log(1-Po)/t.val))[which(t.val!=t.ma)]
  } else if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
    Po[which(t.val!=t.ma)] <- NA # Omit studies where follow-up duration is improprer
  }
  
  if (model == "normal/identity") {
    return(-(exp(citl)*(Po)-exp(citl)-(Po)))
  }
  if (model %in% c("normal/log", "poisson/log")) {
    return (log(-(exp(citl)*(Po)-exp(citl)-(Po))))
  }
  return (NA)
}

restore.oe.var.seOE1 <- function(se, OE, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  if (model == "normal/identity") {
    out <- se**2
    
    if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  } else if (model %in% c("normal/log", "poisson/log")) {
    out <- (se**2)/(OE**2) # Equation 16 in appendix BMJ paper
    
    if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  }
  return(out)
}

restore.oe.var.seOE2 <- function(se, O, E, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.var.seOE1(se=se, OE=(O/E), t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}
  
#In most situations, O and E are reported separately without any estimate of uncertainty.
#In the following derivations, we regard E as a fixed constant. We treat O as a binomially distributed
#variable since O is given as the number of successes (events) from N subjects
restore.oe.var.OEN <- function(O, E, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  if (model == "normal/identity") {
    cc <- which(E==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    N[cc] <- N[cc]+correction
    Po <- O/N
    
    out <- O*(1-Po)/(E**2)
    
    if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  } else if (model %in% c("normal/log", "poisson/log")) {
    cc <- which(E==0 | O==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    N[cc] <- N[cc]+correction
    Po <- O/N
    
    out <- (1-Po)/(O)
    
    if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  }
  return(out)
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
                                                                                                                               