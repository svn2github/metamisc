validation <- function(x) {
  UseMethod("validation")
}


validation.default <- function(x, ds.ipd, time.calibration=NA) {
  
  cal.intercept <- function(y, lp, family) {
    y <- as.numeric(y)
    lp <- as.numeric(lp)
    predictions <- as.data.frame(cbind(y,lp))
    m.ll <- glm("y~1", offset=lp, data=predictions, family=family)
    m.ci <- confint(m.ll)
    results = array(NA, dim=4)
    names(results) = c("estimate", "se", "2.5%", "97.5%")
    results["estimate"] <- coefficients(m.ll)[1]
    results["se"] <- sqrt(vcov(m.ll)[1,1])
    results[c("2.5%", "97.5%")] <- m.ci
    return(results)
  }
  cal.slope <- function(y, lp, family) {
    predictions <- as.data.frame(cbind(y,lp))
    m.ll = glm(as.formula("y~lp"), data=predictions, family=family)
    m.ci <- confint(m.ll, "lp")
    results = array(NA, dim=4)
    names(results) = c("estimate", "se", "2.5%", "97.5%")
    results["estimate"] <- coefficients(m.ll)[2]
    results["se"] <- sqrt(vcov(m.ll)[2,2])
    results[c("2.5%", "97.5%")] <- m.ci
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
    print(out)
    return(out)
  }
  
  if (!"pm" %in% class(x))
    stop("Model type not supported!!")
  
  
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
  } else {
    outcome <- all.vars(x$formula)[1]
    family <- x$family
    coefs <- x$coefficients
    predictions <- calc.lp(coefs, ds.ipd, x$formula)
    
    
    #calibration slope and intercept
    m.intercept <- cal.intercept(predictions$y, predictions$lp, x$family)
    m.slope <- cal.slope(predictions$y, predictions$lp, x$family)
  } 
  
  if (sum(c("pm", "glm") %in% c("glm", "lm"))>0)
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
    names(OEresults) <- c("estimate", "2.5%", "97.5%")
    
    
    
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
  cat("\nPerformance\n*************************************\n")
  ci.roc <- signif(ci(x$roc), digits = 3)
  cat(paste("Area under the ROC curve: ", round(x$roc$auc,3), " (95% CI: ",ci.roc[1], "; ", ci.roc[3],")\n", sep=""))
  cat(paste("Observed versus expected: "), round(x$cal$OE["estimate"],3), "\n", sep="")
  cat(paste("Calibration-in-the-large: ", round(x$cal$intercept["estimate"],3), 
            " (95% CI: ", round(x$cal$intercept["2.5%"],3), "; ", round(x$cal$intercept["97.5%"],3),")\n", sep=""))
  cat(paste("Calibration slope: ", round(x$cal$slope["estimate"],3), 
            " (95% CI: ", round(x$cal$slope["2.5%"],3), "; ", round(x$cal$slope["97.5%"],3),")\n", sep=""))
}

plot.validation <- function (x, type="discrimination", ...) {
  
  #This code was copied from package "rms" to avoid unneeded loading of other functionalities
  cut2 <- function (x, cuts, m = 150, g, levels.mean = FALSE, digits, minmax = TRUE, 
                    oneval = TRUE, onlycuts = FALSE) 
  {
    method <- 1
    x.unique <- sort(unique(c(x[!is.na(x)], if (!missing(cuts)) cuts)))
    min.dif <- min(diff(x.unique))/2
    min.dif.factor <- 1
    if (missing(digits)) 
      digits <- if (levels.mean) 
        5
    else 3
    oldopt <- options(digits = digits)
    on.exit(options(oldopt))
    xlab <- attr(x, "label")
    if (missing(cuts)) {
      nnm <- sum(!is.na(x))
      if (missing(g)) 
        g <- max(1, floor(nnm/m))
      if (g < 1) 
        stop("g must be >=1, m must be positive")
      options(digits = 15)
      n <- table(x)
      xx <- as.double(names(n))
      options(digits = digits)
      cum <- cumsum(n)
      m <- length(xx)
      y <- as.integer(ifelse(is.na(x), NA, 1))
      labs <- character(g)
      cuts <- approx(cum, xx, xout = (1:g) * nnm/g, method = "constant", 
                     rule = 2, f = 1)$y
      cuts[length(cuts)] <- max(xx)
      lower <- xx[1]
      upper <- 1e+45
      up <- low <- double(g)
      i <- 0
      for (j in 1:g) {
        cj <- if (method == 1 || j == 1) 
          cuts[j]
        else {
          if (i == 0) 
            stop("program logic error")
          s <- if (is.na(lower)) 
            FALSE
          else xx >= lower
          cum.used <- if (all(s)) 
            0
          else max(cum[!s])
          if (j == m) 
            max(xx)
          else if (sum(s) < 2) 
            max(xx)
          else approx(cum[s] - cum.used, xx[s], xout = (nnm - 
                                                          cum.used)/(g - j + 1), method = "constant", 
                      rule = 2, f = 1)$y
        }
        if (cj == upper) 
          next
        i <- i + 1
        upper <- cj
        y[x >= (lower - min.dif.factor * min.dif)] <- i
        low[i] <- lower
        lower <- if (j == g) 
          upper
        else min(xx[xx > upper])
        if (is.na(lower)) 
          lower <- upper
        up[i] <- lower
      }
      low <- low[1:i]
      up <- up[1:i]
      variation <- logical(i)
      for (ii in 1:i) {
        r <- range(x[y == ii], na.rm = TRUE)
        variation[ii] <- diff(r) > 0
      }
      if (onlycuts) 
        return(unique(c(low, max(xx))))
      flow <- format(low)
      fup <- format(up)
      bb <- c(rep(")", i - 1), "]")
      labs <- ifelse(low == up | (oneval & !variation), flow, 
                     paste("[", flow, ",", fup, bb, sep = ""))
      ss <- y == 0 & !is.na(y)
      if (any(ss)) 
        stop(paste("categorization error in cut2.  Values of x not appearing in any interval:\n", 
                   paste(format(x[ss], digits = 12), collapse = " "), 
                   "\nLower endpoints:", paste(format(low, digits = 12), 
                                               collapse = " "), "\nUpper endpoints:", paste(format(up, 
                                                                                                   digits = 12), collapse = " ")))
      y <- structure(y, class = "factor", levels = labs)
    }
    else {
      if (minmax) {
        r <- range(x, na.rm = TRUE)
        if (r[1] < cuts[1]) 
          cuts <- c(r[1], cuts)
        if (r[2] > max(cuts)) 
          cuts <- c(cuts, r[2])
      }
      l <- length(cuts)
      k2 <- cuts - min.dif
      k2[l] <- cuts[l]
      y <- cut(x, k2)
      if (!levels.mean) {
        brack <- rep(")", l - 1)
        brack[l - 1] <- "]"
        fmt <- format(cuts)
        labs <- paste("[", fmt[1:(l - 1)], ",", fmt[2:l], 
                      brack, sep = "")
        if (oneval) {
          nu <- table(cut(x.unique, k2))
          if (length(nu) != length(levels(y))) 
            stop("program logic error")
          levels(y) <- ifelse(nu == 1, c(fmt[1:(l - 2)], 
                                         fmt[l]), labs)
        }
        else levels(y) <- labs
      }
    }
    if (levels.mean) {
      means <- tapply(x, y, function(w) mean(w, na.rm = TRUE))
      levels(y) <- format(means)
    }
    attr(y, "class") <- "factor"
    if (length(xlab)) 
      label(y) <- xlab
    y
  }
  
  plot.calibration <- function(X, 
                               main="",
                               m=20, #number of observations per group
                               shaded.bg=T,
                               width.bins=10, #height of the bars
                               cex=2,  #text size
                               cex.lab=1, 
                               cex.axis=1, 
                               cex.main=1, 
                               cex.sub=1,
                               shade.col="gray")
  {
    #if (! X$family$link %in% c("logit", "log", "identity")) stop("Model family not supported!")
    if (! X$family$link %in% c("logit")) stop(paste("Calibration curves for prediction models with link \"", X$family$link, "\" currently not supported!")) #currently only support for logistic regression models
    
    #params
    line.par = list(col = "black")
    logmin = 0.03
    logmax = 0.95
    
    
    #c.logit = x%*%beta
    lp = X$predictions$lp #as.vector(c.logit)
    y = X$predictions$y
    
    
    #Calculate performance
    data <- data.frame(y = X$predictions$y, lp = X$predictions$lp)
    gam1 <- glm(y ~ lp, data = data, family=X$family)
    
    x <- seq(min(X$predictions$lp), max(X$predictions$lp), length = 500)
    yy <- predict(gam1, newdata=data.frame(lp = x), type="link", se.fit = TRUE)
    #x <- x[!is.na(yy$fit)]
    #yy$fit <- yy$fit[!is.na(yy$fit)]
    
    se.lower <- (yy$fit + qnorm(0.025) * yy$se.fit)
    se.upper <- (yy$fit + qnorm(0.975) * yy$se.fit)
    xlab <- "Predicted value"
    ylab <- "Observed value"
    xp = x #identity link
    if (X$family$link=="logit") {
      se.lower <- inv.logit(yy$fit + qnorm(0.025) * yy$se.fit)
      se.upper <- inv.logit(yy$fit + qnorm(0.975) * yy$se.fit)
      xp = inv.logit(x) #predicted probability
      xlab = "Predicted probability"
      ylab = "Actual probability"
      ylim = c(-0.15,1)
      xlim = c(0,1)
    } else if (X$family$link=="log") {
      se.lower <- exp(yy$fit + qnorm(0.025) * yy$se.fit)
      se.upper <- exp(yy$fit + qnorm(0.975) * yy$se.fit)
      xp = exp(x)
      xlab = "Predicted number of events"
      ylab = "Observed number of events"
    }
    
    # Plot parameters
    ylim = c(-0.15,1)
    xlim = c(0,1)
    xpol = c(xp, rev(xp), xp[1])
    ypol = c(se.lower, rev(se.upper), se.lower[1])
    xval = xp
    yval = predict(gam1, newdata=data.frame(lp = x), type="response")
    xref = c(0,1)
    yref = c(0,1)
    xaxt = "s"
    yaxt = "s"
    
    plot(0, 0, type = "n",  xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, main=main, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, cex.sub=cex.sub,xaxt=xaxt,yaxt=yaxt)
    
    
    #Draw a background shade
    if (shaded.bg) {
      rect(-5, 0, 5, 1, col=gray.colors(10)[10], border=NA)
      abline(h=seq(0.1, 1, by=0.1), col="white", lwd=2)
    }
    
    
    polygon(xpol, ypol, col = shade.col, border = NA, density = NULL)
    lines(xval, yval, col = line.par$col)
    lines(xref,yref,lty=2)
    
    # Triangles
    p = X$predictions$yhat
    q <- cut2(x=p, levels.mean = TRUE,g=7)
    means <- as.double(levels(q))
    prop <- tapply(y, q, function(x) mean(x, na.rm = TRUE))
    points(means, prop, pch = 2)  
    
    abline(a=0,b=0,lty=3)
    
    lim = c(0, 1)
    #f1 = glm(y~lp, data=data.frame(y=X$predictions$y, lp=X$predictions$lp), family=X$family)
    x1 <- X$predictions$yhat #predict(f1, type="response")
    #x1[p == 0] <- 0
    #x1[p == 1] <- 1
    bins1 <- seq(lim[1], lim[2], length = dim(X$predictions)[1]/width.bins)
    x1 <- x1[x1 >= lim[1] & x1 <= lim[2]]
    f1 <- table(cut(x1, bins1))
    
    #make sure the highest value for pcty equals 0.12
    pcty = as.numeric(f1/sum(f1))#*h.bars
    pcty <- pcty * (0.12/max(pcty))
    
    pctx = rep(0,length(pcty))
    diff = (bins1[2]-bins1[1])/2
    for (j in 1:length(pcty)) {
      if (pcty[j]>0)    rect(bins1[j],-0.15,bins1[j+1],-0.15+pcty[j],col="black")
    }
    box()
  }
  
  if (type=="discrimination") {
    plot(0, 0, type = "n",  xlab = "1-specificity", ylab = "Sensitivity", xlim = c(0,1), ylim = c(0,1), main="Discrimination")
    lines(rev(1-x$roc$sp), rev(x$roc$se), lwd=2)
    lines(c(0,1), c(0,1), lty=2)
  } else if (type=="calibration") {
    plot.calibration(x, main="Calibration", ...)
  }
}