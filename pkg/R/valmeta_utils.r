plotForest <- function (vmasum, xlab, refline, ...) {
  inv.logit <- function(x) {1/(1+exp(-x)) }
  
  k <- dim(vmasum$data)[1]
  slab <- c(as.character(vmasum$slab))
  yi <- c(vmasum$data[,"theta"])
  ci.lb <- c(vmasum$data[,"theta.95CIl"])
  ci.ub <- c(vmasum$data[,"theta.95CIu"])
  
  if (vmasum$model=="normal/logit") {
    yi <- sapply(yi, inv.logit)
    ci.lb <- sapply(ci.lb, inv.logit)
    ci.ub <- sapply(ci.ub, inv.logit)
    refline <- NA
  } else if (vmasum$model == "normal/log" | vmasum$model == "poisson/log") {
    yi <- sapply(yi, exp)
    ci.lb <- sapply(ci.lb, exp)
    ci.ub <- sapply(ci.ub, exp)
    refline <- 1
  }
  
  #Sort data
  i.index <- order(yi)
  
  # Add meta-analysis results
  slab <- c(slab[i.index], "Summary Est.")
  yi <- c(yi[i.index], vmasum$results["estimate"])
  ci.lb <- c(ci.lb[i.index], vmasum$results["95CIl"])
  ci.ub <- c(ci.ub[i.index], vmasum$results["95CIu"])
  
  ALL <- data.frame(study=slab, mean=yi, m.lower=ci.lb, m.upper=ci.ub, order=length(yi):1)

    
  # reorder factor levels based on another variable (HPD.mean)
  ALL$study.ES_order <- reorder(ALL$study, ALL$order, mean)       
 
  p <- with(ALL, ggplot(ALL[!is.na(ALL$mean), ], 
              aes(x = study.ES_order, y = mean, ymin = m.lower, ymax = m.upper)) + 
              geom_pointrange() + 
              theme_bw() + 
              coord_flip() + 
              ylab(xlab) + 
              xlab(""))
  
  # Add refline
  if (is.numeric(refline)) {
    p <- p + geom_hline(yintercept = refline,  linetype = "dotted") 
  }
  if (vmasum$model=="normal/logit") {
    p <- p + scale_y_continuous(limits=c(0, 1))
  } else if (vmasum$model == "normal/log" | vmasum$model == "poisson/log") {
    p <- p + scale_y_continuous(trans = "log", breaks = c(0, 0.5, 1, 1.5, 2, 3, 5, 10, 25))
  }
  # Add meta-analysis summary
  g2 <- with(ALL, subset(ALL, study == "Summary Est."))
  g2$pi.upper <- vmasum$results["95PIu"]
  g2$pi.lower <- vmasum$results["95PIl"]

  # Prediction interval
  p <- p + with(g2, geom_segment(data=g2, aes(y = pi.lower, x="Summary Est.", xend = "Summary Est.", yend = pi.upper), linetype=2))

  # Confidence interval
  p <- p + with(g2, geom_segment(data=g2, aes(y = m.lower, x="Summary Est.", xend = "Summary Est.", yend = m.upper), size=1.4))
  
  # Summary estimate
  p <- p + with(g2, geom_point(data=g2, shape=23, size=3, fill="white"))

  p
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
                                                                                                                               