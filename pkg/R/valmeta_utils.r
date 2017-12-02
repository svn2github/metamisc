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
                                                                                                                               