#' Forest plot
#' 
#' Generate a forest plot by specifying the various effect sizes, confidence intervals and summary estimate.
#' @param theta Numeric vector with effect size for each study
#' @param theta.ci Two-dimensional array specifying the lower bound (first column) and upper bound (second column) of the 
#' confidence interval of the effect sizes
#' @param theta.slab Character vector specifying the study labels
#' @param theta.summary Meta-analysis summary estimate of the effect sizes
#' @param theta.summary.ci Numeric vector specifying the lower bound (first item) and upper bound (second item) of the 
#' confidence interval of the summary estimate
#' @param theta.summary.pi Numeric vector specifying the lower bound (first item) and upper bound (second item) of the 
#' prediction interval of the summary estimate. 
#' @param theme Theme to generate the forest plot. By default, the classic dark-on-light ggplot2 theme is used. 
#' See \link[ggplot2]{theme_bw} for more information.
#' @param xlim The \code{x} limits \code{(x1, x2)} of the forest plot
#' @param xlab Optional character string specifying the X label
#' @param refline Optional numeric specifying a reference line
#' @param label.summary Optional character string specifying the label for the summary estimate
#' @param \dots Additional arguments, which are currently ignored.
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
forest <- function (theta, theta.ci, theta.slab, theta.summary, 
                    theta.summary.ci, 
                    theta.summary.pi=c(NA, NA),
                    theme = theme_bw(),
                    xlim,
                    xlab="", refline=0, label.summary="Summary Estimate", ...) {

  if (missing(theta)) stop("Study effect sizes are missing!")
  if (missing(theta.ci)) stop("Confidence intervals of effect sizes missing!")
  if (missing(theta.slab)) stop("Study labels are missing!")
  
  num.studies <- unique(c(length(theta), dim(theta.ci)[1], length(theta.slab)))
  if (length(num.studies)>1) stop(paste("Mismatch in data dimensions!"))
  
  #Extract data
  yi <- theta
  slab <- theta.slab
  
  #Sort data
  i.index <- order(yi)
  
  # Add meta-analysis results
  slab <- c(slab[i.index], label.summary)
  yi <- c(yi[i.index], theta.summary)
  ci.lb <- c(theta.ci[i.index,1], theta.summary.ci[1])
  ci.ub <- c(theta.ci[i.index,2], theta.summary.ci[2])
  
  ALL <- data.frame(study=slab, mean=yi, m.lower=ci.lb, m.upper=ci.ub, order=length(yi):1)
  
  # reorder factor levels based on another variable (HPD.mean)
  ALL$study.ES_order <- reorder(ALL$study, ALL$order, mean)    
  
  p <- with(ALL, ggplot(ALL[!is.na(ALL$mean), ], 
                        aes(x = study.ES_order, y = mean, ymin = m.lower, ymax = m.upper)) + 
              geom_pointrange() + 
              coord_flip() + 
              theme +
              ylab(xlab) + 
              xlab(""))
  
  if (!missing(xlim)) {
    p <- p + ylim(xlim)
  }
  
  # Add refline
  if (is.numeric(refline)) {
    p <- p + geom_hline(yintercept = refline,  linetype = "dotted") 
  }
  
  # Add meta-analysis summary
  g2 <- with(ALL, subset(ALL, study == label.summary))
  g2$pi.upper <- theta.summary.pi[2]
  g2$pi.lower <- theta.summary.pi[1]
  g2$ci.upper <- theta.summary.ci[2]
  g2$ci.lower <- theta.summary.ci[1]
  
  # Prediction interval
  #p <- p + with(g2, geom_segment(data=g2, aes(y = pi.lower, x="Summary Est.", xend = "Summary Est.", yend = pi.upper), linetype=2))
  if (!is.na(g2$pi.lower) & !is.na(g2$pi.upper)) {
    p <- p +with(g2, geom_errorbar(data=g2, aes(ymin = pi.lower, ymax = ci.lower, x=label.summary), width = 0.5, linetype=2))
    p <- p +with(g2, geom_errorbar(data=g2, aes(ymin = ci.upper, ymax = pi.upper, x=label.summary), width = 0.5, linetype=2))
  }

  # Confidence interval
  p <- p +with(g2, geom_errorbar(data=g2, aes(ymin = ci.lower, ymax = ci.upper, x=label.summary), width = 0.5, size=1.3))
  #p <- p + with(g2, geom_segment(data=g2, aes(y = m.lower, x="Summary Est.", xend = "Summary Est.", yend = m.upper), size=1.4))
  
  # Summary estimate
  p <- p + with(g2, geom_point(data=g2, shape=23, size=3, fill="white"))
  
  p
}