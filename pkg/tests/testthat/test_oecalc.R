context("Testing calculation of AD (O:E ratio)")
skip_on_cran()

library(lme4)
data(EuroSCORE)

test_that ("Coversion between CITL and log O:E ratio", {
  
  # Generate a test dataset
  N <- 1000
  x1 <- rnorm(N, mean=0, sd=1)
  x2 <- rnorm(N, mean=0, sd=2)
  x3 <- rnorm(N, mean=0, sd=1)
  py <- inv.logit(cbind(1, x1, x2, x3) %*% c(-2, 1, 1, 1.5))
  y <- rbinom(N, prob=py, size = 1)
  
  # Fit model
  fit1 <- glm(y ~ x1 + x2 + x3, family = binomial())
  
  # Calculate LP
  pred.lp <- predict(fit1)
  
  # Introduce mis-calibration-in-the-large
  intercept.add <- 1.2
  pred.lp2 <- pred.lp + intercept.add
  
  # Calculate calibration-in-the-large
  fit2 <- glm(y ~ 1, offset = pred.lp2, family = binomial())
  
  # Check if estimate for calibration-in-the-large is correct
  citl <- as.numeric(coefficients(fit2)[1])
  expect_equal(citl, -intercept.add)
  
  # Calculate total O:E ratio
  O <- sum(y)
  E <- sum(inv.logit(pred.lp2))
  OE <- O/E
  
  # Restored Po
  Po.restored <- mean(inv.logit(pred.lp2 + citl))
  expect_equal(O/N, Po.restored)
  
  # Derive Pe from Po and citl
  Po.restored2 <- mean(1/(1+(exp(-pred.lp2)*exp(citl))))
  Po.restored3 <- mean(1/2 + 1/2 * tanh((pred.lp2 + citl)/2))
  
  
  # Po <- 1/2 + 1/2 * mean(tanh((pred.lp2 + citl)/2))
  # (Po-0.5)*2 <- mean(tanh((pred.lp2 + citl)/2))
  # tan((Po-0.5)*2) <- pred.lp2 + citl
  
  
})


test_that("Standard error of log O:E ratio", {
  
  # Poisson distribution for total number of observed events
  logoese1 <- sqrt(1/EuroSCORE$n.events - 1/EuroSCORE$n)
  logoese2 <- sqrt((resoe.O.E.N(O=EuroSCORE$n.events, E=EuroSCORE$e.events, N=EuroSCORE$n, correction = 1/2, g="log(OE)"))[,2])
  logoese3 <- (oecalc(O=n.events, E=e.events, N=n, data=EuroSCORE, g="log(OE)"))$theta.se
  expect_equal(logoese1, logoese2, logoese3)
  
  # Binomial distribution for total number of observed events
  logoese4 <- sqrt(1/EuroSCORE$n.events)
  logoese5 <- sqrt((resoe.O.E(O=EuroSCORE$n.events, E=EuroSCORE$e.events, correction=1/2, g="log(OE)"))[,2])
  logoese6 <- (oecalc(O=EuroSCORE$n.events, E=EuroSCORE$e.events, g="log(OE)"))$theta.se
  expect_equal(logoese4, logoese5, logoese6)
  
  # Mixture of Poisson and Binomial
  n.total <- EuroSCORE$n
  n.total[1:10] <- NA
  logoese7 <- sqrt(1/EuroSCORE$n.events - 1/EuroSCORE$n)
  logoese7[1:10] <- sqrt(1/EuroSCORE$n.events)[1:10]
  logoese8 <- sqrt((resoe.O.E.N(O=EuroSCORE$n.events, E=EuroSCORE$e.events, N=EuroSCORE$n, correction = 1/2, g="log(OE)"))[,2])
  logoese8[1:10] <- sqrt((resoe.O.E(O=EuroSCORE$n.events, E=EuroSCORE$e.events, correction=1/2, g="log(OE)"))[,2])[1:10]
  logoese9 <- (oecalc(O=EuroSCORE$n.events, E=EuroSCORE$e.events, n=n.total, g="log(OE)"))$theta.se
  expect_equal(logoese7, logoese8, logoese8)
  
  # Derive from the 95% confidence interval
  OE <- EuroSCORE$n.events / EuroSCORE$e.events
  tOE.se <- sqrt(1/EuroSCORE$n.events - 1/EuroSCORE$n)
  OE.cilb <- exp(log(OE) + qnorm(0.025) * tOE.se)
  OE.ciub <- exp(log(OE) + qnorm(0.975) * tOE.se)
  tOE.deriv <- resoe.OE.ci(OE=OE, OE.cilb=OE.cilb, OE.ciub=OE.ciub, OE.cilv=rep(0.95,length(OE)), g="log(OE)")
  expect_equal(tOE.se, sqrt(tOE.deriv[,2]))
})


test_that("Study names", {
  oe.est1 <- oecalc(O=n.events, E=e.events, N=n, data=EuroSCORE, g="log(OE)")
  oe.est2 <- oecalc(O=n.events, E=e.events, N=n, data=EuroSCORE, g="log(OE)", slab=Study)
  expect_equal(oe.est1$theta, oe.est2$theta)
})


