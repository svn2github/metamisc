context("Testing calculation of AD (O:E ratio)")
skip_on_cran()

library(lme4)
data(EuroSCORE)


test_that("Standard error of log O:E ratio", {
  
  # Poisson distribution for total number of observed events
  logoese1 <- sqrt(1/EuroSCORE$n.events - 1/EuroSCORE$n)
  logoese2 <- (oecalc(O=n.events, E=e.events, N=n, data=EuroSCORE, g="log(OE)"))$theta.se
  expect_equal(logoese1, logoese2)
  
  # Binomial distribution for total number of observed events
  logoese3 <- sqrt(1/EuroSCORE$n.events)
  logoese4 <- sqrt((resoe.O.E(O=EuroSCORE$n.events, E=EuroSCORE$e.events, correction=1/2, g="log(OE)"))[,2])
  logoese5 <- (oecalc(O=EuroSCORE$n.events, E=EuroSCORE$e.events, g="log(OE)"))$theta.se
  expect_equal(logoese3, logoese4, logoese5)
  
  # Derive from the 95% confidence interval
  OE <- EuroSCORE$n.events / EuroSCORE$e.events
  tOE.se <- sqrt(1/EuroSCORE$n.events - 1/EuroSCORE$n)
  OE.cilb <- exp(log(OE) + qnorm(0.025) * tOE.se)
  OE.ciub <- exp(log(OE) + qnorm(0.975) * tOE.se)
  tOE.deriv <- resoe.OE.ci(OE=OE, OE.cilb=OE.cilb, OE.ciub=OE.ciub, OE.cilv=rep(0.95,length(OE)), g="log(OE)")
  expect_equal(tOE.se, sqrt(tOE.deriv[,2]))
})


