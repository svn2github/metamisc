context("Testing calculation of AD (O:E ratio)")
skip_on_cran()

library(lme4)
data(EuroSCORE)


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


