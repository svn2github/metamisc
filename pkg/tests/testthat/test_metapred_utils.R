context("utility functions for metapred")


### Some stuff necessary for testing
n <- 2
td <- data.frame(y = rep(0, n * 4), x = rep(0, n * 4), z = c(rep(0, n * 2), rep(1, n * 2)), s = rep(c(rep(0, n), rep(1, n)), 2))

test_that("Centering within clusters works.", {
  expect_true(is.data.frame(cd <- metamisc:::centerCovs(data = td, y.name = "y", cluster.name = "s") ))
  expect_identical(cd$y, td$y)
  expect_identical(cd$x, rep(0, n * 4))
  expect_identical(cd$z, c(rep(-.5, n * 2), rep(.5, n * 2)))
  expect_identical(cd$s, td$s)
})

test_that("Centering within a single works.", {
  expect_true(is.data.frame(cd <- metamisc:::centerCovs(data = td, y.name = "y", cluster.name = "x") ))
  expect_identical(cd$y, td$y)
  expect_identical(cd$x, td$x)
  expect_identical(cd$z, c(rep(-.5, n * 2), rep(.5, n * 2)) )
  expect_identical(cd$s, rep(c(-.5, -.5, .5, .5), 2) )
})

# Deprecated
# test_that("asDataList and Reduce are complements.", {
#   expect_true(is.list(dl <- metamisc:::asDataList(td, td$z)))
#   expect_identical(td, Reduce(rbind, dl)) # Note that this is not always true. But with these parameters it should.
# })

tds <- 1:20

test_that("l1o produces val and dev", {
  expect_true(is.list(cv.l1o <- metamisc:::l1o(tds)))
  expect_true(is.list(cv.l1o$val))
  expect_true(is.list(cv.l1o$dev))
})

test_that("bootstrap produces val and dev", {
  expect_true(is.list(cv.bootstrap <- metamisc:::bootstrap(tds)))
  expect_true(is.list(cv.bootstrap$val))
  expect_true(is.list(cv.bootstrap$dev))
})

test_that("fixed produces val and dev", {
  expect_true(is.list(cv.fixed <- metamisc:::fixed(tds)))
  expect_true(is.list(cv.fixed$val))
  expect_true(is.list(cv.fixed$dev))
})

test_that("successive produces val and dev", {
  expect_true(is.list(cv.successive <- metamisc:::successive(tds)))
  expect_true(is.list(cv.successive$val))
  expect_true(is.list(cv.successive$dev))
})

### Some stuff necessary for testing
set.seed(8092017)
n <- 100
n.cov <- 3
td <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))
td.ig <- td + 1 # For inverse gaussian and Gamma.

### To be included:
# one-stage
# predFUN.
# Tests for options of predict.metapred

test_that("The predict functions predict accurately.", {
  m.bi <- glm(td, family = binomial)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.bi, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.bi, td, coef(m.bi))) == unlist(m.bi$fitted.values))) # == intentionally ignores names.
  
  m.lm <- lm(td)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.lm, two.stage = TRUE)))
  expect_true(all.equal(unlist(pm(m.lm, td, coef(m.lm))) ,as.matrix(unlist(m.lm$fitted.values)),
                        use.names = F, check.attributes = F)) # also intentionally ignores names.
  
  m.no <- glm(td)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.no, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.no, td, coef(m.no))) == unlist(m.no$fitted.values))) # == intentionally ignores names.
  
  m.gm <- glm(td.ig, family = Gamma)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.gm, two.stage = TRUE)))
  expect_true(all.equal(unlist(pm(m.gm, td.ig, coef(m.gm))) ,as.matrix(unlist(m.gm$fitted.values)),
                        use.names = F, check.attributes = F)) # also intentionally ignores names.
  
  m.ig <- glm(td.ig, family = inverse.gaussian)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.ig, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.ig, td.ig, coef(m.ig))) == unlist(m.ig$fitted.values))) # == intentionally ignores names.
  
  m.po <- glm(td, family = poisson)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.po, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.po, td, coef(m.po))) == unlist(m.po$fitted.values))) # == intentionally ignores names.
  
  m.q <- glm(td, family = quasi)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.q, two.stage = TRUE)))
  expect_true(all.equal(unlist(pm(m.q, td, coef(m.q))) ,as.matrix(unlist(m.q$fitted.values)),
                        use.names = F, check.attributes = F)) # also intentionally ignores names.
  
  m.qb <- glm(td, family = quasibinomial)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.qb, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.qb, td, coef(m.qb))) == unlist(m.qb$fitted.values))) # == intentionally ignores names.
  
  m.qp <- glm(td, family = quasipoisson)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.qp, two.stage = TRUE)))
  expect_true(all.equal(unlist(pm(m.qp, td, coef(m.qp))) ,as.matrix(unlist(m.qp$fitted.values)),
                        use.names = F, check.attributes = F)) # also intentionally ignores names.
})
