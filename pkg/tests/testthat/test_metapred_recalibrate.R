### Some stuff necessary for testing
set.seed(8092017)
n <- 100
n.cov <- 3
d  <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .7), ncol = n.cov + 1, nrow = n))
td <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))
d3 <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))



### Tests
# # Deprecated
# test_that("computeInt recalibrates", {
#   g <- glm(d, family = binomial)
#   expect_true(is.numeric(r.int <- metamisc:::computeInt(g, newdata = td, estFUN = glm, family = binomial)))
#   expect_true(r.int[[1]] <= coef(g)[1])
# })

test_that("computeRecal recalibrates", {
  g <- glm(d, family = binomial)
  expect_true(is.numeric(r.int <- metamisc:::computeRecal(g, newdata = td, estFUN = glm)))
  expect_true(r.int[[1]] <= coef(g)[1])
})

test_that("recalibrate recalibrates", {
  ### basically recalibrate.glm
  g <- glm(d, family = binomial)
  expect_true(is.list(r <- metamisc:::recalibrate(g, td)))
  expect_true(inherits(r, "glm"))
  r.int <- metamisc:::computeRecal(r, newdata = td, estFUN = glm)
  expect_equal(r.int[[1]], coef(r)[[1]]) # ComputeInt after recalibrate should have no added effect.

  ### basically recalibrate.metapred
  mp <- metapred(td, strata = "X4", recal.int = FALSE)
  mp.r <- recalibrate(mp, newdata = td) # ignores clustering. Hence, intercept changes.
  expect_true(is.list(mp.r))
  expect_true(inherits(mp.r, "metapred"))
})

test_that("metapred with recal.int = T works.", {
  mp.t <- metapred(d, strata = "X4", recal.int = TRUE)
  expect_true(is.list(mp.t$stepwise$coefficients.recal))

  mp.f <- metapred(d, strata = "X4", recal.int = FALSE)
  expect_true(is.null(mp.f$stepwise$coefficients.recal))
})


test_that("predict.metapred with recal.int = T works.", {
  mp.f <- metapred(d, strata = "X4", recal.int = FALSE)
  expect_true(is.numeric(p <- predict(object = mp.f, newdata = td, recal.int = TRUE)))
  expect_true(all(p <= 1))
  expect_true(all(p >= 0))

  mp.f <- metapred(d3, strata = "X4", recal.int = FALSE)
  expect_true(is.numeric(p <- predict(object = mp.f, newdata = td, recal.int = TRUE)))
  expect_true(all(p <= 1))
  expect_true(all(p >= 0))

  expect_true(is.numeric(p2 <- predict(object = mp.f, newdata = td, recal.int = FALSE)))
  expect_true(all(p2 <= 1))
  expect_true(all(p2 >= 0))

  expect_false(all(p == p2))
  expect_true(all(p != p2)) # With this seed.
})


######################

sampleBinary <- function(n = 50, J = 1, b = rep(log(2), J), alpha = NULL, col.names = NULL ) {
  J <- length(b)
  if (is.null(alpha)) alpha <- -log(sqrt(prod(exp(b))))
  if (is.null(col.names)) col.names <- c("Y", paste("X", 1:J, sep = ""))
  coefs <- c(alpha, b)
  x  <- cbind(1, matrix(rbinom(n * J, size = 1, prob = .5), nrow = n, ncol = J))
  lp <- coefs %*% t(x)
  p  <- inv.logit(lp)
  y  <- stats::rbinom(length(lp), size = 1, prob = p)

  out <- data.frame(cbind(y,x[ , -1]))
  colnames(out) <- col.names
  out
}

# library(MASS)
#
# sampleBinaryMVN <- function(n, b, s = NULL) {
#   if (is.null(s)) {
#     s <- matrix(0, ncol = length(b) - 1, nrow = length(b) - 1)
#     diag(s) <- 1
#   }
#   x <- cbind(1, MASS::mvrnorm(n = n, rep(0, nrow(s)), s, empirical = TRUE))
#   lp <- x %*% b
#   p <- metamisc:::inv.logit(lp)
#   y <- stats::rbinom(n, 1, p)
#   cbind(y,data.frame(x[ , -1]))
# }

#####
set.seed(2039420)
n <- 1000
d1 <- sampleBinary(n, b = c(.3, .4, -.5, -.6), alpha = .2)
d2 <- sampleBinary(n, b = c(.3, .4, -.5, -.6), alpha = .7)

g1 <- glm(d1, family = binomial)
g2 <- glm(d2, family = binomial)


test_that("recalibrate recalibrates glm accurately.", {
  g2$coefficients[c(2,4)] <- g1$coefficients[c(2,4)]
  expect_true(!isTRUE(all.equal(coef(g1), coef(g2)))) # As they are different models

  g2.recal.int <- metamisc:::recalibrate(g2, newdata = d1)
  expect_true(!isTRUE(all.equal(coef(g1), coef(g2.recal.int)))) # As only intercept should be updated
  expect_true(!isTRUE(all.equal(coef(g1)[1], coef(g2.recal.int)[1]))) # Due to different coefs, intercept should still be different.

  g2.recal.all <- metamisc:::recalibrate(g2, newdata = d1, f = formula(d1))
  expect_true(all.equal(coef(g1), coef(g2.recal.all))) # NOW, all have been updated.

  g2$coefficients[-1] <- g1$coefficients[-1]
  g2.recal.int2 <- metamisc:::recalibrate(g2, newdata = d1)
  expect_true(all.equal(coef(g1), coef(g2.recal.int2))) # Now, only the intercept needs to be updated.
})

