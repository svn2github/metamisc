context("utility functions for metapred")


### Some stuff necessary for testing
n <- 2
td <- data.frame(y = rep(0, n * 4), x = rep(0, n * 4), z = c(rep(0, n * 2), rep(1, n * 2)), s = rep(c(rep(0, n), rep(1, n)), 2))

test_that("Centering works.", {
  expect_true(is.data.frame(cd <- metamisc:::centerData(data = td, cluster.vec = td$s, center.which = FALSE) ) )
  expect_identical(td, cd) # Nothing should change if center.which = FALSE
  
  expect_true(is.data.frame(cd <- metamisc:::centerData(data = td, cluster.vec = td$s, center.which = seq_len(ncol(td)))))
  expect_identical(cd$x, rep(0, n * 4))
  expect_identical(cd$z, c(rep(-.5, n * 2), rep(.5, n * 2)))
  expect_identical(cd$y, rep(0, n * 4))
  # Note that the center indicator is also centered. This is not intended, but can be worked around.
  
  expect_true(is.data.frame(cd <- metamisc:::centerData(data = td, cluster.vec = td$s))) # used to be:center.1st = FALSE, center.rest = TRUE (is now the default)
  expect_identical(cd$x, rep(0, n * 4))
  expect_identical(cd$z, c(rep(-.5, n * 2), rep(.5, n * 2)))
  expect_identical(cd$y, td$y)
  
  expect_true(is.data.frame(cd <- metamisc:::centerData(data = td, cluster.vec = td$s, center.which = 0)))
  expect_identical(cd$x, td$x)
  expect_identical(cd$z, td$z)
  expect_identical(cd$y, td$y)
  
  
  # Or using the cluster variable in the data.frame:
  expect_true(is.data.frame(cd <- metamisc:::centerData(data = td, cluster.var = "s", center.which = 0)))
  expect_identical(cd$x, td$x)
  expect_identical(cd$z, td$z)
  expect_identical(cd$y, td$y)
})


test_that("asDataList and Reduce are complements.", {
  expect_true(is.list(dl <- metamisc:::asDataList(td, td$z)))
  expect_identical(td, Reduce(rbind, dl)) # Note that this is not always true. But with these parameters it should.
})

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