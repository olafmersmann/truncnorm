context("reg-recycle")

test_that("recylcing x in dtruncnorm", {
  means <- c(-1, 0, 1)
  x <- 1
  r1 <- dtruncnorm(x, 0, Inf, means, 1)
  expect_equal(r1[1], dtruncnorm(x, 0, Inf, means[1], 1))
  expect_equal(r1[2], dtruncnorm(x, 0, Inf, means[2], 1))
  expect_equal(r1[3], dtruncnorm(x, 0, Inf, means[3], 1))
})

test_that("recylcing p in ptruncnorm", {
  means <- c(-1, 0, 1)
  q <- 0.5
  r1 <- ptruncnorm(q, 0, Inf, means, 1)
  expect_equal(r1[1], ptruncnorm(q, 0, Inf, means[1], 1))
  expect_equal(r1[2], ptruncnorm(q, 0, Inf, means[2], 1))
  expect_equal(r1[3], ptruncnorm(q, 0, Inf, means[3], 1))
})

test_that("recylcing p in qtruncnorm", {
  means <- c(-1, 0, 1)
  p <- 0.5
  r1 <- qtruncnorm(p, 0, Inf, means, 1)
  expect_equal(r1[1], qtruncnorm(p, 0, Inf, means[1], 1))
  expect_equal(r1[2], qtruncnorm(p, 0, Inf, means[2], 1))
  expect_equal(r1[3], qtruncnorm(p, 0, Inf, means[3], 1))
})
