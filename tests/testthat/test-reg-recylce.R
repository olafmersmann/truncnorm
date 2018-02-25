context("reg-recycle")
test_that("recylcing 1", {
  means <- c(-1, 0, 1)
  r1 <- dtruncnorm(1, 0, Inf, means, 1)
  expect_equal(r1[1], dtruncnorm(1, 0, Inf, means[1], 1))
  expect_equal(r1[2], dtruncnorm(1, 0, Inf, means[2], 1))
  expect_equal(r1[3], dtruncnorm(1, 0, Inf, means[3], 1))
})
