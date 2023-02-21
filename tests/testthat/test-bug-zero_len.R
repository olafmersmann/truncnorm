context("bug-zero_len")

test_that("dtruncnorm", {
  n = numeric(0)
  expect_equal(dtruncnorm(n, 0, Inf, 0, 1), NULL)
  expect_equal(dtruncnorm(0, n, Inf, 0, 1), NULL)
  expect_equal(dtruncnorm(0, 0, n, 0, 1), NULL)
  expect_equal(dtruncnorm(0, 0, Inf, n, 1), NULL)
  expect_equal(dtruncnorm(0, 0, Inf, 0, n), NULL)
 
  n = NULL
  expect_equal(dtruncnorm(n, 0, Inf, 0, 1), NULL)
  expect_equal(dtruncnorm(0, n, Inf, 0, 1), NULL)
  expect_equal(dtruncnorm(0, 0, n, 0, 1), NULL)
  expect_equal(dtruncnorm(0, 0, Inf, n, 1), NULL)
  expect_equal(dtruncnorm(0, 0, Inf, 0, n), NULL)
})

test_that("ptruncnorm", {
  n = numeric(0)
  expect_equal(qtruncnorm(n, 0, Inf, 0, 1), NULL)
  expect_equal(qtruncnorm(0, n, Inf, 0, 1), NULL)
  expect_equal(qtruncnorm(0, 0, n, 0, 1), NULL)
  expect_equal(qtruncnorm(0, 0, Inf, n, 1), NULL)
  expect_equal(qtruncnorm(0, 0, Inf, 0, n), NULL)
 
  n = NULL
  expect_equal(qtruncnorm(n, 0, Inf, 0, 1), NULL)
  expect_equal(qtruncnorm(0, n, Inf, 0, 1), NULL)
  expect_equal(qtruncnorm(0, 0, n, 0, 1), NULL)
  expect_equal(qtruncnorm(0, 0, Inf, n, 1), NULL)
  expect_equal(qtruncnorm(0, 0, Inf, 0, n), NULL)
})

test_that("qtruncnorm", {
  n = numeric(0)
  expect_equal(qtruncnorm(n, 0, Inf, 0, 1), NULL)
  expect_equal(qtruncnorm(0, n, Inf, 0, 1), NULL)
  expect_equal(qtruncnorm(0, 0, n, 0, 1), NULL)
  expect_equal(qtruncnorm(0, 0, Inf, n, 1), NULL)
  expect_equal(qtruncnorm(0, 0, Inf, 0, n), NULL)
 
  n = NULL
  expect_equal(qtruncnorm(n, 0, Inf, 0, 1), NULL)
  expect_equal(qtruncnorm(0, n, Inf, 0, 1), NULL)
  expect_equal(qtruncnorm(0, 0, n, 0, 1), NULL)
  expect_equal(qtruncnorm(0, 0, Inf, n, 1), NULL)
  expect_equal(qtruncnorm(0, 0, Inf, 0, n), NULL)
})

test_that("etruncnorm", {
  n = numeric(0)
  expect_equal(etruncnorm(n, Inf, 0, 1), NULL)
  expect_equal(etruncnorm(0, n, 0, 1), NULL)
  expect_equal(etruncnorm(0, Inf, n, 1), NULL)
  expect_equal(etruncnorm(0, Inf, 0, n), NULL)
 
  n = NULL
  expect_equal(etruncnorm(n, Inf, 0, 1), NULL)
  expect_equal(etruncnorm(0, n, 0, 1), NULL)
  expect_equal(etruncnorm(0, Inf, n, 1), NULL)
  expect_equal(etruncnorm(0, Inf, 0, n), NULL)
})

test_that("vtruncnorm", {
  n = numeric(0)
  expect_equal(vtruncnorm(n, Inf, 0, 1), NULL)
  expect_equal(vtruncnorm(0, n, 0, 1), NULL)
  expect_equal(vtruncnorm(0, Inf, n, 1), NULL)
  expect_equal(vtruncnorm(0, Inf, 0, n), NULL)
 
  n = NULL
  expect_equal(vtruncnorm(n, Inf, 0, 1), NULL)
  expect_equal(vtruncnorm(0, n, 0, 1), NULL)
  expect_equal(vtruncnorm(0, Inf, n, 1), NULL)
  expect_equal(vtruncnorm(0, Inf, 0, n), NULL)
})
