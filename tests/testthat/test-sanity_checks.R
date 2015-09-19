context("sanity checks")

################################################################################
## Check d/e/vtruncnorm all in one function:
check_dev <- function(a, b, mean=0, sd=1) {
  prefix <- sprintf("DEV: a=%f, b=%f, mean=%f, sd=%f", a, b, mean, sd)
  e <- etruncnorm(a, b, mean, sd)
  v <- vtruncnorm(a, b, mean, sd)

  id <- integrate(function(x) dtruncnorm(x, a, b, mean, sd), a, b)$value
  ee <- integrate(function(x) x * dtruncnorm(x, a, b, mean, sd), a, b)$value
  ev <- integrate(function(x) (x-ee)^2 * dtruncnorm(x, a, b, mean, sd), a, b)$value

  test_that(prefix, {
    expect_equal(id, 1.0, tolerance=0.00005)
    expect_equal(e, ee, tolerance=0.00005)
    expect_equal(v, ev, tolerance=0.00005)
  })
}

## Left truncated:
check_dev(-3, Inf, 0, 1)
check_dev(-2, Inf, 1, 1)
check_dev( 2, Inf, 0, 1)
check_dev( 3, Inf, 1, 1)

check_dev(-3, Inf, 0, 2)
check_dev(-2, Inf, 1, 2)
check_dev( 2, Inf, 0, 2)
check_dev( 3, Inf, 1, 2)

## Doubly truncated:
check_dev(-3.0, -2.5, 0, 1)
check_dev(-3.0, -1.5, 0, 1)
check_dev(-3.0, -0.5, 0, 1)
check_dev(-3.0,  0.5, 0, 1)

check_dev(0.0, 0.5, 0, 1)
check_dev(0.0, 1.5, 0, 1)
check_dev(0.0, 2.5, 0, 1)
check_dev(0.0, 3.5, 0, 1)

## Extreme cases:
check_dev( 0.0, 1.0,  0.0, 10.0)
check_dev( 0.0, 1.0,  5.0,  1.0)
check_dev(-1.0, 0.0,  0.0, 10.0)
check_dev( 0.0, 1.0, -5.0,  1.0)
check_dev( 0.0, 1.0,  5.0,  0.1)

## Integer arguments:
check_dev(0L, 1L, 0L, 10L)

################################################################################
## Sanity checks on random number generators
check_r <- function(a, b, mean, sd, n=10000) {
  prefix <- sprintf("R: a=%f, b=%f, mean=%f, sd=%f", a, b, mean, sd)
  x <- rtruncnorm(n, a, b, mean, sd)
  e.x <- mean(x)
  e <- etruncnorm(a, b, mean, sd)
  true_sd <- sqrt(vtruncnorm(a, b, mean, sd)) 

  ## FIXME: Really sample from open intervall?
  test_that(prefix, {
    expect_true(all(x > a))
    expect_true(all(x < b))
    expect_equal(mean(x), e, tolerance=0.05, scale=sd)
    expect_equal(sd(x), true_sd, tolerance=0.05, scale=sd)
  })
}

## rtruncnorm == rnorm:
check_r(-Inf, Inf, 0, 1)

## 0 in (a, b):
check_r(-1, 1, 0, 1)
check_r(-1, 1, 1, 1)
check_r(-1, 1, 0, 2)

## 0 < (a, b):
check_r(1, 2, 0, 1)
check_r(1, 2, 1, 1)
check_r(1, 2, 0, 2)

## 0 > (a, b):
check_r(-2, -1, 0, 1)
check_r(-2, -1, 1, 1)
check_r(-2, -1, 0, 2)

## left truncation:
check_r(-2, Inf, 0, 1)
check_r(-2, Inf, 1, 1)
check_r(-2, Inf, 0, 2)
check_r( 0, Inf, 0, 1)
check_r( 0, Inf, 1, 1)
check_r( 0, Inf, 0, 2)
check_r( 2, Inf, 0, 1)
check_r( 2, Inf, 1, 1)
check_r( 2, Inf, 0, 2)

check_r(-0.2, Inf, 0, 1)
check_r(-0.2, Inf, 1, 1)
check_r(-0.2, Inf, 0, 2)
check_r( 0.0, Inf, 0, 1)
check_r( 0.0, Inf, 1, 1)
check_r( 0.0, Inf, 0, 2)
check_r( 0.2, Inf, 0, 1)
check_r( 0.2, Inf, 1, 1)
check_r( 0.2, Inf, 0, 2)

## Right truncation:
check_r(-Inf, -2, 0, 1)
check_r(-Inf, -2, 1, 1)
check_r(-Inf, -2, 0, 2)
check_r(-Inf,  0, 0, 1)
check_r(-Inf,  0, 1, 1)
check_r(-Inf,  0, 0, 2)
check_r(-Inf,  2, 0, 1)
check_r(-Inf,  2, 1, 1)
check_r(-Inf,  2, 0, 2)

check_r(-Inf, -0.2, 0, 1)
check_r(-Inf, -0.2, 1, 1)
check_r(-Inf, -0.2, 0, 2)
check_r(-Inf,  0.0, 0, 1)
check_r(-Inf,  0.0, 1, 1)
check_r(-Inf,  0.0, 0, 2)
check_r(-Inf,  0.2, 0, 1)
check_r(-Inf,  0.2, 1, 1)
check_r(-Inf,  0.2, 0, 2)

## Extreme examples:
check_r(-5, -4, 0, 1)

## Integer examples:
check_r(-5L, -4L, 0L, 1L)

################################################################################
check_pq <- function(a, b, mean, sd) {
  prefix <- sprintf("PQ: a=%f, b=%f, mean=%f, sd=%f", a, b, mean, sd)
  test_that(prefix, {
    for (p in runif(500)) {
      q <- qtruncnorm(p, a, b, mean, sd)
      pp <- ptruncnorm(q, a, b, mean, sd)
      expect_equal(pp, p, tolerance=0.00001)
    }
  })
}

check_pq(-1, 0, 0, 1)
check_pq(-1, 1, 0, 1)
check_pq( 1, 2, 0, 1)
check_pq(-1, 0, 4, 1)
check_pq(-1, 1, 4, 1)
check_pq( 1, 2, 4, 1)
check_pq(-1, 0, 0, 3)
check_pq(-1, 1, 0, 3)
check_pq( 1, 2, 0, 3)
check_pq(-1, Inf, 0, 1)
check_pq(-1, Inf, 4, 1)
check_pq(-1, Inf, 0, 3)
check_pq(-Inf, 1, 0, 1)
check_pq(-Inf, 1, 4, 1)
check_pq(-Inf, 1, 0, 3)

## Integer examples:
check_pq(1L, 2L, 0L, 3L)
