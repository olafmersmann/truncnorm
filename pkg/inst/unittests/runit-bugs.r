##
## runit-bugs.r - Unit tests for reported bugs to avoid regressions
##

## Cast n argument from numeric to integer.
test.001 <- function() {
  data <- rtruncnorm(10, 0, 1, 0, 1)
  checkEquals(length(data), 10L)
  checkException(rtruncnorm("a", 0, 1, 0, 1))
}

## If n is a vector, return length(n) elements
test.002 <- function() {
  data <- rtruncnorm(1:10, 0, 1, 0, 1)
  checkEquals(length(data), length(1:10))
}

## These two regression tests are not checked by default because the
## are stochastic in nature and might fail even though everything is
## OK. If you change anything in rtruncnorm.c, be sure to run these
## and investigate if they fail.
if (FALSE) {
  ## Wrong sampling strategy:
  test.003 <- function() {
    x <- rtruncnorm(n=10000000, -0.75, 0.75, 0, 1)
    mx <- mean(x)
    vx <- var(x)
    ev <- etruncnorm(-0.75, 0.75, 0, 1)
    vv <- vtruncnorm(-0.75, 0.75, 0, 1)
    checkEqualsNumeric(ev, mx, tolerance=1e-3)
    checkEqualsNumeric(vv, vx, tolerance=1e-3)
  }
  
  ## Wrong sampling strategy:
  test.004 <- function() {
    x <- rtruncnorm(n=10000000, -1, 4, 0, 1)
    mx <- mean(x)
    vx <- var(x)
    ev <- etruncnorm(-1, 4, 0, 1)
    vv <- vtruncnorm(-1, 4, 0, 1)
    checkEqualsNumeric(ev, mx, tolerance=1e-3)
    checkEqualsNumeric(vv, vx, tolerance=1e-3)
  }
}
