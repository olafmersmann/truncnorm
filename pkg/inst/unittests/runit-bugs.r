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
