test_that("sen.ci.mar validates p and returns expected structure", {
  set.seed(123)
  n <- 80
  T <- abs(rnorm(n))
  A <- abs(rnorm(n))
  score <- 0.3 * T + 0.3 * A + rnorm(n, sd = 1)
  D <- as.logical(score > stats::quantile(score, 0.7))
  D[sample(n, 20)] <- NA

  # invalid p should error
  expect_error(sen.ci.mar(T, D, A, p = -0.1, n.boot = 5, plot = FALSE))
  expect_error(sen.ci.mar(T, D, A, p =  1.1, n.boot = 5, plot = FALSE))

  out <- sen.ci.mar(T, D, A, p = 0.9, n.boot = 10, plot = FALSE)

  expect_type(out, "list")
  expect_true(all(c("n.total", "n.case", "n.control", "p.missing", "pt.est") %in% names(out)))

  expect_equal(out$n.total, n)
  expect_true(is.numeric(out$pt.est))
  expect_length(out$pt.est, 4)
})
