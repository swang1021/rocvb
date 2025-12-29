test_that("yi.ci.mar returns expected structure", {
  set.seed(123)
  n <- 80
  T <- abs(rnorm(n))
  A <- abs(rnorm(n))
  score <- 0.3 * T + 0.3 * A + rnorm(n, sd = 1)
  D <- as.logical(score > stats::quantile(score, 0.7))
  D[sample(n, 20)] <- NA

  out <- yi.ci.mar(T, D, A, n.boot = 10, plot = FALSE)

  expect_type(out, "list")
  expect_true(all(c("n.total", "n.case", "n.control", "p.missing", "pt.est", "optimal.cutoff") %in% names(out)))

  expect_equal(out$n.total, n)
  expect_true(is.numeric(out$pt.est))
  expect_length(out$pt.est, 4)

  expect_true(is.numeric(out$optimal.cutoff) || is.matrix(out$optimal.cutoff) || is.data.frame(out$optimal.cutoff))
})
