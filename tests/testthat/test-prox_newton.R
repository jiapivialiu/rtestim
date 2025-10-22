# ==============================================================================
# Tests: prox_newton_testing() C++ function
# Tests the proximal Newton method for Poisson trend filtering
# ==============================================================================

test_that("prox_newton_testing with QR solver produces valid output structure", {
  # Create synthetic Poisson data
  set.seed(123)
  n <- 50
  x <- 1:n
  korder <- 2

  # Generate counts with smooth trend
  true_theta <- 100 + 30 * sin(2 * pi * x / 20)
  y <- rpois(n, lambda = true_theta)
  w <- rep(1, n)

  # Test parameters
  lambda <- 10.0
  M <- 100
  Minner <- 50
  Mline <- 20
  tol <- 1e-6

  result <- prox_newton_testing(
    M = M, Minner = Minner, Mline = Mline, korder = korder,
    y = y, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1,  # QR solver
    tol = tol
  )

  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("lambda", "theta", "niter"))

  # Check lambda is preserved
  expect_equal(result$lambda, lambda)

  # Check theta dimensions and positivity
  expect_length(result$theta, n)
  expect_true(all(result$theta > 0))
  expect_true(all(is.finite(result$theta)))

  # Check iterations
  expect_type(result$niter, "integer")
  expect_true(result$niter >= 0)
  expect_true(result$niter <= M * Minner)  # Total possible iterations
})

test_that("prox_newton_testing with different korder values works correctly", {
  set.seed(456)
  n <- 40
  x <- 1:n
  y <- rpois(n, lambda = 50 + 20 * (x / n))
  w <- rep(1, n)
  lambda <- 5.0

  # Test korder = 1 (linear trend)
  result_k1 <- prox_newton_testing(
    M = 50, Minner = 30, Mline = 10, korder = 1,
    y = y, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )
  expect_length(result_k1$theta, n)
  expect_true(all(is.finite(result_k1$theta)))

  # Test korder = 2 (quadratic trend)
  result_k2 <- prox_newton_testing(
    M = 50, Minner = 30, Mline = 10, korder = 2,
    y = y, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )
  expect_length(result_k2$theta, n)
  expect_true(all(is.finite(result_k2$theta)))

  # Test korder = 3 (cubic trend)
  result_k3 <- prox_newton_testing(
    M = 50, Minner = 30, Mline = 10, korder = 3,
    y = y, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )
  expect_length(result_k3$theta, n)
  expect_true(all(is.finite(result_k3$theta)))

  # Higher korder should generally allow more flexibility
  # (though exact relationship depends on data and lambda)
  expect_true(all(result_k3$theta > 0))
})

test_that("prox_newton_testing with different lambda values produces expected smoothness", {
  set.seed(789)
  n <- 60
  x <- 1:n
  # Data with noise
  true_signal <- 80 + 40 * sin(2 * pi * x / 30)
  y <- rpois(n, lambda = true_signal)
  w <- rep(1, n)
  korder <- 2

  # Small lambda (less smoothing)
  result_small <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y, x = x, w = w, lambda = 1.0,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )

  # Large lambda (more smoothing)
  result_large <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y, x = x, w = w, lambda = 100.0,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )

  # Calculate smoothness (sum of squared differences)
  smoothness_small <- sum(diff(result_small$theta, differences = korder)^2)
  smoothness_large <- sum(diff(result_large$theta, differences = korder)^2)

  # Larger lambda should produce smoother result
  expect_true(smoothness_large < smoothness_small)

  # Both should be valid
  expect_true(all(is.finite(result_small$theta)))
  expect_true(all(is.finite(result_large$theta)))
})

test_that("prox_newton_testing with weighted data works correctly", {
  set.seed(321)
  n <- 50
  x <- 1:n
  korder <- 2
  y <- rpois(n, lambda = 60 + 20 * cos(2 * pi * x / 25))

  # Uniform weights
  w_uniform <- rep(1, n)
  result_uniform <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y, x = x, w = w_uniform, lambda = 10.0,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )

  # Varying weights (emphasize middle)
  w_varying <- dnorm(x, mean = n/2, sd = n/4)
  w_varying <- w_varying / sum(w_varying) * n  # Normalize
  result_varying <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y, x = x, w = w_varying, lambda = 10.0,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )

  # Results should differ
  expect_false(identical(result_uniform$theta, result_varying$theta))

  # Both should be valid
  expect_true(all(is.finite(result_uniform$theta)))
  expect_true(all(is.finite(result_varying$theta)))
  expect_length(result_uniform$theta, n)
  expect_length(result_varying$theta, n)
})

test_that("prox_newton_testing with unequally spaced x works correctly", {
  set.seed(654)
  n <- 40
  # Unequally spaced observation times
  x <- cumsum(runif(n, 0.5, 1.5))
  korder <- 2

  true_signal <- 70 + 25 * sin(2 * pi * x / max(x))
  y <- rpois(n, lambda = true_signal)
  w <- rep(1, n)

  result <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y, x = x, w = w, lambda = 10.0,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )

  # Check valid output
  expect_length(result$theta, n)
  expect_true(all(result$theta > 0))
  expect_true(all(is.finite(result$theta)))
  expect_type(result$niter, "integer")
})

test_that("prox_newton_testing with different line search parameters works", {
  set.seed(987)
  n <- 50
  x <- 1:n
  korder <- 2
  y <- rpois(n, lambda = 90 + 30 * (x / n)^2)
  w <- rep(1, n)
  lambda <- 15.0

  # Conservative line search (smaller alpha, larger gamma)
  result_conservative <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y, x = x, w = w, lambda = lambda,
    ls_alpha = 0.3, ls_gamma = 0.9,
    linear_solver = 1, tol = 1e-6
  )

  # Aggressive line search (larger alpha, smaller gamma)
  result_aggressive <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y, x = x, w = w, lambda = lambda,
    ls_alpha = 0.7, ls_gamma = 0.7,
    linear_solver = 1, tol = 1e-6
  )

  # Both should produce valid results
  expect_true(all(is.finite(result_conservative$theta)))
  expect_true(all(is.finite(result_aggressive$theta)))

  # Iteration counts might differ
  expect_type(result_conservative$niter, "integer")
  expect_type(result_aggressive$niter, "integer")
})

test_that("prox_newton_testing convergence with tight tolerance", {
  set.seed(111)
  n <- 30
  x <- 1:n
  korder <- 2
  y <- rpois(n, lambda = 50)  # Constant signal with noise
  w <- rep(1, n)
  lambda <- 20.0

  # Tight tolerance
  result_tight <- prox_newton_testing(
    M = 200, Minner = 100, Mline = 30, korder = korder,
    y = y, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-8
  )

  # Loose tolerance
  result_loose <- prox_newton_testing(
    M = 200, Minner = 100, Mline = 30, korder = korder,
    y = y, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-4
  )

  # Tight tolerance might require more iterations
  # (though not guaranteed depending on convergence behavior)
  expect_true(result_tight$niter >= 0)
  expect_true(result_loose$niter >= 0)

  # Both should be close
  expect_true(cor(result_tight$theta, result_loose$theta) > 0.99)
})

test_that("prox_newton_testing handles edge cases correctly", {
  set.seed(222)
  n <- 50
  x <- 1:n
  korder <- 2
  w <- rep(1, n)
  lambda <- 10.0

  # Case 1: Very small counts
  y_small <- rpois(n, lambda = 2)
  result_small <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y_small, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )
  expect_true(all(is.finite(result_small$theta)))
  expect_true(all(result_small$theta > 0))

  # Case 2: Large counts
  y_large <- rpois(n, lambda = 100)
  result_large <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y_large, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )
  expect_true(all(is.finite(result_large$theta)))
  expect_true(all(result_large$theta > 0))

  # Case 3: Some zeros in data (valid for Poisson)
  y_zeros <- c(0, 0, rpois(n-4, lambda = 50), 0, 0)
  result_zeros <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y_zeros, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 1, tol = 1e-6
  )
  expect_true(all(is.finite(result_zeros$theta)))
  expect_true(all(result_zeros$theta > 0))
})

test_that("prox_newton_testing with KF solver produces valid results", {
  set.seed(333)
  n <- 50
  x <- 1:n
  korder <- 2
  y <- rpois(n, lambda = 80 + 30 * sin(2 * pi * x / 20))
  w <- rep(1, n)
  lambda <- 10.0

  # Test KF solver
  result_kf <- prox_newton_testing(
    M = 100, Minner = 50, Mline = 20, korder = korder,
    y = y, x = x, w = w, lambda = lambda,
    ls_alpha = 0.5, ls_gamma = 0.8,
    linear_solver = 2,  # KF solver
    tol = 1e-6
  )

  # Check output structure
  expect_type(result_kf, "list")
  expect_named(result_kf, c("lambda", "theta", "niter"))
  expect_length(result_kf$theta, n)
  expect_true(all(is.finite(result_kf$theta)))
  expect_true(all(result_kf$theta > 0))
})

