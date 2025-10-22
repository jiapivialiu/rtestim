# ==============================================================================
# Tests for linear_single_solve_test() function
# Tests the linear system solve step with different solvers
# ==============================================================================

test_that("linear_single_solve_test works with QR solver", {
  # Create simple test data
  set.seed(42)
  n <- 50
  k <- 2  # korder
  x <- 1:n
  
  # Generate test data
  y <- rnorm(n, mean = 5, sd = 1)
  weights <- rep(1, n) 
  rho <- 1.0
  adj_mean <- rnorm(n - k, mean = 0, sd = 0.1)
  
  # Test QR solver (linear_solver = 1)
  result_qr <- linear_single_solve_test(
    linear_solver = 1,
    y = y,
    weights = weights,
    x = x,
    rho = rho,
    adj_mean = adj_mean
  )
  
  # Check result properties
  expect_type(result_qr, "double")
  expect_length(result_qr, n)
  
  # Check for no NAs or Infs
  expect_false(any(is.na(result_qr)),
               label = "QR solver should not produce NAs")
  expect_true(all(is.finite(result_qr)),
              label = "QR solver should not produce infinite values")
})

test_that("linear_single_solve_test works with KF solver", {
  # Create simple test data
  set.seed(42)
  n <- 50
  k <- 2  # korder
  x <- 1:n
  
  # Generate test data
  y <- rnorm(n, mean = 5, sd = 1)
  weights <- rep(1, n)
  rho <- 1.0
  adj_mean <- rnorm(n - k, mean = 0, sd = 0.1)
  
  # Test KF solver (linear_solver = 2)
  result_kf <- linear_single_solve_test(
    linear_solver = 2,
    y = y,
    weights = weights,
    x = x,
    rho = rho,
    adj_mean = adj_mean
  )
  
  # Check result properties
  expect_type(result_kf, "double")
  expect_length(result_kf, n)
  
  # Check for no NAs or Infs
  expect_false(any(is.na(result_kf)),
               label = "KF solver should not produce NAs")
  expect_true(all(is.finite(result_kf)),
              label = "KF solver should not produce infinite values")
})

test_that("QR and KF solvers produce similar results in linear_single_solve_test", {
  # Create simple test data
  set.seed(42)
  n <- 50
  k <- 2
  x <- 1:n
  
  # Generate test data
  y <- rnorm(n, mean = 5, sd = 1)
  weights <- rep(1, n)
  rho <- 1.0
  adj_mean <- rnorm(n - k, mean = 0, sd = 0.1)
  
  # Test both solvers
  result_qr <- linear_single_solve_test(1, y, weights, x, rho, adj_mean)
  result_kf <- linear_single_solve_test(2, y, weights, x, rho, adj_mean)
  
  # Both should be valid
  expect_false(any(is.na(result_qr)))
  expect_false(any(is.na(result_kf)))
  expect_true(all(is.finite(result_qr)))
  expect_true(all(is.finite(result_kf)))
  
  # Solvers should produce correlated results
  cor_val <- cor(result_qr, result_kf)
  expect_gt(cor_val, 0.9,
            label = "QR and KF solvers should produce correlated results (r > 0.9)")
})

test_that("linear_single_solve_test handles different korder values", {
  set.seed(123)
  n <- 60
  x <- 1:n
  y <- rnorm(n, mean = 10, sd = 2)
  weights <- rep(1, n)
  rho <- 1.0
  
  # Test with korder = 1
  k1 <- 1
  adj_mean1 <- rnorm(n - k1)
  result1 <- linear_single_solve_test(1, y, weights, x, rho, adj_mean1)
  expect_length(result1, n)
  expect_false(any(is.na(result1)))
  
  # Test with korder = 2
  k2 <- 2
  adj_mean2 <- rnorm(n - k2)
  result2 <- linear_single_solve_test(1, y, weights, x, rho, adj_mean2)
  expect_length(result2, n)
  expect_false(any(is.na(result2)))
  
  # Test with korder = 3
  k3 <- 3
  adj_mean3 <- rnorm(n - k3)
  result3 <- linear_single_solve_test(1, y, weights, x, rho, adj_mean3)
  expect_length(result3, n)
  expect_false(any(is.na(result3)))
})

test_that("linear_single_solve_test handles varying weights", {
  set.seed(456)
  n <- 40
  k <- 2
  x <- 1:n
  y <- rnorm(n, mean = 5, sd = 1)
  rho <- 1.0
  adj_mean <- rnorm(n - k)
  
  # Test with uniform weights
  weights_uniform <- rep(1, n)
  result_uniform <- linear_single_solve_test(1, y, weights_uniform, x, rho, adj_mean)
  expect_false(any(is.na(result_uniform)))
  
  # Test with varying weights
  weights_varying <- runif(n, min = 0.5, max = 2.0)
  result_varying <- linear_single_solve_test(1, y, weights_varying, x, rho, adj_mean)
  expect_false(any(is.na(result_varying)))
  
  # Results should be different with different weights
  expect_false(identical(result_uniform, result_varying))
})

test_that("linear_single_solve_test handles different rho values", {
  set.seed(789)
  n <- 40
  k <- 2
  x <- 1:n
  y <- rnorm(n, mean = 5, sd = 1)
  weights <- rep(1, n)
  adj_mean <- rnorm(n - k)
  
  # Test with small rho
  result_small_rho <- linear_single_solve_test(1, y, weights, x, 0.1, adj_mean)
  expect_false(any(is.na(result_small_rho)))
  
  # Test with medium rho
  result_med_rho <- linear_single_solve_test(1, y, weights, x, 1.0, adj_mean)
  expect_false(any(is.na(result_med_rho)))
  
  # Test with large rho
  result_large_rho <- linear_single_solve_test(1, y, weights, x, 10.0, adj_mean)
  expect_false(any(is.na(result_large_rho)))
  
  # Results should differ with different rho
  expect_false(identical(result_small_rho, result_large_rho))
})

test_that("linear_single_solve_test handles unequally spaced x", {
  set.seed(321)
  n <- 50
  k <- 2
  
  # Create unequally spaced x
  x <- cumsum(runif(n, min = 0.5, max = 2.0))
  
  y <- rnorm(n, mean = 5, sd = 1)
  weights <- rep(1, n)
  rho <- 1.0
  adj_mean <- rnorm(n - k)
  
  # Test with QR solver
  result <- linear_single_solve_test(1, y, weights, x, rho, adj_mean)
  
  expect_length(result, n)
  expect_false(any(is.na(result)))
  expect_true(all(is.finite(result)))
})

test_that("linear_single_solve_test handles edge cases", {
  # Small dataset
  n <- 20
  k <- 2
  x <- 1:n
  y <- rnorm(n, mean = 1, sd = 0.1)
  weights <- rep(1, n)
  rho <- 1.0
  adj_mean <- rnorm(n - k)
  
  expect_no_error(
    result <- linear_single_solve_test(1, y, weights, x, rho, adj_mean)
  )
  expect_false(any(is.na(result)))
  
  # Very small values
  y_small <- rep(0.01, n)
  expect_no_error(
    result_small <- linear_single_solve_test(1, y_small, weights, x, rho, adj_mean)
  )
  expect_false(any(is.na(result_small)))
  
  # Large values
  y_large <- rep(1000, n)
  expect_no_error(
    result_large <- linear_single_solve_test(1, y_large, weights, x, rho, adj_mean)
  )
  expect_false(any(is.na(result_large)))
})

test_that("linear_single_solve_test produces consistent results", {
  # Test that same inputs give same outputs (deterministic)
  set.seed(999)
  n <- 40
  k <- 2
  x <- 1:n
  y <- rnorm(n)
  weights <- rep(1, n)
  rho <- 1.0
  adj_mean <- rnorm(n - k)
  
  result1 <- linear_single_solve_test(1, y, weights, x, rho, adj_mean)
  result2 <- linear_single_solve_test(1, y, weights, x, rho, adj_mean)
  
  expect_equal(result1, result2,
               label = "Function should be deterministic")
})

test_that("linear_single_solve_test with correct dimensions runs successfully", {
  # Test that function works with correct dimensions
  n <- 30
  k <- 2
  x <- 1:n
  y <- rnorm(n)
  weights <- rep(1, n)
  rho <- 1.0
  adj_mean <- rnorm(n - k)
  
  # Should work with correct dimensions
  expect_no_error(
    result <- linear_single_solve_test(1, y, weights, x, rho, adj_mean)
  )
  
  expect_length(result, n)
  expect_false(any(is.na(result)))
  expect_true(all(is.finite(result)))
})
