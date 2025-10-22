# ==============================================================================
# Tests: admm_gauss_testing() C++ function
# Tests the Gaussian ADMM solver step directly
# This helps diagnose differences between QR and KF solvers
# ==============================================================================

test_that("admm_gauss_testing with QR solver produces valid theta output", {
  # Create synthetic test data
  set.seed(123)
  n <- 50
  x <- 1:n
  korder <- 2
  
  # Generate smooth signals (centered data for Gaussian ADMM)
  y <- 20 * sin(2 * pi * x / 20) + rnorm(n, 0, 2)
  
  # Initialize variables
  w <- rep(1, n)
  theta <- log(abs(y) + 1) + rnorm(n, 0, 0.1)
  z <- rep(0, n - korder)
  u <- rep(0, n - korder)
  rho <- 1.0
  lam_z <- 0.5
  
  # Test QR solver
  result_qr <- admm_gauss_testing(
    M = 100,
    korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 1,  # QR
    tol = 1e-6
  )
  
  # Check structure
  expect_type(result_qr, "list")
  expect_named(result_qr, c("theta", "z", "u", "theta_init", "z_init", "u_init",
                            "niter", "converged", "num_na_theta", "num_inf_theta",
                            "num_na_z", "num_inf_z", "max_theta_change", 
                            "mean_theta_change", "solver"))
  
  # Check solver type
  expect_equal(result_qr$solver, "QR")
  
  # Check for valid results (no NAs or Infs)
  expect_equal(result_qr$num_na_theta, 0)
  expect_equal(result_qr$num_inf_theta, 0)
  expect_equal(result_qr$num_na_z, 0)
  expect_equal(result_qr$num_inf_z, 0)
  
  # Check dimensions
  expect_length(result_qr$theta, n)
  expect_length(result_qr$z, n - korder)
  expect_length(result_qr$u, n - korder)
  
  # Check that theta actually changed
  expect_true(result_qr$max_theta_change > 0)
  expect_true(result_qr$mean_theta_change > 0)
  
  # Check convergence
  expect_true(result_qr$converged)
  expect_true(result_qr$niter < 100)
})

test_that("admm_gauss_testing with KF solver produces valid theta output", {
  # Create synthetic test data
  set.seed(123)
  n <- 50
  x <- 1:n
  korder <- 2
  
  # Generate smooth signals
  y <- 20 * sin(2 * pi * x / 20) + rnorm(n, 0, 2)
  
  # Initialize variables
  w <- rep(1, n)
  theta <- log(abs(y) + 1) + rnorm(n, 0, 0.1)
  z <- rep(0, n - korder)
  u <- rep(0, n - korder)
  rho <- 1.0
  lam_z <- 0.5
  
  # Test KF solver
  result_kf <- admm_gauss_testing(
    M = 100,
    korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 2,  # KF
    tol = 1e-6
  )
  
  # Check structure
  expect_type(result_kf, "list")
  expect_equal(result_kf$solver, "KF")
  
  # Check for valid results (THIS IS WHERE KF FAILS - produces NAs)
  expect_equal(result_kf$num_na_theta, 0, 
               info = sprintf("KF solver produced %d NA values in theta", result_kf$num_na_theta))
  expect_equal(result_kf$num_inf_theta, 0,
               info = sprintf("KF solver produced %d Inf values in theta", result_kf$num_inf_theta))
  expect_equal(result_kf$num_na_z, 0)
  expect_equal(result_kf$num_inf_z, 0)
  
  # Check dimensions
  expect_length(result_kf$theta, n)
  expect_length(result_kf$z, n - korder)
  
  # Check that theta changed
  expect_true(result_kf$max_theta_change > 0)
  expect_true(result_kf$mean_theta_change > 0)
})

test_that("admm_gauss_testing: QR and KF solvers produce similar theta", {  
  # Create synthetic test data
  set.seed(456)
  n <- 50
  x <- 1:n
  korder <- 2
  
  # Generate smooth signals
  y <- 30 * sin(2 * pi * x / 25) + rnorm(n, 0, 3)
  
  # Initialize variables (same for both)
  w <- rep(1, n)
  theta <- log(abs(y) + 1) + rnorm(n, 0, 0.1)
  z <- rep(0, n - korder)
  u <- rep(0, n - korder)
  rho <- 1.0
  lam_z <- 0.5
  
  # Test QR solver
  result_qr <- admm_gauss_testing(
    M = 100, korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  
  # Test KF solver (same initial conditions)
  result_kf <- admm_gauss_testing(
    M = 100, korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 2, tol = 1e-6
  )
  
  # Both should be valid
  expect_equal(result_qr$num_na_theta, 0)
  expect_equal(result_kf$num_na_theta, 0)
  
  # Compare theta outputs
  expect_equal(length(result_qr$theta), length(result_kf$theta))
  
  # Results should be highly correlated
  correlation <- cor(result_qr$theta, result_kf$theta)
  expect_gt(correlation, 0.85, 
            label = sprintf("QR-KF correlation: %.4f", correlation))
  
  # Relative difference should be small
  rel_diff <- abs(result_qr$theta - result_kf$theta) / (abs(result_qr$theta) + 1e-10)
  mean_rel_diff <- mean(rel_diff)
  expect_lt(mean_rel_diff, 0.05,
            label = sprintf("Mean relative difference: %.4f", mean_rel_diff))
})

test_that("admm_gauss_testing with different korder values", {
  set.seed(789)
  n <- 40
  x <- 1:n
  y <- 15 + 10 * (x / n)^2 + rnorm(n, 0, 2)
  w <- rep(1, n)
  rho <- 1.0
  lam_z <- 0.5
  
  # Test korder = 1
  theta1 <- log(abs(y) + 1)
  z1 <- rep(0, n - 1)
  u1 <- rep(0, n - 1)
  result_k1 <- admm_gauss_testing(
    M = 100, korder = 1,
    y = y, x = x, w = w,
    theta = theta1, z = z1, u = u1,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  expect_equal(result_k1$num_na_theta, 0)
  expect_length(result_k1$theta, n)
  expect_length(result_k1$z, n - 1)
  
  # Test korder = 2
  theta2 <- log(abs(y) + 1)
  z2 <- rep(0, n - 2)
  u2 <- rep(0, n - 2)
  result_k2 <- admm_gauss_testing(
    M = 100, korder = 2,
    y = y, x = x, w = w,
    theta = theta2, z = z2, u = u2,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  expect_equal(result_k2$num_na_theta, 0)
  expect_length(result_k2$theta, n)
  expect_length(result_k2$z, n - 2)
  
  # Test korder = 3
  theta3 <- log(abs(y) + 1)
  z3 <- rep(0, n - 3)
  u3 <- rep(0, n - 3)
  result_k3 <- admm_gauss_testing(
    M = 100, korder = 3,
    y = y, x = x, w = w,
    theta = theta3, z = z3, u = u3,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  expect_equal(result_k3$num_na_theta, 0)
  expect_length(result_k3$theta, n)
  expect_length(result_k3$z, n - 3)
})

test_that("admm_gauss_testing with different rho values", {
  set.seed(321)
  n <- 50
  x <- 1:n
  korder <- 2
  y <- 25 * sin(2 * pi * x / 30) + rnorm(n, 0, 3)
  w <- rep(1, n)
  theta <- log(abs(y) + 1)
  z <- rep(0, n - korder)
  u <- rep(0, n - korder)
  lam_z <- 0.5
  
  # Small rho
  result_small <- admm_gauss_testing(
    M = 100, korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = 0.1, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  expect_equal(result_small$num_na_theta, 0)
  
  # Large rho
  result_large <- admm_gauss_testing(
    M = 100, korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = 10.0, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  expect_equal(result_large$num_na_theta, 0)
  
  # Both should converge
  expect_true(result_small$converged)
  expect_true(result_large$converged)
})

test_that("admm_gauss_testing with varying weights", {
  set.seed(654)
  n <- 50
  x <- 1:n
  korder <- 2
  y <- 30 + 15 * cos(2 * pi * x / 20) + rnorm(n, 0, 2)
  theta <- log(abs(y) + 1)
  z <- rep(0, n - korder)
  u <- rep(0, n - korder)
  rho <- 1.0
  lam_z <- 0.5
  
  # Uniform weights
  w_uniform <- rep(1, n)
  result_uniform <- admm_gauss_testing(
    M = 100, korder = korder,
    y = y, x = x, w = w_uniform,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  expect_equal(result_uniform$num_na_theta, 0)
  
  # Varying weights
  w_varying <- exp(-(x - n/2)^2 / (2 * (n/4)^2))
  result_varying <- admm_gauss_testing(
    M = 100, korder = korder,
    y = y, x = x, w = w_varying,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  expect_equal(result_varying$num_na_theta, 0)
  
  # Results should differ due to different weights
  expect_false(isTRUE(all.equal(result_uniform$theta, result_varying$theta)))
})

test_that("admm_gauss_testing with unequally spaced x", {
  set.seed(987)
  n <- 40
  # Unequally spaced observations
  x <- cumsum(runif(n, 0.5, 1.5))
  korder <- 2
  y <- 20 + 10 * sin(2 * pi * x / max(x)) + rnorm(n, 0, 2)
  w <- rep(1, n)
  theta <- log(abs(y) + 1)
  z <- rep(0, n - korder)
  u <- rep(0, n - korder)
  rho <- 1.0
  lam_z <- 0.5
  
  result <- admm_gauss_testing(
    M = 100, korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  
  expect_equal(result$num_na_theta, 0)
  expect_equal(result$num_inf_theta, 0)
  expect_length(result$theta, n)
  expect_true(result$converged)
})

test_that("admm_gauss_testing convergence behavior", {
  set.seed(111)
  n <- 50
  x <- 1:n
  korder <- 2
  y <- 40 * sin(2 * pi * x / 25) + rnorm(n, 0, 3)
  w <- rep(1, n)
  theta <- log(abs(y) + 1)
  z <- rep(0, n - korder)
  u <- rep(0, n - korder)
  rho <- 1.0
  lam_z <- 0.5
  
  # Tight tolerance
  result_tight <- admm_gauss_testing(
    M = 200, korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-8
  )
  
  # Loose tolerance
  result_loose <- admm_gauss_testing(
    M = 200, korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-4
  )
  
  # Both should be valid
  expect_equal(result_tight$num_na_theta, 0)
  expect_equal(result_loose$num_na_theta, 0)
  
  # Tight tolerance might require more iterations
  expect_true(result_tight$niter >= result_loose$niter || 
              (result_tight$niter < 200 && result_loose$niter < 200))
  
  # Results should be similar
  expect_gt(cor(result_tight$theta, result_loose$theta), 0.99)
})

test_that("admm_gauss_testing diagnostic output for KF solver NA issue", {
  # Simple test case to diagnose KF solver
  set.seed(42)
  n <- 30
  x <- 1:n
  korder <- 2
  y <- rep(10, n) + rnorm(n, 0, 1)  # Nearly constant signal
  w <- rep(1, n)
  theta <- log(abs(y) + 1)
  z <- rep(0, n - korder)
  u <- rep(0, n - korder)
  rho <- 1.0
  lam_z <- 0.5
  
  cat("\n=== Diagnostic: Simple constant signal ===\n")
  cat("Input: n =", n, ", korder =", korder, "\n")
  cat("y range:", range(y), "\n")
  cat("theta_init range:", range(theta), "\n\n")
  
  # Test QR solver
  cat("Testing QR solver...\n")
  result_qr <- admm_gauss_testing(
    M = 100, korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 1, tol = 1e-6
  )
  cat("QR: niter =", result_qr$niter, 
      ", converged =", result_qr$converged, "\n")
  cat("QR: NA count =", result_qr$num_na_theta,
      ", Inf count =", result_qr$num_inf_theta, "\n")
  cat("QR: theta range:", range(result_qr$theta), "\n\n")
  
  # Test KF solver
  cat("Testing KF solver...\n")
  result_kf <- admm_gauss_testing(
    M = 100, korder = korder,
    y = y, x = x, w = w,
    theta = theta, z = z, u = u,
    rho = rho, lam_z = lam_z,
    linear_solver = 2, tol = 1e-6
  )
  cat("KF: niter =", result_kf$niter,
      ", converged =", result_kf$converged, "\n")
  cat("KF: NA count =", result_kf$num_na_theta,
      ", Inf count =", result_kf$num_inf_theta, "\n")
  
  if (result_kf$num_na_theta == 0) {
    cat("KF: theta range:", range(result_kf$theta), "\n")
    cat("Correlation QR-KF:", cor(result_qr$theta, result_kf$theta), "\n")
  } else {
    cat("KF: PRODUCED NA VALUES - cannot compute range\n")
    cat("KF: First 10 theta values:", result_kf$theta[1:10], "\n")
  }
  
  # This test always passes - it's just for diagnostic output
  expect_true(TRUE)
})
