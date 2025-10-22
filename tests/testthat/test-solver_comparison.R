test_that("estimate_rt works with Sparse QR solver", {
  # Create synthetic data
  set.seed(123)
  n <- 100
  x <- 1:n
  
  # Generate smooth Rt curve
  true_rt <- 1.15 + 0.35 * sin(2 * pi * x / 30)
  
  # Generate infection counts
  infections <- numeric(n)
  infections[1] <- 50
  for (i in 2:n) {
    expected <- infections[i-1] * true_rt[i] * 0.98
    infections[i] <- rpois(1, max(1, expected))
    infections[i] <- min(infections[i], 500)
  }
  
  # Simple delay distribution
  delay_distn <- discretize_gamma(1:15, shape = 5, scale = 1)
  
  # Calculate observed counts
  observed <- stats::filter(infections, rev(delay_distn), sides = 1)
  observed[is.na(observed)] <- 0
  observed <- round(pmax(observed, 1))
  
  # Test with Sparse QR solver
  result_qr <- estimate_rt(
    observed_counts = observed,
    x = x,
    delay_distn = delay_distn,
    nsol = 20,
    korder = 2,
    linear_solver = "qr"
  )
  
  # Check that result is valid
  expect_s3_class(result_qr, "poisson_rt")
  
  # Get fitted values - estimate_rt returns path of solutions
  fitted_vals <- fitted(result_qr)
  expect_false(any(is.na(fitted_vals)),
               label = "QR solver should not produce NAs in fitted values")
  
  # Check dimensions
  expect_equal(nrow(fitted_vals), n)
  
  # Predict to get Rt estimates
  rt_vals <- predict(result_qr)
  expect_false(any(is.na(rt_vals)),
               label = "QR solver should not produce NAs in Rt predictions")
  
  # Check that Rt values are reasonable (positive and finite)
  expect_true(all(rt_vals > 0 & is.finite(rt_vals)))
})

test_that("estimate_rt works with Kalman Filter solver", {
  # Create synthetic data
  set.seed(123)
  n <- 100
  x <- 1:n
  
  # Generate smooth Rt curve
  true_rt <- 1.15 + 0.35 * sin(2 * pi * x / 30)
  
  # Generate infection counts
  infections <- numeric(n)
  infections[1] <- 50
  for (i in 2:n) {
    expected <- infections[i-1] * true_rt[i] * 0.98
    infections[i] <- rpois(1, max(1, expected))
    infections[i] <- min(infections[i], 500)
  }
  
  # Simple delay distribution
  delay_distn <- discretize_gamma(1:15, shape = 5, scale = 1)
  
  # Calculate observed counts
  observed <- stats::filter(infections, rev(delay_distn), sides = 1)
  observed[is.na(observed)] <- 0
  observed <- round(pmax(observed, 1))
  
  # Test with Kalman Filter solver
  result_kf <- estimate_rt(
    observed_counts = observed,
    x = x,
    delay_distn = delay_distn,
    nsol = 20,
    korder = 1,
    linear_solver = "kf"
  )
  
  # Check that result is valid
  expect_s3_class(result_kf, "poisson_rt")
  
  # Get fitted values
  fitted_vals <- fitted(result_kf)
  expect_false(any(is.na(fitted_vals)),
               label = "KF solver should not produce NAs in fitted values")
  
  # Check dimensions
  expect_equal(nrow(fitted_vals), n)
  
  # Predict to get Rt estimates
  rt_vals <- predict(result_kf)
  expect_false(any(is.na(rt_vals)),
               label = "KF solver should not produce NAs in Rt predictions")
  
  # Check that Rt values are reasonable (positive and finite)
  expect_true(all(rt_vals > 0 & is.finite(rt_vals)))
})

test_that("QR and KF solvers produce similar results in estimate_rt", {
  # Create synthetic data
  set.seed(123)
  n <- 100
  x <- 1:n
  
  # Generate smooth Rt curve
  true_rt <- 1.15 + 0.35 * sin(2 * pi * x / 30)
  
  # Generate infection counts
  infections <- numeric(n)
  infections[1] <- 50
  for (i in 2:n) {
    expected <- infections[i-1] * true_rt[i] * 0.98
    infections[i] <- rpois(1, max(1, expected))
    infections[i] <- min(infections[i], 500)
  }
  
  # Simple delay distribution
  delay_distn <- discretize_gamma(1:15, shape = 5, scale = 1)
  
  # Calculate observed counts
  observed <- stats::filter(infections, rev(delay_distn), sides = 1)
  observed[is.na(observed)] <- 0
  observed <- round(pmax(observed, 1))
  
  # Test with both solvers
  result_qr <- estimate_rt(
    observed_counts = observed,
    x = x,
    delay_distn = delay_distn,
    nsol = 20,
    korder = 1,
    linear_solver = "qr"
  )
  
  result_kf <- estimate_rt(
    observed_counts = observed,
    x = x,
    delay_distn = delay_distn,
    nsol = 20,
    korder = 1,
    linear_solver = "kf"
  )
  
  # Both should produce valid results with no NAs
  fitted_qr <- fitted(result_qr)[, 10]
  fitted_kf <- fitted(result_kf)[, 10]
  
  expect_false(any(is.na(fitted_qr)))
  expect_false(any(is.na(fitted_kf)))
  
  # Fitted values should be highly correlated
  cor_val <- cor(fitted_qr, fitted_kf)
  expect_gt(cor_val, 0.9,
            label = "Fitted values from QR and KF should be highly correlated (r > 0.9)")
})

test_that("Sparse QR solver handles small datasets", {
  set.seed(456)
  n <- 30
  x <- 1:n
  y <- rpois(n, 10)
  delay_distn <- discretize_gamma(1:5, shape = 2, scale = 1)
  
  expect_no_error(
    result <- estimate_rt(
      observed_counts = y,
      x = x,
      delay_distn = delay_distn,
      nsol = 10,
      korder = 2,
      linear_solver = "qr"
    )
  )
  
  expect_equal(sum(is.na(result$rt)), 0)
})
