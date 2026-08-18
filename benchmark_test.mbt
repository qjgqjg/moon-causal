///|
test "benchmark deterministic sample sizes" {
  for size in [10, 25, 50, 75] {
    let report = run_synthetic_benchmark(
      SyntheticConfig::new(size, 2, 1.0, 1234UL),
    )
    assert_eq(report.sample_size, size)
    assert_true(report.effective_sample_size > 0.0)
    assert_true(report.standard_error >= 0.0)
    assert_true(report.absolute_error >= 0.0)
  }
}

///|
test "benchmark treatment effects" {
  let zero = run_synthetic_benchmark(SyntheticConfig::new(80, 2, 0.0, 5UL))
  let positive = run_synthetic_benchmark(SyntheticConfig::new(80, 2, 2.0, 5UL))
  assert_true(is_finite(zero.estimated_ate))
  assert_true(is_finite(positive.estimated_ate))
  assert_true(positive.true_ate > zero.true_ate)
}

///|
test "benchmark dimensions" {
  for columns in [1, 2, 3, 5] {
    let data = simulate_observational_data(
      SyntheticConfig::new(20, columns, 1.0, 55UL),
    )
    assert_eq(data.dataset.p(), columns)
    assert_eq(data.dataset.covariates.length(), 20)
  }
}

///|
test "benchmark bootstrap" {
  let data = simulate_observational_data(SyntheticConfig::new(40, 2, 1.5, 33UL))
  let model = fit_logistic_regression(
    data.dataset.covariates,
    data.dataset.treatment,
  )
  let p = predict_propensity(model, data.dataset.covariates)
  let summary = bootstrap_ipw_ate(
    data.dataset.outcome,
    data.dataset.treatment,
    p,
    50,
    22UL,
  )
  assert_eq(summary.replicates, 50)
  assert_true(summary.successful_replicates > 0)
  assert_true(summary.upper >= summary.lower)
}

///|
test "benchmark cross validation" {
  let data = simulate_observational_data(SyntheticConfig::new(60, 3, 1.0, 12UL))
  let accuracy = cross_validate_propensity(
    data.dataset.covariates,
    data.dataset.treatment,
    5,
    4UL,
  )
  let rmse = cross_validate_outcome_rmse(
    data.dataset.covariates,
    data.dataset.outcome,
    5,
    4UL,
  )
  assert_true(accuracy >= 0.0 && accuracy <= 1.0)
  assert_true(rmse >= 0.0)
}
