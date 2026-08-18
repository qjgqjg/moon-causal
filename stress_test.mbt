///|
test "stress moderate sample" {
  let config = SyntheticConfig::new(250, 5, 1.25, 777UL)
  let data = simulate_observational_data(config)
  let model = fit_logistic_regression(
    data.dataset.covariates,
    data.dataset.treatment,
    max_iterations=250,
  )
  let p = predict_propensity(model, data.dataset.covariates)
  let estimate = estimate_ipw_ate(
    data.dataset.outcome,
    data.dataset.treatment,
    p,
  )
  assert_eq(p.length(), 250)
  assert_true(is_finite(estimate.estimate))
  assert_true(estimate.effective_sample_size > 0.0)
  assert_true(auc_score(p, data.dataset.treatment) >= 0.0)
}

///|
test "stress repeated analysis" {
  for seed in [1UL, 2UL, 3UL, 4UL, 5UL] {
    let data = simulate_observational_data(
      SyntheticConfig::new(30, 2, 0.75, seed),
    )
    let analysis = run_observational_analysis(data.dataset)
    assert_true(is_finite(analysis.ipw_estimate.estimate))
    assert_true(analysis.balance.length() == 2)
    assert_true(analysis.overlap.overlap_fraction >= 0.0)
  }
}

///|
test "stress estimand curve" {
  let data = simulate_observational_data(
    SyntheticConfig::new(100, 2, 1.0, 100UL),
  )
  let model = fit_logistic_regression(
    data.dataset.covariates,
    data.dataset.treatment,
  )
  let p = predict_propensity(model, data.dataset.covariates)
  let curve = effect_curve_by_trim(
    data.dataset.outcome,
    data.dataset.treatment,
    p,
    [0.0, 0.01, 0.05, 0.1, 0.2],
  )
  assert_eq(curve.length(), 5)
  for estimate in curve {
    assert_true(is_finite(estimate.estimate))
  }
}

///|
test "stress bootstrap moderate" {
  let data = simulate_observational_data(
    SyntheticConfig::new(50, 2, 1.0, 200UL),
  )
  let model = fit_logistic_regression(
    data.dataset.covariates,
    data.dataset.treatment,
  )
  let p = predict_propensity(model, data.dataset.covariates)
  let result = bootstrap_ipw_ate(
    data.dataset.outcome,
    data.dataset.treatment,
    p,
    100,
    300UL,
  )
  assert_eq(result.replicates, 100)
  assert_true(result.successful_replicates > 0)
}

///|
test "stress matrix dimension" {
  for size in 1..<=8 {
    let matrix = identity_matrix(size)
    assert_eq(matrix.length(), size)
    assert_eq(matrix[0].length(), size)
    let inverse = matrix_inverse(matrix)
    assert_eq(inverse.length(), size)
  }
}

///|
test "stress folds sizes" {
  for size in [1, 2, 3, 10, 31, 64] {
    let folds = make_folds(size, 7, 42UL)
    let total = folds.fold(init=0, (sum, fold) => {
      sum + fold.validation_indices.length()
    })
    assert_eq(total, size)
  }
}

///|
test "stress survival sizes" {
  let times = Array::new(capacity=40)
  let events = Array::new(capacity=40)
  for i in 0..<40 {
    times.push((i % 10).to_double())
    events.push(i % 3 != 0)
  }
  let curve = kaplan_meier(times, events)
  assert_true(curve.length() <= 10)
  assert_true(restricted_mean_survival(curve, 10.0) >= 0.0)
}
