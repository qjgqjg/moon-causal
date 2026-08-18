///|
test "simulation sizes and treatment rates" {
  for n in [0, 1, 5, 20, 100] {
    let config = SyntheticConfig::new(n, 2, 1.5, (n + 1).to_uint64())
    let simulated = simulate_observational_data(config)
    assert_eq(simulated.dataset.n(), n)
    assert_eq(simulated.true_propensity.length(), n)
    assert_true(simulated.dataset.is_valid())
  }
}

///|
test "simulation covariate dimensions" {
  for columns in [1, 2, 4, 8] {
    let simulated = simulate_observational_data(
      SyntheticConfig::new(30, columns, 2.0, 10UL),
    )
    assert_eq(simulated.dataset.p(), columns)
    for row in simulated.dataset.covariates {
      assert_eq(row.length(), columns)
    }
  }
}

///|
test "simulation benchmark is finite" {
  for effect in [-2.0, 0.0, 1.0, 3.0] {
    let report = run_synthetic_benchmark(
      SyntheticConfig::new(120, 3, effect, 88UL),
    )
    assert_true(is_finite(report.estimated_ate))
    assert_true(is_finite(report.absolute_error))
    assert_true(report.effective_sample_size > 0.0)
    assert_true(report.propensity_auc >= 0.0 && report.propensity_auc <= 1.0)
  }
}

///|
test "pipeline exposes diagnostics" {
  let simulated = simulate_observational_data(
    SyntheticConfig::new(50, 2, 1.0, 9UL),
  )
  let analysis = run_observational_analysis(simulated.dataset)
  assert_eq(analysis.propensity_scores.length(), 50)
  assert_eq(analysis.balance.length(), 2)
  assert_true(analysis.overlap.minimum >= 0.0)
  assert_true(analysis.overlap.maximum <= 1.0)
  let vector = analysis_summary_vector(analysis)
  for value in vector {
    assert_true(is_finite(value))
  }
}

///|
test "train validation split preserves dimensions" {
  let x = [[1.0], [2.0], [3.0], [4.0], [5.0], [6.0]]
  let data = CausalDataset::new(x, [true, false, true, false, true, false], [
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
  ])
  for fraction in [0.0, 0.2, 0.5, 0.9, 1.0] {
    let split = train_validation_split(data, fraction, 3UL)
    assert_eq(split.training.n() + split.validation.n(), data.n())
    assert_eq(split.training.p(), data.p())
    if split.validation.n() > 0 {
      assert_eq(split.validation.p(), data.p())
    }
  }
}

///|
test "model cross validation metrics are finite" {
  let x = [[0.0], [1.0], [2.0], [3.0], [4.0], [5.0], [6.0], [7.0]]
  let t = [false, false, false, true, true, true, true, true]
  let y = [0.0, 1.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  assert_true(is_finite(cross_validate_propensity(x, t, 4, 1UL)))
  assert_true(is_finite(cross_validate_outcome_rmse(x, y, 4, 1UL)))
}
