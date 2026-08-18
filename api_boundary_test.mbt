///|
test "dataset boundary transformations" {
  let empty = CausalDataset::new([], [], [])
  assert_true(empty.is_valid())
  assert_eq(empty.p(), 0)
  assert_eq(empty.treatment_rate(), 0.0)
  assert_eq(empty.outcome_range(), (0.0, 0.0))
  assert_eq(empty.column(0), [])
  assert_eq(empty.filter_by_score([], 0.0, 1.0).n(), 0)
  let malformed = CausalDataset::new([[1.0], [2.0, 3.0]], [true, false], [
    1.0, 2.0,
  ])
  assert_true(!malformed.is_valid())
}

///|
test "histogram and support boundaries" {
  assert_eq(propensity_histogram([], 3), [0, 0, 0])
  assert_eq(propensity_histogram([0.0, 1.0], 1), [2])
  assert_eq(propensity_histogram([0.0, 0.5, 1.0], 2).length(), 2)
  assert_eq(stratum_counts([], 3), [0, 0, 0])
  assert_eq(stratum_counts([0, 3, -1, 1], 2), [1, 1])
  assert_eq(common_support_indices([], 0.0, 1.0), [])
  assert_eq(common_support_indices([0.1, 0.9], 0.2, 0.8), [])
}

///|
test "odds boundary stability" {
  for probability in [0.0, 0.00000001, 0.5, 0.99999999, 1.0] {
    let logit = safe_logit(probability)
    let odds = odds_from_probability(probability)
    assert_true(is_finite(logit))
    assert_true(is_finite(odds))
    assert_true(probability_from_odds(odds) >= 0.0)
    assert_true(probability_from_odds(odds) <= 1.0)
  }
}

///|
test "model fit edge sizes" {
  let empty_model = fit_linear_outcome_model([], [])
  assert_true(empty_model.has_converged() == false)
  assert_eq(empty_model.coefficients, [])
  let single_model = fit_linear_outcome_model([[1.0]], [2.0])
  assert_eq(predict_outcomes(single_model, [[1.0]]).length(), 1)
  let logistic = fit_logistic_regression(
    [[1.0], [1.0]],
    [true, false],
    max_iterations=10,
  )
  assert_eq(predict_propensity(logistic, [[1.0], [1.0]]).length(), 2)
}

///|
test "estimate boundary intervals" {
  let empty = estimate_difference_in_means([], [])
  assert_eq(empty.sample_size, 0)
  assert_eq(empty.estimate, 0.0)
  let group_missing = estimate_ipw_ate([1.0], [true], [0.5])
  assert_true(is_finite(group_missing.estimate))
  let strata_missing = estimate_from_strata([1.0], [true], [0.5], 2)
  assert_true(is_finite(strata_missing.estimate))
}
