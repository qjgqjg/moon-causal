///|
test "convenience summaries" {
  let x = [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]]
  assert_eq(feature_means(x), [2.0, 20.0])
  assert_eq(feature_minima(x), [1.0, 10.0])
  assert_eq(feature_maxima(x), [3.0, 30.0])
  assert_eq(feature_standard_deviations(x).length(), 2)
  assert_eq(binned_counts([0.0, 0.5, 1.0], 2), [1, 2])
  assert_eq(binned_counts([], 2), [0, 0])
}

///|
test "convenience pairwise" {
  let values = [1.0, 3.0, 5.0]
  assert_eq(pairwise_difference(values), [
    [0.0, -2.0, -4.0],
    [2.0, 0.0, -2.0],
    [4.0, 2.0, 0.0],
  ])
  assert_eq(absolute_pairwise_difference(values), [
    [0.0, 2.0, 4.0],
    [2.0, 0.0, 2.0],
    [4.0, 2.0, 0.0],
  ])
  assert_eq(pairwise_difference([]), [])
}

///|
test "convenience treatment diagnostics" {
  let y = [1.0, 2.0, 4.0, 5.0]
  let t = [false, false, true, true]
  assert_true(treatment_outcome_covariance(y, t) > 0.0)
  assert_true(crude_effect_size(y, t) > 0.0)
  assert_true(is_finite(treatment_outcome_covariance([], [])))
  assert_eq(overlap_fraction([0.2, 0.3, 0.7, 0.8], t), 0.0)
  assert_close(common_support_width([0.2, 0.3, 0.7, 0.8], t), -0.4, 1.0e-12)
}

///|
test "convenience quality flags" {
  let estimate = Estimate::from_standard_error(1.0, 0.1, 100, 90.0, "ATE")
  let flags = estimate_quality_flags(estimate)
  assert_eq(flags.length(), 5)
  for flag in flags {
    assert_true(flag)
  }
  let model = fit_logistic_regression([[0.0], [1.0], [2.0], [3.0]], [
    false, false, true, true,
  ])
  let model_flags = model_quality_flags(model)
  assert_eq(model_flags.length(), 4)
  assert_true(model_flags[1])
  assert_true(model_flags[2])
}

///|
test "convenience probability groups" {
  let scores = [0.1, 0.2, 0.8, 0.9]
  let t = [false, false, true, true]
  assert_true(overlap_fraction(scores, t) == 0.0)
  assert_true(common_support_width(scores, t) < 0.0)
  let support = overlap_report([0.2, 0.4, 0.6, 0.8], t)
  assert_true(support.overlap_fraction >= 0.0)
  assert_true(support.overlap_fraction <= 1.0)
}

///|
test "convenience data summary repeated" {
  for size in [1, 2, 5, 10] {
    let row = Array::make(2, 1.0)
    let data = Array::make(size, row)
    assert_eq(feature_means(data), [1.0, 1.0])
    assert_eq(feature_minima(data), [1.0, 1.0])
    assert_eq(feature_maxima(data), [1.0, 1.0])
    assert_eq(binned_counts(feature_means(data), 3).length(), 3)
  }
}
