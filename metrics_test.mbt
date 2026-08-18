///|
test "metric identities" {
  let y = [1.0, 2.0, 3.0, 4.0]
  let p = [1.0, 2.0, 3.0, 4.0]
  assert_close(covariance(y, p), variance(y), 1.0e-12)
  assert_close(correlation(y, p), 1.0, 1.0e-12)
  assert_close(r_squared(y, p), 1.0, 1.0e-12)
  assert_close(root_mean_squared_error_raw(y, p), 0.0, 1.0e-12)
  assert_close(mean_absolute_error_raw(y, p), 0.0, 1.0e-12)
  assert_close(
    standardized_effect_size(2.0, y),
    2.0 / std_dev(y, sample=false),
    1.0e-12,
  )
}

///|
test "probability metrics" {
  let probabilities = [0.2, 0.4, 0.6, 0.8]
  let treatment = [false, false, true, true]
  assert_true(log_loss(probabilities, treatment) > 0.0)
  assert_true(brier_score(probabilities, treatment) >= 0.0)
  assert_true(calibration_error(probabilities, treatment, 2) >= 0.0)
  assert_close(treatment_t_statistic(2.0, 1.0), 2.0, 1.0e-12)
  assert_eq(treatment_t_statistic(2.0, 0.0), 0.0)
}

///|
test "robust statistic ordering" {
  let values = [100.0, 1.0, 2.0, 3.0, 4.0]
  assert_true(trimmed_mean(values, 0.2) < mean(values))
  assert_true(winsorized_mean(values, 0.2) < mean(values))
  assert_true(median(values) < mean(values))
  assert_true(median_absolute_deviation(values) > 0.0)
  let z = robust_z_scores(values)
  assert_eq(z.length(), values.length())
  assert_true(outlier_indices(values, 3.0).length() >= 1)
}

///|
test "rank and empirical distribution" {
  let values = [3.0, 1.0, 2.0, 2.0]
  let ranks = rank_transform(values)
  assert_true(ranks[1] < ranks[0])
  assert_eq(empirical_cdf(values, [0.0, 2.0, 4.0]), [0.0, 0.75, 1.0])
  assert_close(
    weighted_quantile(values, [1.0, 1.0, 1.0, 1.0], 0.5),
    2.0,
    1.0e-12,
  )
}
