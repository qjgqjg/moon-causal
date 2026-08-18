///|
test "quality estimator agreement" {
  let x = [
    [-2.0, 0.0],
    [-1.0, 1.0],
    [0.0, 0.0],
    [1.0, 1.0],
    [2.0, 0.0],
    [3.0, 1.0],
  ]
  let t = [false, false, false, true, true, true]
  let y = [0.0, 1.0, 1.5, 3.0, 4.0, 5.0]
  let model = fit_logistic_regression(x, t, max_iterations=500)
  let p = predict_propensity(model, x)
  let ipw = estimate_ipw_ate(y, t, p)
  let overlap = estimate_overlap_ate(y, t, p)
  let atc = estimate_atc(y, t, p)
  assert_true(is_finite(ipw.estimate))
  assert_true(is_finite(overlap.estimate))
  assert_true(is_finite(atc.estimate))
  assert_true(ipw.effective_sample_size > 0.0)
  assert_true(overlap.effective_sample_size > 0.0)
}

///|
test "quality balance after trimming" {
  let t = [false, false, false, true, true, true, true, false]
  let p = [0.02, 0.15, 0.25, 0.45, 0.55, 0.75, 0.9, 0.98]
  let raw = smd([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], t)
  let weight = overlap_weights(t, p)
  let adjusted = weighted_smd(
    [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
    t,
    weight,
  )
  assert_true(is_finite(raw) && is_finite(adjusted))
  assert_true(positivity_violations(p, 0.05) == 2)
}

///|
test "quality reports use stable fields" {
  let data = CausalDataset::new(
    [[1.0], [2.0], [3.0], [4.0]],
    [true, false, true, false],
    [2.0, 1.0, 4.0, 3.0],
  )
  let analysis = run_observational_analysis(data)
  let report = overlap_report(analysis.propensity_scores, data.treatment)
  assert_true(report.minimum >= 0.0)
  assert_true(report.maximum <= 1.0)
  assert_eq(
    smd_profile(data.covariates, data.treatment, data.feature_names, 0.1).before.length(),
    1,
  )
  assert_eq(analysis_summary_vector(analysis).length(), 5)
}

///|
test "quality confidence intervals contain point" {
  let estimates = [
    Estimate::from_standard_error(0.0, 0.0, 1, 1.0, "zero"),
    Estimate::from_standard_error(1.0, 0.1, 10, 9.0, "positive"),
    Estimate::from_standard_error(-1.0, 0.1, 10, 9.0, "negative"),
  ]
  for estimate in estimates {
    assert_true(estimate.lower <= estimate.estimate)
    assert_true(estimate.upper >= estimate.estimate)
    assert_true(estimate.as_vector().length() == 5)
  }
}

///|
test "quality continuous effects" {
  let dose = [0.0, 1.0, 2.0, 3.0, 4.0]
  let outcome = [0.0, 1.0, 4.0, 9.0, 16.0]
  let x = [[0.0], [1.0], [2.0], [3.0], [4.0]]
  let curve = dose_response_curve(dose, outcome, 5)
  assert_eq(curve.length(), 5)
  let gps = generalized_propensity_scores(dose, x)
  assert_eq(gps.length(), 5)
  for value in gps {
    assert_true(value >= 0.0)
  }
}

///|
test "quality event study" {
  let y = [1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 2.0, 2.0]
  let g = [true, true, true, true, false, false, false, false]
  let time = [-2, -1, 0, 1, -2, -1, 0, 1]
  let result = event_study(y, g, time)
  assert_eq(result.length(), 4)
  for point in result {
    assert_true(point.count > 0)
  }
}

///|
test "quality imputation report" {
  let missing = [[1.0], [(-1.0).sqrt()], [3.0], [(-1.0).sqrt()]]
  let (filled, summary) = impute_non_finite(missing)
  assert_eq(summary.rows, 4)
  assert_eq(summary.columns, 1)
  assert_eq(summary.imputed_values, 2)
  assert_close(filled[1][0], 2.0, 1.0e-12)
  assert_close(filled[3][0], 2.0, 1.0e-12)
}
