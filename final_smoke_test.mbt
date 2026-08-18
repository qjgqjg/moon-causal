///|
test "final smoke numerical branches" {
  let values = [-10.0, -1.0, 0.0, 1.0, 10.0]
  for value in values {
    assert_true(is_finite(sigmoid(value)))
    assert_true(is_finite(safe_logit(sigmoid(value))))
  }
  let matrix = [[1.0, 2.0, 3.0], [0.0, 1.0, 4.0], [5.0, 6.0, 0.0]]
  let inverse = matrix_inverse(matrix)
  assert_eq(inverse.length(), 3)
  assert_close(matrix_multiply(matrix, inverse)[0][0], 1.0, 1.0e-8)
}

///|
test "final smoke estimation branches" {
  let y = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
  let t = [false, true, false, true, false, true]
  let p = [0.2, 0.4, 0.5, 0.6, 0.7, 0.8]
  let estimates = [
    estimate_difference_in_means(y, t),
    estimate_ipw_ate(y, t, p),
    estimate_ipw_att(y, t, p),
    estimate_atc(y, t, p),
    estimate_overlap_ate(y, t, p),
    estimate_from_strata(y, t, p, 3),
  ]
  for estimate in estimates {
    assert_true(is_finite(estimate.estimate))
    assert_true(estimate.lower <= estimate.upper)
  }
}

///|
test "final smoke preprocessing branches" {
  let x = [[0.0, 1.0], [1.0, 2.0], [2.0, 3.0], [3.0, 4.0]]
  let transformed = [
    center_covariates(x),
    scale_covariates(x),
    interaction_features(x),
    polynomial_features(x, 2),
  ]
  for matrix in transformed {
    assert_eq(matrix.length(), 4)
  }
  assert_eq(feature_means(x), [1.5, 2.5])
  assert_eq(feature_minima(x), [0.0, 1.0])
  assert_eq(feature_maxima(x), [3.0, 4.0])
}

///|
test "final smoke diagnostics branches" {
  let data = CausalDataset::new(
    [[0.0], [1.0], [2.0], [3.0]],
    [false, false, true, true],
    [1.0, 2.0, 4.0, 5.0],
  )
  let analysis = run_observational_analysis(data)
  let report = audit_analysis(analysis)
  assert_true(report.valid_dataset)
  assert_true(report.adequate_sample)
  assert_true(report.has_treatment_variation)
  assert_true(report.has_outcome_variation)
  assert_true(report.quality_score >= 0.0)
  assert_true(report.quality_score <= 1.0)
}

///|
test "final smoke resampling branches" {
  let data = [1.0, 2.0, 3.0, 4.0]
  let shuffled = shuffled_indices(4, 123UL)
  assert_eq(shuffled.length(), 4)
  assert_eq(systematic_sample(4, 2, 123UL).length(), 2)
  assert_eq(bootstrap_counts(4, 4, 123UL).length(), 4)
  assert_eq(delete_one_mean(data).length(), 4)
  assert_eq(influence_function(data), [-1.5, -0.5, 0.5, 1.5])
}

///|
test "final smoke structured outputs" {
  let graph = CausalGraph::new(["x", "t", "y"])
    .add_edge("x", "t")
    .add_edge("t", "y")
  assert_eq(graph.topological_order(), ["x", "t", "y"])
  let mediation = mediation_decomposition([4.0, 5.0], [3.0, 4.0], [2.0, 2.0], [
    1.0, 1.0,
  ])
  assert_true(mediation_balance_check(mediation, 1.0e-12))
  let risk = risk_measure([true, false, true, false], [true, true, false, false])
  assert_close(risk.risk_difference, 0.0, 1.0e-12)
}
