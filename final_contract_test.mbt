///|
test "public contract arrays" {
  let arrays = [
    estimate_propensity_scores([[0.0], [1.0]], [1.0]),
    ipw_weights_ate([true, false], [0.5, 0.5]),
    ipw_weights_att([true, false], [0.5, 0.5]),
    stabilized_ipw_weights([true, false], [0.5, 0.5]),
    overlap_weights([true, false], [0.5, 0.5]),
    overlap_population_weights([0.5, 0.5]),
  ]
  for array in arrays {
    assert_eq(array.length(), 2)
  }
}

///|
test "public contract scalar stability" {
  let scalars = [
    sigmoid(0.0),
    mean([1.0, 2.0]),
    variance([1.0, 2.0]),
    std_dev([1.0, 2.0]),
    rosenbaum_bounds_gamma(0.05, 2.0),
    safe_logit(0.5),
    e_value(1.5),
  ]
  for scalar in scalars {
    assert_true(is_finite(scalar))
  }
}

///|
test "public contract table outputs" {
  let covariates = [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0], [4.0, 8.0]]
  let treatment = [true, false, true, false]
  assert_eq(balance_table(covariates, treatment, ["a", "b"]).length(), 2)
  assert_eq(covariance_matrix(covariates).length(), 2)
  assert_eq(
    weighted_covariance_matrix(covariates, [1.0, 1.0, 1.0, 1.0]).length(),
    2,
  )
  assert_eq(summarize_covariates(covariates).columns, 2)
  assert_eq(make_design_matrix(covariates, true).column_names.length(), 3)
}

///|
test "public contract model outputs" {
  let x = [[0.0], [1.0], [2.0], [3.0]]
  let t = [false, false, true, true]
  let y = [1.0, 2.0, 4.0, 5.0]
  let propensity = fit_logistic_regression(x, t)
  let outcome = fit_linear_outcome_model(x, y)
  assert_eq(predict_propensity(propensity, x).length(), 4)
  assert_eq(predict_outcomes(outcome, x).length(), 4)
  assert_eq(model_quality_flags(propensity).length(), 4)
  assert_eq(marginal_effects(outcome, x).length(), 1)
}

///|
test "public contract graph outputs" {
  let graph = CausalGraph::new(["confounder", "treatment", "outcome"])
    .add_edge("confounder", "treatment")
    .add_edge("confounder", "outcome")
    .add_edge("treatment", "outcome")
  assert_true(graph.edge_count() == 3)
  assert_true(graph.topological_order().length() == 3)
  assert_true(graph.ancestors("outcome").contains("confounder"))
  assert_true(graph.children("confounder").length() == 2)
}

///|
test "public contract intervals" {
  let values = [1.0, 2.0, 3.0, 4.0]
  let interval = bootstrap_percentile_interval(values, 0.95)
  assert_true(interval.0 <= interval.1)
  let estimate = g_computation_ate(values, [0.0, 1.0, 2.0, 3.0])
  assert_true(estimate.lower <= estimate.upper)
  assert_eq(estimate.as_vector().length(), 5)
}

///|
test "public contract data transformations" {
  let data = CausalDataset::new(
    [[1.0], [2.0], [3.0], [4.0]],
    [true, false, true, false],
    [2.0, 1.0, 4.0, 3.0],
  )
  let split = train_validation_split(data, 0.5, 18UL)
  assert_eq(split.training.n() + split.validation.n(), 4)
  assert_eq(data.select([0, 2]).treatment_rate(), 1.0)
  assert_eq(data.select([1, 3]).treatment_rate(), 0.0)
  assert_eq(complete_case_dataset(data).n(), 4)
}

///|
test "public contract diagnostics" {
  let p = [0.1, 0.2, 0.7, 0.8]
  let t = [false, false, true, true]
  let report = overlap_report(p, t)
  assert_true(report.minimum <= report.maximum)
  assert_true(positivity_violations(p, 0.01) == 0)
  assert_eq(propensity_histogram(p, 4).length(), 4)
  assert_eq(propensity_calibration(p, t, 4).length(), 4)
}

///|
test "public contract multi arm" {
  let outcomes = [1.0, 2.0, 4.0, 5.0, 7.0, 8.0]
  let arms = [0, 0, 1, 1, 2, 2]
  assert_eq(arm_levels(arms), [0, 1, 2])
  assert_eq(arm_sample_counts(arms), [(0, 2), (1, 2), (2, 2)])
  assert_eq(arm_estimates(outcomes, arms, 0).length(), 3)
  assert_true(is_finite(pairwise_arm_difference(outcomes, arms, 2, 0).estimate))
}

///|
test "public contract reports" {
  let data = CausalDataset::new(
    [[0.0], [1.0], [2.0], [3.0]],
    [false, false, true, true],
    [1.0, 2.0, 4.0, 5.0],
  )
  let analysis = run_observational_analysis(data)
  let audit = audit_analysis(analysis)
  assert_true(audit.valid_dataset)
  assert_eq(analysis_summary_vector(analysis).length(), 5)
  assert_eq(audit_model(analysis.propensity_model).length(), 4)
  assert_eq(audit_estimate(analysis.ipw_estimate).length(), 5)
}

///|
test "public contract survival" {
  let curve = kaplan_meier([1.0, 2.0, 3.0, 4.0], [true, false, true, true])
  assert_true(curve.length() > 0)
  assert_true(survival_at(curve, 0.0) == 1.0)
  assert_true(survival_at(curve, 5.0) >= 0.0)
  assert_true(restricted_mean_survival(curve, 4.0) <= 4.0)
}
