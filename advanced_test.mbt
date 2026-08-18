///|
fn assert_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit raise {
  if (actual - expected).abs() > tolerance {
    fail("expected {expected}, got {actual}")
  }
}

///|
test "linear algebra primitives and singular boundaries" {
  assert_close(dot([1.0, 2.0], [3.0, 4.0]), 11.0, 1.0e-12)
  assert_eq(dot([1.0], []), 0.0)
  assert_close(l2_norm([3.0, 4.0]), 5.0, 1.0e-12)
  assert_close(l1_norm([-2.0, 3.0]), 5.0, 1.0e-12)
  assert_eq(matrix_multiply([], []), [])
  assert_eq(
    solve_linear_system([[2.0, 1.0], [1.0, 3.0]], [5.0, 6.0]).length(),
    2,
  )
  assert_close(
    solve_linear_system([[2.0, 1.0], [1.0, 3.0]], [5.0, 6.0])[0],
    1.8,
    1.0e-10,
  )
  assert_eq(solve_linear_system([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0]), [])
  assert_close(matrix_inverse([[2.0, 0.0], [0.0, 4.0]])[1][1], 0.25, 1.0e-12)
  assert_eq(matrix_inverse([]), [])
}

///|
test "stable numerical guards" {
  assert_close(sigmoid(1000.0), 1.0, 1.0e-12)
  assert_close(sigmoid(-1000.0), 0.0, 1.0e-12)
  assert_close(safe_probability(-2.0), 1.0e-8, 1.0e-14)
  assert_close(safe_probability(2.0), 1.0 - 1.0e-8, 1.0e-14)
  assert_eq(quantile([], 0.5), 0.0)
  assert_close(quantile([4.0, 1.0, 3.0, 2.0], 0.5), 2.5, 1.0e-12)
  assert_close(effective_sample_size([1.0, 1.0, 1.0]), 3.0, 1.0e-12)
  assert_close(effective_sample_size([1.0, 0.0, 0.0]), 1.0, 1.0e-12)
  assert_eq(standardize([2.0, 2.0]), [0.0, 0.0])
}

///|
test "dataset construction and selection" {
  let dataset = CausalDataset::new([[1.0, 2.0], [3.0, 4.0]], [true, false], [
    4.0, 5.0,
  ])
  assert_true(dataset.is_valid())
  assert_eq(dataset.n(), 2)
  assert_eq(dataset.p(), 2)
  assert_eq(dataset.treated_count(), 1)
  assert_eq(dataset.control_count(), 1)
  assert_eq(dataset.feature_names, ["x1", "x2"])
  assert_eq(dataset.select([1]).outcome, [5.0])
  assert_true(
    CausalDataset::new([[1.0], [2.0]], [true], [1.0]).is_valid() == false,
  )
}

///|
test "logistic model learns deterministic treatment" {
  let x = [[-3.0], [-2.0], [-1.0], [1.0], [2.0], [3.0]]
  let t = [false, false, false, true, true, true]
  let model = fit_logistic_regression(
    x,
    t,
    learning_rate=0.4,
    max_iterations=1000,
  )
  let scores = predict_propensity(model, x)
  assert_true(scores[0] < scores[5])
  assert_true(classification_accuracy(scores, t) >= 0.8)
  assert_true(model.loss >= 0.0)
  assert_eq(predict_propensity(model, []).length(), 0)
}

///|
test "linear outcome model and metrics" {
  let x = [[0.0], [1.0], [2.0], [3.0]]
  let y = [1.0, 3.0, 5.0, 7.0]
  let model = fit_linear_outcome_model(x, y)
  let predictions = predict_outcomes(model, x)
  assert_close(predictions[0], 1.0, 1.0e-6)
  assert_close(root_mean_squared_error_raw(y, predictions), 0.0, 1.0e-6)
  assert_close(r_squared(y, predictions), 1.0, 1.0e-6)
  assert_close(correlation(y, predictions), 1.0, 1.0e-6)
  assert_eq(residuals(model, x, y).length(), 4)
}

///|
test "estimators handle overlap and extreme scores" {
  let y = [3.0, 4.0, 1.0, 2.0]
  let t = [true, true, false, false]
  let p = [0.99, 0.8, 0.2, 0.01]
  let weights = stabilized_ipw_weights(t, p)
  assert_eq(weights.length(), 4)
  assert_true(weight_design_effect(weights) >= 1.0)
  let estimate = estimate_ipw_ate(y, t, p)
  assert_true(is_finite(estimate.estimate))
  let att = estimate_ipw_att(y, t, p)
  assert_true(is_finite(att.standard_error))
  assert_eq(trim_propensity_scores([-1.0, 0.5, 2.0]), [0.01, 0.5, 0.99])
  assert_eq(ipw_weights_ate([true], [0.0]).length(), 1)
  assert_true(is_finite(ipw_weights_ate([true], [0.0])[0]))
}

///|
test "matching and balance diagnostics" {
  let x = [[0.0], [0.1], [0.2], [2.0], [2.1], [2.2]]
  let t = [true, true, true, false, false, false]
  let pairs = nearest_neighbor_match(x, t, caliper=3.0)
  assert_eq(pairs.length(), 3)
  assert_close(matching_mean_distance(pairs), 1.3433209885, 1.0e-9)
  let y = [5.0, 5.1, 4.9, 3.0, 3.1, 2.9]
  assert_true(estimate_matched_ate(y, pairs).estimate > 1.0)
  let metrics = balance_table(x, t, ["x1"])
  assert_eq(metrics.length(), 1)
  assert_true(maximum_absolute_smd(metrics) > 0.0)
  assert_eq(balanced_covariate_count(metrics), 0)
  assert_close(
    weighted_smd([1.0, 2.0], [true, false], [1.0, 1.0]),
    0.0,
    1.0e-12,
  )
}

///|
test "graph traversal and cycle detection" {
  let graph = CausalGraph::new(["age", "treatment", "outcome"])
  let graph = graph
    .add_edge("age", "treatment")
    .add_edge("age", "outcome")
    .add_edge("treatment", "outcome")
  assert_true(graph.has_edge("age", "treatment"))
  assert_eq(graph.edge_count(), 3)
  assert_eq(graph.parents("outcome"), ["age", "treatment"])
  assert_eq(graph.ancestors("outcome").length(), 2)
  assert_eq(graph.topological_order().length(), 3)
  let cyclic = graph.add_edge("outcome", "age")
  assert_eq(cyclic.topological_order(), [])
  assert_eq(graph.backdoor_candidates("treatment", "outcome"), ["age"])
}

///|
test "preprocessing and feature engineering boundaries" {
  let x = [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]]
  let summary = summarize_covariates(x)
  assert_eq(summary.rows, 3)
  assert_eq(summary.columns, 2)
  assert_close(summary.means[0], 2.0, 1.0e-12)
  assert_close(scale_covariates(x)[0][0], -1.224744871, 1.0e-6)
  assert_close(center_covariates(x)[1][1], 0.0, 1.0e-12)
  assert_eq(column_values(x, 8), [])
  assert_close(median([1.0, 9.0, 2.0]), 2.0, 1.0e-12)
  assert_close(median_absolute_deviation([1.0, 2.0, 3.0]), 1.0, 1.0e-12)
  assert_true(winsorize([0.0, 1.0, 100.0], 0.2, 0.8)[2] < 100.0)
  assert_eq(polynomial_features(x, 2)[0].length(), 4)
  assert_eq(interaction_features(x)[0].length(), 3)
  assert_eq(rank_transform([1.0, 2.0, 3.0])[0], 0.0)
  assert_eq(
    train_validation_split(
      CausalDataset::new(x, [true, false, true], [1.0, 2.0, 3.0]),
      0.5,
      3UL,
    ).validation.n(),
    1,
  )
}

///|
test "diagnostics, metrics, and binary outcomes" {
  let p = [0.1, 0.2, 0.8, 0.9]
  let t = [false, false, true, true]
  let overlap = overlap_report(p, t)
  assert_close(overlap.overlap_fraction, 0.0, 1.0e-12)
  assert_eq(positivity_violations([0.0, 0.5, 1.0], 0.01), 2)
  assert_close(auc_score(p, t), 1.0, 1.0e-12)
  assert_close(brier_score(p, t), 0.025, 1.0e-12)
  assert_true(log_loss(p, t) > 0.0)
  assert_close(mean_absolute_error_raw([1.0, 2.0], [2.0, 4.0]), 1.5, 1.0e-12)
  let risk = risk_measure([true, false, true, false], [true, true, false, false])
  assert_close(risk.treated_risk, 0.5, 1.0e-12)
  assert_close(risk.control_risk, 0.5, 1.0e-12)
  assert_close(risk.risk_difference, 0.0, 1.0e-12)
  assert_eq(binary_outcome_to_double([true, false]), [1.0, 0.0])
  assert_eq(threshold_outcome([0.1, 0.9], 0.5), [false, true])
}

///|
test "uncertainty and sensitivity" {
  let y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let t = [false, false, false, true, true, true]
  let bootstrap = bootstrap_difference_in_means(y, t, replicates=100, seed=11UL)
  assert_eq(bootstrap.replicates, 100)
  assert_true(bootstrap.successful_replicates > 0)
  assert_true(bootstrap.upper >= bootstrap.lower)
  let p = rosenbaum_probability_bounds(0.05, 2.0)
  assert_true(p.upper_probability > p.lower_probability)
  assert_close(e_value(1.0), 1.0, 1.0e-12)
  assert_true(critical_gamma(0.01, 0.05, 10.0) > 1.0)
  let mediation = mediation_decomposition([5.0], [3.0], [2.0], [1.0])
  assert_close(mediation.total_effect, 4.0, 1.0e-12)
  assert_true(mediation_balance_check(mediation, 1.0e-12))
}

///|
test "panel and dose response" {
  let y = [1.0, 2.0, 1.0, 1.0, 2.0, 4.0, 1.0, 2.0]
  let g = [true, true, false, false, true, true, false, false]
  let post = [false, false, false, false, true, true, true, true]
  let did = difference_in_differences(y, g, post)
  assert_close(did.estimate, 1.0, 1.0e-12)
  let events = event_study(y, g, [-1, 0, 1, 2, -1, 0, 1, 2])
  assert_eq(events.length(), 4)
  let curve = dose_response_curve([0.0, 0.5, 1.0, 1.5], [1.0, 2.0, 2.5, 4.0], 2)
  assert_eq(curve.length(), 2)
  assert_true(curve[1].mean_outcome > curve[0].mean_outcome)
  let weights = synthetic_control_weights([1.0, 2.0], [[1.0, 2.0], [4.0, 5.0]])
  assert_close(sum(weights), 1.0, 1.0e-12)
}

///|
test "synthetic benchmark and full pipeline" {
  let config = SyntheticConfig::new(80, 3, 2.0, 42UL)
  let simulated = simulate_observational_data(config)
  assert_eq(simulated.dataset.n(), 80)
  assert_true(simulated.dataset.is_valid())
  let report = run_synthetic_benchmark(config)
  assert_eq(report.sample_size, 80)
  assert_true(report.effective_sample_size > 0.0)
  let analysis = run_observational_analysis(simulated.dataset)
  assert_eq(analysis_summary_vector(analysis).length(), 5)
  assert_true(analysis.quality_score() >= 0.0)
  assert_true(is_finite(doubly_robust_from_dataset(analysis).estimate))
  assert_true(
    cross_validate_propensity(
      simulated.dataset.covariates,
      simulated.dataset.treatment,
      4,
      5UL,
    ) >=
    0.0,
  )
  assert_true(
    cross_validate_outcome_rmse(
      simulated.dataset.covariates,
      simulated.dataset.outcome,
      4,
      5UL,
    ) >=
    0.0,
  )
}
