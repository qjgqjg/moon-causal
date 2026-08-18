///|
test "acceptance matrix validates dataset contracts" {
  let data = CausalDataset::new(
    [[0.0, 1.0], [1.0, 0.0], [2.0, 1.0], [1.0, 1.0], [2.0, 0.0], [0.0, 0.0]],
    [true, false, true, false, true, false],
    [2.0, 1.0, 4.0, 2.0, 5.0, 3.0],
  )
  assert_true(data.is_valid())
  assert_true(data.validate_alignment())
  assert_eq(data.n(), 6)
  assert_eq(data.p(), 2)
  assert_eq(data.treated_count(), 3)
  assert_eq(data.control_count(), 3)
  assert_eq(data.treatment_rate(), 0.5)
  assert_eq(data.outcome_range(), (1.0, 5.0))
  assert_eq(data.column(0), [0.0, 1.0, 2.0, 1.0, 2.0, 0.0])
  assert_eq(data.outcome_mean(true), 11.0 / 3.0)
  let selected = data.select([0, 2, 4])
  assert_eq(selected.n(), 3)
  assert_eq(selected.treated_count(), 3)
}

///|
test "acceptance matrix covers preprocessing boundaries" {
  let values = [1.0, 2.0, 100.0, 4.0, 5.0]
  let covariates = [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]]
  let summary = summarize_covariates(covariates)
  assert_eq(summary.columns, 2)
  assert_eq(summary.rows, 3)
  assert_eq(summary.means, [2.0, 4.0])
  assert_eq(column_values(covariates, 1), [2.0, 4.0, 6.0])
  assert_eq(column_values(covariates, 8), [])
  let centered = center_covariates(covariates)
  assert_eq(centered[0], [-1.0, -2.0])
  let scaled = scale_covariates(covariates)
  assert_true(scaled[0][0] < 0.0)
  let winsorized = winsorize(values, 0.2, 0.8)
  assert_true(winsorized[2] < 100.0)
  let ranks = rank_transform(values)
  assert_eq(ranks[0], 0.0)
  assert_eq(ranks[2], 0.8)
  let expanded = polynomial_features([[2.0, 3.0]], 3)
  assert_eq(expanded[0], [2.0, 3.0, 4.0, 9.0, 8.0, 27.0])
  let interactions = interaction_features([[2.0, 3.0]])
  assert_eq(interactions[0], [2.0, 3.0, 6.0])
}

///|
test "acceptance matrix validates graph operations" {
  let nodes = ["A", "B", "C", "D", "E"]
  let graph = CausalGraph::new(nodes)
    .add_edge("A", "B")
    .add_edge("B", "C")
    .add_edge("A", "D")
    .add_edge("D", "E")
  assert_eq(graph.edge_count(), 4)
  assert_true(graph.has_edge("A", "B"))
  assert_eq(graph.topological_order(), ["A", "B", "D", "C", "E"])
  assert_eq(graph.parents("C"), ["B"])
  assert_eq(graph.children("A"), ["B", "D"])
  assert_eq(graph.ancestors("E"), ["D", "A"])
  let cycle = graph.add_edge("C", "A")
  assert_eq(cycle.topological_order(), [])
  assert_eq(graph.backdoor_candidates("B", "C"), ["A"])
}

///|
test "acceptance matrix matches and balances covariates" {
  let covariates = [[0.0], [1.0], [2.0], [0.1], [1.1], [2.1]]
  let treatment = [true, true, true, false, false, false]
  let outcomes = [5.0, 6.0, 7.0, 3.0, 4.0, 5.0]
  let names = ["age"]
  let metrics = balance_table(covariates, treatment, names)
  assert_eq(metrics.length(), 1)
  assert_true(maximum_absolute_smd(metrics) >= 0.0)
  assert_eq(balanced_covariate_count(metrics), 0)
  let weighted = weighted_smd(covariates.map(fn(row) { row[0] }), treatment, [
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
  ])
  assert_true(weighted.abs() >= 0.0)
  let pairs = nearest_neighbor_match(covariates, treatment)
  assert_eq(pairs.length(), 3)
  let matched = estimate_matched_ate(outcomes, pairs)
  assert_true(matched.estimate > 1.0)
  assert_true(matching_mean_distance(pairs) >= 0.0)
}

///|
test "acceptance matrix covers outcome and propensity models" {
  let x = [[0.0], [1.0], [2.0], [3.0], [4.0], [5.0]]
  let y = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
  let linear = fit_linear_outcome_model(x, y)
  let predictions = predict_outcomes(linear, x)
  assert_true(r_squared(y, predictions) > 0.98)
  assert_true(linear.has_converged())
  assert_true(linear.coefficient_l2() >= 0.0)
  let residual = residuals(linear, x, y)
  assert_eq(residual.length(), 6)
  let treatment = [false, false, false, true, true, true]
  let logistic = fit_logistic_regression(
    x,
    treatment,
    learning_rate=0.1,
    max_iterations=300,
  )
  let scores = predict_propensity(logistic, x)
  assert_true(scores[0] < scores[5])
  assert_true(scores.all(fn(value) { value > 0.0 && value < 1.0 }))
  assert_true(log_loss(scores, treatment) < 0.7)
  assert_true(classification_accuracy(scores, treatment) >= 0.5)
}

///|
test "acceptance matrix covers estimators and uncertainty" {
  let treatment = [true, true, true, false, false, false]
  let outcomes = [4.0, 5.0, 6.0, 2.0, 3.0, 4.0]
  let propensity = [0.7, 0.7, 0.7, 0.3, 0.3, 0.3]
  let ipw = estimate_ipw_ate(outcomes, treatment, propensity)
  let att = estimate_ipw_att(outcomes, treatment, propensity)
  let atc = estimate_atc(outcomes, treatment, propensity)
  let overlap = estimate_overlap_ate(outcomes, treatment, propensity)
  assert_true(ipw.estimate > 1.0 && ipw.estimate < 3.0)
  assert_true(att.estimate > 1.0 && att.estimate < 3.0)
  assert_true(atc.estimate > 1.0 && atc.estimate < 3.0)
  assert_true(overlap.estimate > 1.0 && overlap.estimate < 3.0)
  let gcomp = g_computation_ate([4.0, 5.0, 6.0], [2.0, 3.0, 4.0])
  assert_true(gcomp.estimate > 1.5 && gcomp.estimate < 2.5)
  assert_true(ipw.effective_sample_size > 0.0)
  let flags = audit_estimate(ipw)
  assert_eq(flags.length(), 5)
  let bootstrap = bootstrap_ipw_ate(outcomes, treatment, propensity, 24, 7UL)
  assert_eq(bootstrap.replicates, 24)
  assert_true(bootstrap.standard_error >= 0.0)
  assert_true(bootstrap_percentile_interval([1.0, 2.0, 3.0, 4.0], 0.8).0 <= 3.0)
}

///|
test "acceptance matrix covers heterogeneity and policy" {
  let treatment = [true, true, true, false, false, false]
  let outcomes = [5.0, 6.0, 7.0, 2.0, 3.0, 4.0]
  let covariate = [0.1, 0.2, 0.9, 0.1, 0.8, 0.9]
  let effects = threshold_subgroup_effects(covariate, outcomes, treatment, 0.5)
  assert_eq(effects.length(), 2)
  let quantiles = quantile_subgroup_effects(covariate, outcomes, treatment, 3)
  assert_eq(quantiles.length(), 3)
  let cate = conditional_average_treatment_effect([2.0, 4.0], [1.0, 1.5])
  assert_eq(cate, [1.0, 2.5])
  let policy = policy_value(outcomes, treatment, [
    1.0, -1.0, 0.5, -1.0, 1.0, -0.5,
  ])
  assert_true(policy >= 0.0)
  let curve = uplift_curve(outcomes, treatment, covariate, 3)
  assert_eq(curve.length(), 3)
  assert_true(uplift_auc(outcomes, treatment, covariate) >= 0.0)
  assert_true(policy_treatment_rate([1.0, -1.0, 0.0], 0.0) > 0.0)
}

///|
test "acceptance matrix covers missing data contracts" {
  let nan = (-1.0).sqrt()
  let covariates = [[1.0, 2.0], [nan, 3.0], [4.0, nan]]
  let (clean, summary) = impute_non_finite(covariates)
  assert_eq(summary.imputed_values, 2)
  assert_eq(summary.rows, 3)
  assert_eq(summary.columns, 2)
  assert_eq(clean[1][0], 2.5)
  assert_eq(clean[2][1], 2.5)
  assert_eq(complete_case_indices(covariates), [0])
  assert_eq(missingness_rate(covariates), 2.0 / 6.0)
  assert_eq(missingness_indicator(covariates)[1], [1.0, 0.0])
  assert_eq(mean_impute_vector([1.0, nan, 3.0]), [1.0, 2.0, 3.0])
}

///|
test "acceptance matrix covers panel and survival" {
  let outcomes = [1.0, 2.0, 5.0, 7.0]
  let groups = [false, false, true, true]
  let periods = [false, true, false, true]
  let did = difference_in_differences(outcomes, groups, periods)
  assert_eq(did.estimate, 1.0)
  let event = event_study(outcomes, groups, [-1, 0, -1, 0])
  assert_eq(event.length(), 2)
  assert_eq(pretrend_difference(outcomes, groups, [-1, -1, 0, 0]), -1.5)
  let donors = synthetic_control_weights([1.0, 2.0], [[1.0, 2.0], [2.0, 3.0]])
  assert_true((donors[0] + donors[1] - 1.0).abs() < 0.0001)
  let curve = kaplan_meier([1.0, 2.0, 3.0, 4.0], [true, false, true, true])
  assert_eq(curve.length(), 3)
  assert_true(survival_at(curve, 2.0) < 1.0)
  assert_true(restricted_mean_survival(curve, 4.0) > 0.0)
  assert_true(hazard_ratio_from_curves(curve, curve, 3.0) > 0.0)
}

///|
test "acceptance matrix covers multi arm and dose response" {
  let arms = [0, 1, 2, 0, 1, 2]
  let outcomes = [1.0, 3.0, 6.0, 2.0, 4.0, 8.0]
  assert_eq(arm_levels(arms), [0, 1, 2])
  let estimates = arm_estimates(outcomes, arms, 0)
  assert_eq(estimates.length(), 3)
  assert_eq(estimates[0].mean_outcome, 1.5)
  assert_eq(arm_sample_counts(arms), [(0, 2), (1, 2), (2, 2)])
  let pair = pairwise_arm_difference(outcomes, arms, 2, 0)
  assert_eq(pair.estimate, 5.5)
  let weighted = weighted_arm_means(outcomes, arms, [
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
  ])
  assert_eq(weighted[1].mean_outcome, 3.5)
  let dose = dose_response_curve([0.0, 1.0, 2.0, 3.0], [1.0, 3.0, 6.0, 8.0], 2)
  assert_eq(dose.length(), 2)
  let gps = generalized_propensity_scores([0.0, 1.0, 2.0], [[0.0], [1.0], [2.0]])
  assert_eq(gps.length(), 3)
  assert_true(gps[0] > 0.0)
  assert_eq(overlap_weights([true, false], [0.25, 0.75]), [0.75, 0.75])
}

///|
test "acceptance matrix covers audit and diagnostics" {
  let data = CausalDataset::new(
    [[0.0], [1.0], [1.0], [0.0]],
    [true, false, true, false],
    [3.0, 1.0, 5.0, 2.0],
  )
  let report = audit_dataset(data)
  assert_true(report.valid_dataset)
  assert_true(report.adequate_sample)
  assert_true(report.has_treatment_variation)
  assert_true(report.has_outcome_variation)
  assert_true(audit_passes(report))
  let scores = [0.2, 0.8, 0.3, 0.7]
  let overlap = overlap_report(scores, data.treatment)
  assert_eq(overlap.overlap_fraction, 0.0)
  assert_eq(positivity_violations([0.0, 0.5, 1.0], 0.01), 2)
  assert_eq(propensity_histogram(scores, 2), [2, 2])
  let calibration = propensity_calibration(scores, data.treatment, 2)
  assert_eq(calibration.length(), 2)
  assert_eq(standardized_residuals([1.0, 2.0], [1.0, 3.0]).length(), 2)
  assert_true(root_mean_squared_error([1.0, 2.0], [1.0, 3.0]) > 0.0)
}

///|
test "acceptance matrix benchmark is reproducible" {
  let config = SyntheticConfig::new(64, 3, 1.5, 20260818UL)
  let first = simulate_observational_data(config)
  let second = simulate_observational_data(config)
  assert_eq(first.dataset.outcome, second.dataset.outcome)
  assert_eq(first.dataset.treatment, second.dataset.treatment)
  assert_eq(first.dataset.covariates, second.dataset.covariates)
  let report = run_synthetic_benchmark(config)
  assert_eq(report.sample_size, 64)
  assert_eq(report.true_ate, 1.5)
  assert_true(is_finite(report.absolute_error))
  assert_true(report.standard_error >= 0.0)
  assert_true(report.effective_sample_size > 0.0)
  assert_true(report.propensity_auc >= 0.0 && report.propensity_auc <= 1.0)
}

///|
test "acceptance matrix protects numeric boundaries" {
  assert_true(sigmoid(-100.0) < 1.0e-30)
  assert_true(sigmoid(100.0) > 1.0 - 1.0e-12)
  assert_eq(rosenbaum_bounds_gamma(0.05, 1.0), 0.05)
  assert_eq(e_value(1.0), 1.0)
  assert_eq(variance([], sample=false), 0.0)
  assert_eq(std_dev([], sample=false), 0.0)
  assert_eq(sum([]), 0.0)
  assert_eq(mean([]), 0.0)
  assert_eq(quantile([], 0.5), 0.0)
  assert_eq(effective_sample_size([1.0, 1.0, 1.0]), 3.0)
  assert_eq(effective_sample_size([0.0, 0.0]), 0.0)
  assert_true(is_finite(1.0))
  assert_true(!is_finite((-1.0).sqrt()))
}

///|
test "acceptance matrix linear algebra remains stable" {
  let matrix = [[4.0, 1.0], [1.0, 3.0]]
  let inverse = matrix_inverse(matrix)
  let product = matrix_multiply(matrix, inverse)
  assert_true((product[0][0] - 1.0).abs() < 0.0001)
  assert_true((product[1][1] - 1.0).abs() < 0.0001)
  assert_true(product[0][1].abs() < 0.0001)
  let solved = solve_linear_system(matrix, [9.0, 7.0])
  assert_true((solved[0] - 1.8181).abs() < 0.01)
  assert_true((solved[1] - 1.7272).abs() < 0.01)
  assert_eq(identity_matrix(3)[1][1], 1.0)
  assert_eq(matrix_transpose([[1.0, 2.0]])[1][0], 2.0)
  assert_eq(dot([1.0, 2.0], [3.0, 4.0]), 11.0)
}

///|
test "acceptance matrix formula and sampling compose" {
  let covariates = [[1.0], [2.0], [3.0]]
  let design = make_design_matrix(covariates, true)
  let predicted = rowwise_linear_predict(design.values, [1.0, 2.0])
  assert_eq(predicted, [3.0, 5.0, 7.0])
  assert_eq(append_column([[1.0], [2.0]], [3.0, 4.0]), [[1.0, 3.0], [2.0, 4.0]])
  assert_eq(add_treatment_interactions(covariates, [true, false, true])[1], [
    2.0, 0.0,
  ])
  let shuffled = shuffled_indices(8, 7UL)
  assert_eq(shuffled.length(), 8)
  assert_true(shuffled.contains(0) && shuffled.contains(7))
  assert_eq(systematic_sample(10, 3, 0UL).length(), 3)
  let bootstrap = bootstrap_counts(5, 4, 13UL)
  assert_eq(bootstrap.length(), 4)
  assert_eq(bootstrap[0].length(), 5)
  assert_eq(delete_one_mean([1.0, 2.0, 3.0]), [2.5, 2.0, 1.5])
  assert_eq(influence_function([1.0, 2.0, 3.0]), [-1.0, 0.0, 1.0])
}
