///|
test "coverage core arithmetic" {
  for values in [[1.0], [1.0, 2.0], [1.0, 2.0, 3.0], [-2.0, 0.0, 4.0]] {
    assert_true(is_finite(sum(values)))
    assert_true(is_finite(mean(values)))
    assert_true(variance(values) >= 0.0)
    assert_true(std_dev(values) >= 0.0)
    assert_true(is_finite(median(values)))
    assert_true(is_finite(trimmed_mean(values, 0.1)))
  }
}

///|
test "coverage vector arithmetic" {
  let a = [1.0, 2.0, 3.0]
  let b = [3.0, 2.0, 1.0]
  assert_eq(add_vectors(a, b), [4.0, 4.0, 4.0])
  assert_eq(subtract_vectors(a, b), [-2.0, 0.0, 2.0])
  assert_eq(scale_vector(a, 2.0), [2.0, 4.0, 6.0])
  assert_close(dot(a, b), 10.0, 1.0e-12)
  assert_close(l1_norm(a), 6.0, 1.0e-12)
  assert_close(l2_norm(a), 14.0.sqrt(), 1.0e-12)
}

///|
test "coverage matrix arithmetic" {
  let a = [[1.0, 2.0], [3.0, 4.0]]
  let b = [[2.0, 0.0], [1.0, 2.0]]
  assert_eq(matrix_transpose(a), [[1.0, 3.0], [2.0, 4.0]])
  assert_eq(matrix_multiply(a, b), [[4.0, 4.0], [10.0, 8.0]])
  assert_eq(matrix_vector_multiply(a, [1.0, 1.0]), [3.0, 7.0])
  assert_eq(identity_matrix(2), [[1.0, 0.0], [0.0, 1.0]])
  assert_eq(add_ridge(a, 1.0), [[2.0, 2.0], [3.0, 5.0]])
}

///|
test "coverage propensity functions" {
  let x = [[0.0], [1.0], [2.0], [-1.0]]
  let t = [false, true, true, false]
  let scores = estimate_propensity_scores(x, [1.0])
  assert_eq(scores.length(), 4)
  for score in scores {
    assert_true(score > 0.0 && score < 1.0)
  }
  let model = fit_logistic_regression(x, t)
  let predictions = predict_propensity(model, x)
  assert_eq(predictions.length(), 4)
  assert_true(auc_score(predictions, t) >= 0.0)
  assert_true(classification_accuracy(predictions, t) >= 0.0)
}

///|
test "coverage adjustment estimators" {
  let y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let t = [false, false, false, true, true, true]
  let p = [0.3, 0.4, 0.5, 0.5, 0.6, 0.7]
  let ate = estimate_ipw_ate(y, t, p)
  let att = estimate_ipw_att(y, t, p)
  assert_true(is_finite(ate.estimate))
  assert_true(is_finite(att.estimate))
  assert_eq(ipw_weights_ate(t, p).length(), 6)
  assert_eq(ipw_weights_att(t, p).length(), 6)
  assert_eq(stabilized_ipw_weights(t, p).length(), 6)
  assert_true(weight_design_effect(stabilized_ipw_weights(t, p)) >= 1.0)
  assert_true(is_finite(g_computation_ate(y, y).estimate))
}

///|
test "coverage matching and stratification" {
  let x = [[0.0], [0.5], [1.0], [1.5], [2.0], [2.5]]
  let t = [true, true, true, false, false, false]
  let scores = [0.1, 0.3, 0.5, 0.6, 0.8, 0.9]
  let pairs = nearest_neighbor_match(x, t, replace=true)
  assert_eq(pairs.length(), 3)
  assert_eq(stratify_by_score(scores, 3).length(), 6)
  let table = balance_table(x, t, ["x"])
  assert_eq(table.length(), 1)
  assert_true(maximum_absolute_smd(table) >= 0.0)
  assert_true(
    stratified_effects([1.0, 2.0, 3.0, 4.0], [true, false, true, false], [
      0, 0, 1, 1,
    ]).length() ==
    2,
  )
}

///|
test "coverage graph and formula" {
  let graph = CausalGraph::new(["x", "t", "y"])
    .add_edge("x", "t")
    .add_edge("t", "y")
  assert_eq(graph.topological_order().length(), 3)
  assert_eq(graph.parents("y"), ["t"])
  let design = make_design_matrix([[1.0, 2.0], [3.0, 4.0]], true)
  assert_eq(design.values.length(), 2)
  assert_eq(design.column_names.length(), 3)
  assert_eq(add_treatment_interactions([[1.0, 2.0]], [false]), [
    [1.0, 2.0, 0.0, 0.0],
  ])
}

///|
test "coverage uncertainty" {
  let y = [1.0, 2.0, 4.0, 5.0]
  let t = [false, false, true, true]
  let bootstrap = bootstrap_difference_in_means(y, t, replicates=30, seed=4UL)
  assert_eq(bootstrap.replicates, 30)
  assert_true(bootstrap.upper >= bootstrap.lower)
  assert_true(permutation_p_value(y, t, permutations=30, seed=4UL) <= 1.0)
  assert_eq(jackknife_influence(y, t).length(), 4)
}

///|
test "coverage preprocessing" {
  let x = [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]]
  assert_eq(center_covariates(x)[1], [0.0, 0.0])
  assert_eq(scale_covariates(x).length(), 3)
  assert_eq(interaction_features(x)[0].length(), 3)
  assert_eq(polynomial_features(x, 3)[0].length(), 6)
  assert_eq(rank_transform([3.0, 1.0, 2.0]).length(), 3)
  assert_eq(winsorize([1.0, 2.0, 100.0], 0.1, 0.9).length(), 3)
}

///|
test "coverage heterogeneous and mediation" {
  let covariate = [0.0, 1.0, 2.0, 3.0]
  let y = [1.0, 2.0, 4.0, 8.0]
  let t = [false, true, false, true]
  assert_eq(threshold_subgroup_effects(covariate, y, t, 2.0).length(), 2)
  assert_eq(quantile_subgroup_effects(covariate, y, t, 2).length(), 2)
  assert_eq(conditional_average_treatment_effect([3.0], [1.0]), [2.0])
  let m = mediation_decomposition([4.0], [3.0], [2.0], [1.0])
  assert_close(m.total_effect, 3.0, 1.0e-12)
  assert_true(mediation_balance_check(m, 1.0e-12))
}

///|
test "coverage binary and dose response" {
  let binary = threshold_outcome([0.0, 0.5, 1.0], 0.5)
  assert_eq(binary, [false, true, true])
  assert_eq(binary_outcome_to_double(binary), [0.0, 1.0, 1.0])
  let risk = risk_measure(binary, [false, true, true])
  assert_true(risk.risk_ratio > 0.0)
  assert_eq(
    dose_response_curve([0.0, 1.0, 2.0], [1.0, 2.0, 4.0], 3).length(),
    3,
  )
  assert_eq(
    generalized_propensity_scores([0.0, 1.0], [[1.0], [2.0]]).length(),
    2,
  )
}
