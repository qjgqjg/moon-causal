///|
test "empty vector boundaries" {
  assert_eq(sum([]), 0.0)
  assert_eq(mean([]), 0.0)
  assert_eq(variance([]), 0.0)
  assert_eq(std_dev([]), 0.0)
  assert_eq(dot([], []), 0.0)
  assert_eq(l1_norm([]), 0.0)
  assert_eq(l2_norm([]), 0.0)
  assert_eq(weighted_mean([], []), 0.0)
  assert_eq(weighted_variance([], []), 0.0)
  assert_eq(effective_sample_size([]), 0.0)
  assert_eq(quantile([], 0.1), 0.0)
  assert_eq(column_values([], 0), [])
  assert_eq(polynomial_features([], 2), [])
  assert_eq(interaction_features([]), [])
  assert_eq(make_folds(0, 5, 1UL).length(), 0)
  assert_eq(shuffled_indices(0, 1UL), [])
  assert_eq(systematic_sample(0, 3, 1UL), [])
  assert_eq(kaplan_meier([], []).length(), 0)
}

///|
test "mismatched vector boundaries" {
  assert_eq(dot([1.0], [1.0, 2.0]), 0.0)
  assert_eq(add_vectors([1.0, 2.0], [3.0]).length(), 1)
  assert_eq(subtract_vectors([1.0, 2.0], []).length(), 0)
  assert_eq(weighted_mean([1.0], [1.0, 2.0]), 0.0)
  assert_eq(weighted_variance([1.0], [1.0, 2.0]), 0.0)
  assert_eq(covariance([1.0], [2.0]).to_string(), "0")
  assert_eq(correlation([1.0], [2.0]), 0.0)
  assert_eq(root_mean_squared_error_raw([1.0], []), 0.0)
  assert_eq(mean_absolute_error_raw([1.0], []), 0.0)
  assert_eq(auc_score([0.1], []), 0.5)
  assert_eq(log_loss([0.1], []).to_string(), "0")
  assert_eq(brier_score([0.1], []).to_string(), "0")
  assert_eq(resample_by_weights([1.0], [1.0, 2.0], 5, 1UL), [])
  assert_eq(sample_variance([1.0], [1.0, 2.0]), 0.0)
}

///|
test "probability and weight boundaries" {
  let clipped = trim_propensity_scores(
    [-100.0, 0.0, 0.5, 1.0, 100.0],
    lower=0.05,
    upper=0.95,
  )
  assert_eq(clipped, [0.05, 0.05, 0.5, 0.95, 0.95])
  let overlap = overlap_weights([true, false, true, false], [0.0, 0.0, 1.0, 1.0])
  for value in overlap {
    assert_true(value > 0.0)
  }
  let entropy = entropy_like_weights([0.0, 1.0], [true, false])
  for value in entropy {
    assert_true(is_finite(value))
  }
  assert_eq(truncate_weights([1.0, 2.0, 3.0], 2.0), [1.0, 2.0, 2.0])
  assert_eq(normalized_weights([1.0, 1.0]), [1.0, 1.0])
  assert_eq(normalized_weights([]), [])
  assert_close(probability_from_odds(odds_from_probability(0.3)), 0.3, 1.0e-8)
  assert_close(safe_logit(0.5), 0.0, 1.0e-12)
  assert_true(propensity_histogram([0.0, 0.5, 1.0], 2)[0] >= 1)
  assert_eq(stratum_counts([0, 1, 2, -1], 2), [1, 1])
}

///|
test "stratification boundaries" {
  assert_eq(stratify_by_score([], 5), [])
  assert_eq(stratify_by_score([0.5, 0.5, 0.5], 3), [0, 0, 0])
  assert_eq(stratify_by_score([0.0, 1.0], 1), [0, 0])
  assert_eq(stratify_by_score([0.0, 1.0], 0), [0, 0])
  let strata = stratify_by_score([0.0, 0.25, 0.5, 0.75, 1.0], 4)
  assert_eq(strata[0], 0)
  assert_eq(strata[4], 3)
  assert_eq(stratum_counts(strata, 4).length(), 4)
  let effects = stratified_effects(
    [1.0, 2.0, 3.0, 4.0],
    [true, false, true, false],
    [0, 0, 1, 1],
  )
  assert_eq(effects.length(), 2)
  assert_close(aggregate_stratified_effects(effects).estimate, -1.0, 1.0e-12)
}

///|
test "matrix boundaries" {
  assert_eq(matrix_zeros(0, 3), [])
  assert_eq(matrix_zeros(3, 0), [[], [], []])
  assert_eq(identity_matrix(0), [])
  assert_eq(matrix_transpose([[], []]), [])
  assert_eq(matrix_multiply([[1.0]], [[], []]), [])
  assert_eq(matrix_vector_multiply([], [1.0]), [])
  assert_eq(solve_linear_system([], []), [])
  assert_eq(solve_linear_system([[1.0, 2.0]], [1.0]), [])
  assert_eq(matrix_inverse([[0.0]]), [])
  assert_eq(covariance_matrix([]), [])
  assert_eq(weighted_covariance_matrix([], []), [])
  assert_eq(add_ridge([[]], 1.0), [[]])
}

///|
test "matching boundaries" {
  assert_eq(nearest_neighbor_match([], []).length(), 0)
  assert_eq(nearest_neighbor_match([[1.0]], [true]).length(), 0)
  assert_eq(nearest_neighbor_match([[1.0]], [false]).length(), 0)
  assert_eq(matching_effects([], []), [])
  assert_eq(matching_mean_distance([]), 0.0)
  assert_eq(estimate_matched_ate([], []).sample_size, 0)
  let pairs = nearest_neighbor_match(
    [[0.0], [0.0], [0.0]],
    [true, false, false],
    replace=true,
  )
  assert_eq(pairs.length(), 1)
  assert_close(pairs[0].distance, 0.0, 1.0e-12)
  assert_eq(
    nearest_neighbor_match([[0.0], [10.0]], [true, false], caliper=0.1).length(),
    0,
  )
}

///|
test "balance boundaries" {
  assert_eq(balance_table([], [], []).length(), 0)
  let metrics = balance_table([[1.0], [1.0]], [true, false], ["x"])
  assert_eq(metrics.length(), 1)
  assert_close(metrics[0].standardized_difference, 0.0, 1.0e-12)
  assert_true(metrics[0].balanced)
  assert_eq(maximum_absolute_smd([]), 0.0)
  assert_eq(balanced_covariate_count([]), 0)
  assert_close(
    weighted_smd([1.0, 2.0], [true, false], [0.0, 0.0]),
    0.0,
    1.0e-12,
  )
}

///|
test "graph boundaries" {
  let graph = CausalGraph::new([])
  assert_eq(graph.topological_order(), [])
  assert_eq(graph.parents("missing"), [])
  assert_eq(graph.children("missing"), [])
  let single = CausalGraph::new(["x"]).add_edge("x", "x")
  assert_eq(single.edge_count(), 0)
  let duplicate = CausalGraph::new(["x", "y"])
    .add_edge("x", "y")
    .add_edge("x", "y")
  assert_eq(duplicate.edge_count(), 1)
}

///|
test "uncertainty boundaries" {
  let all_treated = bootstrap_difference_in_means(
    [1.0, 2.0],
    [true, true],
    replicates=20,
    seed=1UL,
  )
  assert_eq(all_treated.successful_replicates, 0)
  let all_control = bootstrap_difference_in_means(
    [1.0, 2.0],
    [false, false],
    replicates=20,
    seed=1UL,
  )
  assert_eq(all_control.successful_replicates, 0)
  assert_eq(permutation_p_value([], [], permutations=0, seed=1UL), 1.0)
  assert_eq(jackknife_influence([], []), [])
  let curve = sensitivity_curve(0.05, 0.0, 3)
  assert_eq(curve.length(), 3)
  assert_close(curve[0].gamma, 1.0, 1.0e-12)
}

///|
test "mediation boundaries" {
  let empty = mediation_decomposition([], [], [], [])
  assert_eq(empty.total_effect, 0.0)
  assert_true(mediation_balance_check(empty, 0.0))
  let result = mediation_decomposition([4.0, 5.0], [3.0, 4.0], [2.0, 2.0], [
    1.0, 1.0,
  ])
  assert_close(result.direct_effect, 2.5, 1.0e-12)
  assert_close(result.indirect_effect, 1.0, 1.0e-12)
  assert_close(mediation_sensitivity(result, 1.0).indirect_effect, 0.0, 1.0e-12)
}

///|
test "continuous treatment boundaries" {
  assert_eq(dose_response_curve([], [], 3).length(), 3)
  let flat = dose_response_curve([1.0, 1.0], [2.0, 4.0], 2)
  assert_eq(flat.length(), 2)
  assert_eq(generalized_propensity_scores([], []).length(), 0)
  assert_eq(overlap_weights([], []), [])
  assert_eq(dose_weighted_mean([], [], []), 0.0)
}

///|
test "survival boundaries" {
  let curve = kaplan_meier([1.0, 2.0, 3.0], [true, false, true])
  assert_eq(curve.length(), 2)
  assert_close(survival_at(curve, 0.0), 1.0, 1.0e-12)
  assert_true(survival_at(curve, 5.0) < 1.0)
  assert_true(restricted_mean_survival(curve, 3.0) > 0.0)
  assert_eq(restricted_mean_survival([], 5.0), 5.0)
  assert_close(hazard_ratio_from_curves(curve, curve, 3.0), 1.0, 1.0e-12)
}
