///|
test "dataset convenience API" {
  let dataset = CausalDataset::new(
    [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]],
    [true, false, true, false],
    [5.0, 3.0, 7.0, 4.0],
  )
  assert_eq(dataset.column(0), [1.0, 2.0, 3.0, 4.0])
  assert_close(dataset.treatment_rate(), 0.5, 1.0e-12)
  assert_close(dataset.outcome_mean(true), 6.0, 1.0e-12)
  assert_close(dataset.outcome_mean(false), 3.5, 1.0e-12)
  assert_eq(dataset.outcome_range(), (3.0, 7.0))
  assert_true(dataset.validate_alignment())
  assert_eq(dataset.with_outcome([1.0, 1.0, 1.0, 1.0]).outcome_mean(true), 1.0)
  assert_eq(dataset.with_feature_names(["age", "income"]).feature_names, [
    "age", "income",
  ])
  let filtered = dataset.filter_by_score([0.1, 0.3, 0.7, 0.9], 0.25, 0.75)
  assert_eq(filtered.n(), 2)
  assert_eq(common_support_indices([0.0, 0.2, 0.5, 1.0], 0.2, 0.5), [1, 2])
}

///|
test "estimate convenience API" {
  let estimate = Estimate::from_standard_error(1.0, 0.5, 100, 80.0, "ATE")
  assert_true(estimate.is_statistically_significant())
  assert_true(estimate.relative_precision() > 0.0)
  assert_eq(estimate.as_vector().length(), 5)
  assert_true(!estimate.is_statistically_significant(alpha=0.01))
  let nonsignificant = Estimate::from_standard_error(0.1, 1.0, 10, 10.0, "ATE")
  assert_true(!nonsignificant.is_statistically_significant())
}

///|
test "model convenience API" {
  let model = fit_linear_outcome_model([[0.0], [1.0], [2.0]], [1.0, 2.0, 3.0])
  assert_true(model.has_converged())
  assert_close(model.coefficient_l2(), model.coefficients[0].abs(), 1.0e-12)
  assert_eq(marginal_effects(model, [[1.0]]).length(), 1)
  assert_eq(make_design_matrix([[1.0, 2.0]], true).values[0], [1.0, 1.0, 2.0])
  assert_eq(make_design_matrix([[1.0, 2.0]], false).values[0], [1.0, 2.0])
  let interaction = add_treatment_interactions([[2.0, 3.0]], [true])
  assert_eq(interaction[0], [2.0, 3.0, 2.0, 3.0])
  assert_eq(append_column([[1.0], [2.0]], [3.0, 4.0])[1], [2.0, 4.0])
  assert_eq(rowwise_linear_predict([[1.0, 2.0]], [2.0, 3.0]), [8.0])
}

///|
test "sampling utilities" {
  let order = shuffled_indices(10, 7UL)
  assert_eq(order.length(), 10)
  for i in 0..<10 {
    assert_true(order.contains(i))
  }
  let systematic = systematic_sample(10, 4, 7UL)
  assert_eq(systematic.length(), 4)
  for index in systematic {
    assert_true(index >= 0 && index < 10)
  }
  let counts = bootstrap_counts(6, 3, 9UL)
  assert_eq(counts.length(), 3)
  for count in counts {
    assert_eq(sum(count.map(x => x.to_double())), 6.0)
  }
  assert_eq(influence_function([1.0, 2.0, 3.0]), [-1.0, 0.0, 1.0])
  assert_eq(delete_one_mean([1.0, 2.0, 3.0]), [2.5, 2.0, 1.5])
  let sampled = resample_by_weights([1.0, 2.0], [1.0, 0.0], 20, 11UL)
  for value in sampled {
    assert_eq(value, 1.0)
  }
}

///|
test "robust statistics" {
  assert_close(trimmed_mean([0.0, 1.0, 2.0, 100.0], 0.25), 1.5, 1.0e-12)
  assert_true(
    winsorized_mean([0.0, 1.0, 2.0, 100.0], 0.25) < mean([0.0, 1.0, 2.0, 100.0]),
  )
  assert_close(
    weighted_quantile([1.0, 2.0, 3.0], [1.0, 1.0, 1.0], 0.5),
    2.0,
    1.0e-12,
  )
  assert_eq(empirical_cdf([1.0, 2.0, 3.0], [0.0, 2.0, 4.0]), [
    0.0,
    2.0 / 3.0,
    1.0,
  ])
  assert_close(
    standardized_mean_difference_from_groups([2.0, 2.0], [1.0, 1.0]),
    0.0,
    1.0e-12,
  )
  assert_eq(outlier_indices([1.0, 1.0, 1.0], 3.0), [])
  assert_eq(robust_z_scores([2.0, 2.0]), [0.0, 0.0])
}

///|
test "panel API" {
  let outcomes = [10.0, 11.0, 12.0, 13.0, 10.0, 10.0, 11.0, 11.0]
  let groups = [true, true, true, true, false, false, false, false]
  let periods = [false, false, true, true, false, false, true, true]
  let result = difference_in_differences(outcomes, groups, periods)
  assert_close(result.treatment_change, 2.0, 1.0e-12)
  assert_close(result.control_change, 1.0, 1.0e-12)
  assert_close(result.estimate, 1.0, 1.0e-12)
  assert_true(result.standard_error >= 0.0)
  let event = event_study(outcomes, groups, [-1, 0, 1, 2, -1, 0, 1, 2])
  assert_eq(event[0].count, 2)
  assert_close(
    pretrend_difference(outcomes, groups, [-1, 0, 1, 2, -1, 0, 1, 2]),
    0.0,
    1.0e-12,
  )
}

///|
test "risk and survival API" {
  let measure = risk_measure([true, true, false, false, true, false], [
    true, true, false, false, false, true,
  ])
  assert_true(measure.risk_ratio > 0.0)
  assert_true(is_finite(log_risk_ratio(measure)))
  assert_true(relative_risk_reduction(measure).abs() <= 1.0)
  let curve = kaplan_meier([1.0, 2.0, 2.0, 4.0], [true, true, false, true])
  assert_true(curve[0].survival < 1.0)
  assert_true(curve[1].standard_error >= 0.0)
  assert_true(restricted_mean_survival(curve, 5.0) <= 5.0)
}
