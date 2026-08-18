///|
test "end to end observational study" {
  let data = CausalDataset::new(
    [
      [18.0, 0.0],
      [21.0, 1.0],
      [25.0, 0.0],
      [29.0, 1.0],
      [34.0, 0.0],
      [38.0, 1.0],
      [42.0, 0.0],
      [47.0, 1.0],
      [53.0, 0.0],
      [58.0, 1.0],
      [64.0, 0.0],
      [70.0, 1.0],
    ],
    [
      false, true, false, true, false, true, false, true, false, true, false, true,
    ],
    [3.0, 5.0, 4.0, 7.0, 6.0, 8.0, 7.0, 10.0, 8.0, 11.0, 9.0, 12.0],
  )
  let analysis = run_observational_analysis(data)
  assert_true(analysis.dataset.validate_alignment())
  assert_eq(analysis.dataset.n(), 12)
  assert_eq(analysis.propensity_scores.length(), 12)
  assert_eq(analysis.balance.length(), 2)
  assert_true(is_finite(analysis.ipw_estimate.estimate))
  assert_true(is_finite(analysis_quality_score(analysis)))
  let doubly_robust = doubly_robust_from_dataset(analysis)
  assert_true(is_finite(doubly_robust.estimate))
  let filtered = data.filter_by_score(analysis.propensity_scores, 0.05, 0.95)
  assert_true(filtered.n() <= data.n())
}

///|
test "end to end policy analysis" {
  let x = [[0.0], [1.0], [2.0], [3.0], [4.0], [5.0], [6.0], [7.0]]
  let t = [false, false, true, true, false, true, false, true]
  let y = [1.0, 2.0, 5.0, 6.0, 3.0, 7.0, 4.0, 8.0]
  let propensity_model = fit_logistic_regression(x, t, max_iterations=300)
  let propensity = predict_propensity(propensity_model, x)
  let outcome_model = fit_linear_outcome_model(x, y)
  let treated_x = x.map(row => [row[0] + 0.0])
  let control_x = x.map(row => [row[0] + 0.0])
  let treated_predictions = predict_outcomes(outcome_model, treated_x)
  let control_predictions = predict_outcomes(outcome_model, control_x)
  let effects = conditional_average_treatment_effect(
    treated_predictions, control_predictions,
  )
  assert_eq(effects.length(), 8)
  assert_true(is_finite(policy_value(y, t, effects)))
  assert_true(is_finite(estimate_ipw_ate(y, t, propensity).estimate))
  assert_true(is_finite(cross_validate_propensity(x, t, 4, 10UL)))
}

///|
test "end to end panel analysis" {
  let outcomes = [
    10.0, 11.0, 12.0, 14.0, 9.0, 10.0, 10.0, 11.0, 8.0, 8.0, 9.0, 9.0,
  ]
  let group = [
    true, true, true, true, false, false, false, false, false, false, false, false,
  ]
  let post = [
    false, false, true, true, false, false, true, true, false, false, true, true,
  ]
  let result = difference_in_differences(outcomes, group, post)
  assert_true(is_finite(result.estimate))
  assert_true(result.standard_error >= 0.0)
  let event = event_study(outcomes, group, [
    -2, -1, 0, 1, -2, -1, 0, 1, -2, -1, 0, 1,
  ])
  assert_eq(event.length(), 4)
  assert_true(
    is_finite(
      pretrend_difference(outcomes, group, [
        -2, -1, 0, 1, -2, -1, 0, 1, -2, -1, 0, 1,
      ]),
    ),
  )
}

///|
test "end to end survival analysis" {
  let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let events = [true, false, true, false, true, true]
  let treated_curve = kaplan_meier(times, events)
  let control_curve = kaplan_meier([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [
    false, true, false, true, false, true,
  ])
  assert_true(survival_at(treated_curve, 6.0) <= 1.0)
  assert_true(restricted_mean_survival(treated_curve, 6.0) >= 0.0)
  assert_true(
    is_finite(hazard_ratio_from_curves(treated_curve, control_curve, 4.0)),
  )
}

///|
test "end to end multi arm analysis" {
  let outcomes = [1.0, 2.0, 3.0, 4.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let arms = [0, 0, 0, 1, 1, 1, 2, 2, 2]
  let estimates = arm_estimates(outcomes, arms, 0)
  assert_eq(estimates.length(), 3)
  assert_true(estimates[1].difference_from_reference > 0.0)
  assert_true(
    estimates[2].difference_from_reference >
    estimates[1].difference_from_reference,
  )
  let pair = pairwise_arm_difference(outcomes, arms, 2, 0)
  assert_true(pair.estimate > 0.0)
}

///|
test "end to end reproducible benchmark" {
  let config = SyntheticConfig::new(200, 4, 1.75, 20260818UL)
  let first = run_synthetic_benchmark(config)
  let second = run_synthetic_benchmark(config)
  assert_close(first.estimated_ate, second.estimated_ate, 1.0e-12)
  assert_close(first.propensity_auc, second.propensity_auc, 1.0e-12)
  assert_true(first.absolute_error >= 0.0)
  assert_true(first.standard_error >= 0.0)
}

///|
test "end to end uncertainty report" {
  let y = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
  let t = [false, false, false, false, true, true, true, true]
  let p = [0.4, 0.45, 0.5, 0.55, 0.55, 0.6, 0.65, 0.7]
  let summary = bootstrap_ipw_ate(y, t, p, 60, 101UL)
  let estimate = estimate_ipw_ate(y, t, p)
  assert_eq(summary.replicates, 60)
  assert_true(summary.successful_replicates > 0)
  assert_true(summary.lower <= summary.upper)
  assert_true(
    is_finite(bootstrap_bias([summary.lower, summary.upper], estimate.estimate)),
  )
}
