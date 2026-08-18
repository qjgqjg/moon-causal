///|
test "alternative estimands" {
  let y = [2.0, 4.0, 3.0, 5.0, 4.0, 6.0]
  let t = [false, true, false, true, false, true]
  let p = [0.3, 0.4, 0.5, 0.5, 0.6, 0.7]
  let atc = estimate_atc(y, t, p)
  let overlap = estimate_overlap_ate(y, t, p)
  let strata = estimate_from_strata(y, t, p, 3)
  let simple = estimate_difference_in_means(y, t)
  assert_true(is_finite(atc.estimate))
  assert_true(is_finite(overlap.estimate))
  assert_true(is_finite(strata.estimate))
  assert_close(simple.estimate, 2.0, 1.0e-12)
  assert_eq(effect_curve_by_trim(y, t, p, [0.0, 0.1, 0.2]).length(), 3)
}

///|
test "trim curve common support" {
  let y = [1.0, 2.0, 5.0, 6.0]
  let t = [false, false, true, true]
  let p = [0.01, 0.2, 0.8, 0.99]
  let curve = effect_curve_by_trim(y, t, p, [0.0, 0.1, 0.25, 0.5])
  assert_eq(curve.length(), 4)
  for estimate in curve {
    assert_true(is_finite(estimate.estimate))
  }
}

///|
test "smd profile and subgroup curve" {
  let x = [[0.0, 1.0], [1.0, 2.0], [2.0, 3.0], [3.0, 4.0]]
  let t = [false, false, true, true]
  let profile = smd_profile(x, t, ["a", "b"], 0.1)
  assert_eq(profile.before.length(), 2)
  assert_eq(profile.balanced_before, 0)
  let effects = treatment_effect_by_threshold(
    [0.0, 1.0, 2.0, 3.0],
    [1.0, 2.0, 4.0, 8.0],
    t,
    [1.0, 2.0],
  )
  assert_eq(effects.length(), 4)
}

///|
test "policy rate and uplift area" {
  let effects = [-1.0, 0.0, 0.5, 1.0, 2.0]
  assert_close(policy_treatment_rate(effects, 0.5), 0.6, 1.0e-12)
  assert_eq(policy_treatment_rate([], 0.0), 0.0)
  let outcomes = [1.0, 1.0, 3.0, 4.0, 5.0]
  let treatment = [false, false, true, true, true]
  let scores = [0.1, 0.2, 0.4, 0.8, 0.9]
  assert_true(is_finite(uplift_auc(outcomes, treatment, scores)))
}

///|
test "alternate estimate interval" {
  let estimate = estimate_difference_in_means([1.0, 2.0, 4.0, 5.0], [
    false, false, true, true,
  ])
  assert_close(estimate.estimate, 3.0, 1.0e-12)
  assert_true(estimate.lower <= 3.0)
  assert_true(estimate.upper >= 3.0)
  assert_true(estimate.sample_size == 4)
}
