///|
test "imputation replaces only non-finite entries" {
  let x = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let (clean, summary) = impute_non_finite(x)
  assert_eq(clean, x)
  assert_eq(summary.imputed_values, 0)
  let missing = [[1.0, (-1.0).sqrt()], [(-1.0).sqrt(), 4.0], [5.0, 6.0]]
  let (filled, report) = impute_non_finite(missing)
  assert_eq(report.imputed_values, 2)
  assert_close(filled[0][1], 5.0, 1.0e-12)
  assert_close(filled[1][0], 3.0, 1.0e-12)
  assert_close(missingness_rate(missing), 2.0 / 6.0, 1.0e-12)
}

///|
test "complete case and indicators" {
  let x = [[1.0, 2.0], [(-1.0).sqrt(), 4.0], [5.0, 6.0]]
  assert_eq(complete_case_indices(x), [0, 2])
  assert_eq(missingness_indicator(x), [[0.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
  assert_eq(mean_impute_vector([1.0, (-1.0).sqrt(), 3.0]), [1.0, 2.0, 3.0])
  let data = CausalDataset::new(x, [true, false, true], [1.0, 2.0, 3.0])
  assert_eq(complete_case_dataset(data).n(), 2)
  assert_eq(impute_non_finite([]).1.rows, 0)
}

///|
test "multi-arm estimates" {
  let y = [1.0, 2.0, 4.0, 5.0, 8.0, 9.0]
  let arms = [0, 0, 1, 1, 2, 2]
  assert_eq(arm_levels(arms), [0, 1, 2])
  assert_eq(arm_sample_counts(arms), [(0, 2), (1, 2), (2, 2)])
  let estimates = arm_estimates(y, arms, 0)
  assert_eq(estimates.length(), 3)
  assert_close(estimates[0].difference_from_reference, 0.0, 1.0e-12)
  assert_close(estimates[1].difference_from_reference, 3.0, 1.0e-12)
  assert_close(estimates[2].difference_from_reference, 7.0, 1.0e-12)
  let pair = pairwise_arm_difference(y, arms, 2, 1)
  assert_close(pair.estimate, 4.0, 1.0e-12)
  let weighted = weighted_arm_means(y, arms, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  assert_eq(weighted.length(), 3)
}

///|
test "multi-arm imbalance and missing arm" {
  let estimates = arm_estimates([1.0, 2.0], [0, 0], 0)
  assert_eq(estimates.length(), 1)
  assert_eq(pairwise_arm_difference([1.0, 2.0], [0, 0], 0, 1).sample_size, 2)
  assert_eq(arm_levels([]), [])
  assert_eq(arm_sample_counts([]), [])
  assert_eq(weighted_arm_means([], [], []).length(), 0)
}

///|
test "stratified and overlap estimators" {
  let y = [2.0, 3.0, 10.0, 11.0]
  let t = [true, false, true, false]
  let strata = [0, 0, 1, 1]
  let effects = stratified_effects(y, t, strata)
  assert_eq(effects[0].count, 2)
  assert_eq(effects[1].count, 2)
  assert_close(aggregate_stratified_effects(effects).estimate, -1.0, 1.0e-12)
  let p = [0.2, 0.3, 0.7, 0.8]
  let overlap = overlap_population_weights(p)
  assert_close(overlap[0], 0.16, 1.0e-12)
  let weighted_effect = weighted_treatment_effect(
    y,
    t,
    stabilized_ipw_weights(t, p),
  )
  assert_true(is_finite(weighted_effect))
  assert_true(weighted_effect < 0.0)
}
