///|
test "acceptance tail exercises empty and singleton inputs" {
  assert_eq(summarize_covariates([]).rows, 0)
  assert_eq(polynomial_features([], 3), [])
  assert_eq(interaction_features([]), [])
  assert_eq(shuffled_indices(-2, 1UL), [])
  assert_eq(systematic_sample(4, 0, 1UL), [])
  assert_eq(bootstrap_counts(0, 2, 1UL), [[], []])
  assert_eq(complete_case_indices([]), [])
  assert_eq(missingness_rate([]), 0.0)
  assert_eq(kaplan_meier([], []).length(), 0)
  assert_eq(restricted_mean_survival([], 4.0), 4.0)
}

///|
test "acceptance tail exercises deterministic transformations" {
  let state_a = RandomState::new(42UL)
  let state_b = RandomState::new(42UL)
  assert_eq(state_a.next_u64(), state_b.next_u64())
  assert_eq(state_a.next_u64(), state_b.next_u64())
  assert_eq(normalized_weights([1.0, 2.0, 1.0]), [0.75, 1.5, 0.75])
  assert_eq(truncate_weights([0.0, 0.5, 2.0], 0.5), [0.0, 0.5, 0.5])
  assert_eq(binary_outcome_to_double([true, false, true]), [1.0, 0.0, 1.0])
  assert_eq(threshold_outcome([0.0, 1.0, 2.0], 1.0), [false, true, true])
  assert_eq(common_support_indices([0.1, 0.5, 0.9], 0.2, 0.8), [1])
  assert_eq(binned_counts([0.0, 0.1, 0.9, 1.0], 2), [2, 2])
}

///|
test "acceptance tail checks report composition" {
  let dataset = CausalDataset::new(
    [[0.0], [1.0], [0.0], [1.0]],
    [true, false, true, false],
    [2.0, 1.0, 3.0, 2.0],
  )
  let analysis = run_observational_analysis(dataset)
  assert_true(analysis.quality_score() >= 0.0)
  assert_eq(analysis_summary_vector(analysis).length(), 5)
  let dr = doubly_robust_from_dataset(analysis)
  assert_true(dr.standard_error >= 0.0)
  let flags = audit_analysis(analysis)
  assert_true(audit_passes(flags))
}
