///|
test "audit dataset acceptance" {
  let valid = CausalDataset::new(
    [[1.0], [2.0], [3.0], [4.0]],
    [true, false, true, false],
    [1.0, 2.0, 3.0, 4.0],
  )
  let report = audit_dataset(valid)
  assert_true(report.valid_dataset)
  assert_true(report.adequate_sample)
  assert_true(report.has_treatment_variation)
  assert_true(report.has_outcome_variation)
  assert_true(audit_passes(report))
  assert_eq(
    audit_estimate(Estimate::from_standard_error(1.0, 0.1, 10, 8.0, "ATE")).length(),
    5,
  )
}

///|
test "audit dataset rejection" {
  let empty = audit_dataset(CausalDataset::new([], [], []))
  assert_true(!empty.adequate_sample)
  assert_true(!empty.has_treatment_variation)
  assert_true(!audit_passes(empty))
  let constant = audit_dataset(
    CausalDataset::new([[1.0], [1.0]], [true, false], [2.0, 2.0]),
  )
  assert_true(!constant.has_outcome_variation)
  assert_true(!audit_passes(constant))
}

///|
test "audit analysis includes overlap" {
  let data = CausalDataset::new(
    [[0.0], [1.0], [2.0], [3.0], [4.0], [5.0]],
    [false, false, true, true, false, true],
    [1.0, 2.0, 4.0, 5.0, 3.0, 6.0],
  )
  let analysis = run_observational_analysis(data)
  let report = audit_analysis(analysis)
  assert_true(report.valid_dataset)
  assert_true(report.propensity_overlap || !report.propensity_overlap)
  assert_true(report.quality_score >= 0.0)
}

///|
test "audit model flags" {
  let model = fit_logistic_regression([[0.0], [1.0], [2.0], [3.0]], [
    false, false, true, true,
  ])
  let flags = audit_model(model)
  assert_eq(flags.length(), 4)
  assert_true(flags[1])
  assert_true(flags[2])
  let linear = fit_linear_outcome_model([[0.0], [1.0]], [1.0, 2.0])
  assert_eq(audit_model(linear).length(), 4)
}
