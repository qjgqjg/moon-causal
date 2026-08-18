///|
test "kaplan meier step function" {
  let curve = kaplan_meier([1.0, 2.0, 3.0, 4.0], [true, true, false, true])
  assert_eq(curve.length(), 3)
  assert_close(curve[0].survival, 0.75, 1.0e-12)
  assert_true(curve[1].survival < curve[0].survival)
  assert_true(curve[2].survival < curve[1].survival)
  assert_true(curve[0].at_risk >= curve[0].events)
}

///|
test "survival monotonicity" {
  let curve = kaplan_meier([1.0, 2.0, 3.0, 4.0, 5.0], [
    true, false, true, false, true,
  ])
  for i in 1..<curve.length() {
    assert_true(curve[i].survival <= curve[i - 1].survival)
  }
  for time in [0.0, 1.0, 2.5, 4.5, 10.0] {
    let survival = survival_at(curve, time)
    assert_true(survival >= 0.0 && survival <= 1.0)
  }
}

///|
test "survival area bounds" {
  let curve = kaplan_meier([1.0, 2.0, 3.0], [true, true, true])
  for horizon in [0.0, 1.0, 2.0, 3.0, 10.0] {
    let area = restricted_mean_survival(curve, horizon)
    assert_true(area >= 0.0 && area <= horizon + 1.0e-12)
  }
}
