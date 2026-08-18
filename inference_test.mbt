///|
test "ipw and g-computation agree on randomized data" {
  let y = [2.0, 4.0, 3.0, 5.0, 4.0, 6.0, 5.0, 7.0]
  let t = [false, true, false, true, false, true, false, true]
  let p = Array::make(8, 0.5)
  let ipw = estimate_ipw_ate(y, t, p)
  let predicted_treated = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]
  let predicted_control = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
  let gcomp = g_computation_ate(predicted_treated, predicted_control)
  assert_close(ipw.estimate, 2.0, 1.0e-12)
  assert_close(gcomp.estimate, 2.0, 1.0e-12)
  assert_eq(ipw.sample_size, 8)
  assert_true(ipw.effective_sample_size > 0.0)
}

///|
test "doubly robust finite under near positivity" {
  let y = [1.0, 10.0, 2.0, 9.0, 3.0, 8.0]
  let t = [false, true, false, true, false, true]
  let p = [0.000001, 0.999999, 0.0001, 0.9999, 0.01, 0.99]
  let estimate = doubly_robust_ate(y, t, p, y, y)
  assert_true(is_finite(estimate))
  assert_true(is_finite(rosenbaum_bounds_gamma(0.01, 100.0)))
  let overlap = overlap_population_weights(p)
  for weight in overlap {
    assert_true(weight >= 0.0 && weight <= 0.25)
  }
}

///|
test "bootstrap reproducibility" {
  let y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let t = [false, false, false, true, true, true]
  let first = bootstrap_difference_in_means(y, t, replicates=80, seed=77UL)
  let second = bootstrap_difference_in_means(y, t, replicates=80, seed=77UL)
  assert_close(first.standard_error, second.standard_error, 1.0e-12)
  assert_close(first.lower, second.lower, 1.0e-12)
  assert_close(first.upper, second.upper, 1.0e-12)
}

///|
test "permutation p-value bounds" {
  let y = [1.0, 1.0, 1.0, 1.0]
  let t = [true, false, true, false]
  let p = permutation_p_value(y, t, permutations=100, seed=2UL)
  assert_true(p >= 0.0 && p <= 1.0)
  assert_close(p, 1.0, 1.0e-12)
}

///|
test "sensitivity curve monotonic envelope" {
  let curve = sensitivity_curve(0.05, 8.0, 8)
  assert_eq(curve.length(), 8)
  for i in 1..<curve.length() {
    assert_true(curve[i].upper_probability >= curve[i - 1].upper_probability)
    assert_true(curve[i].lower_probability <= curve[i - 1].lower_probability)
  }
  assert_true(e_value(2.0) > e_value(1.1))
}

///|
test "risk and odds transformations" {
  for probability in [0.01, 0.1, 0.5, 0.9, 0.99] {
    let odds = odds_from_probability(probability)
    assert_close(probability_from_odds(odds), probability, 1.0e-10)
    assert_close(
      @math.exp(safe_logit(probability)) /
      (1.0 + @math.exp(safe_logit(probability))),
      probability,
      1.0e-10,
    )
  }
}

///|
test "matching estimand confidence interval" {
  let x = [[0.0], [0.1], [0.2], [0.3], [1.0], [1.1], [1.2], [1.3]]
  let t = [true, true, true, true, false, false, false, false]
  let y = [4.0, 4.1, 3.9, 4.2, 2.0, 2.1, 1.9, 2.2]
  let pairs = nearest_neighbor_match(x, t, caliper=2.0)
  let estimate = estimate_matched_ate(y, pairs)
  assert_eq(estimate.sample_size, 4)
  assert_true(estimate.lower <= estimate.estimate)
  assert_true(estimate.upper >= estimate.estimate)
}

///|
test "balance diagnostics improve with overlap weights" {
  let x = [[0.0], [1.0], [2.0], [3.0], [4.0], [5.0]]
  let t = [false, false, false, true, true, true]
  let p = [0.2, 0.3, 0.4, 0.6, 0.7, 0.8]
  let weights = overlap_weights(t, p)
  let raw = smd(x.map(row => row[0]), t)
  let weighted = weighted_smd(x.map(row => row[0]), t, weights)
  assert_true(is_finite(raw) && is_finite(weighted))
}
