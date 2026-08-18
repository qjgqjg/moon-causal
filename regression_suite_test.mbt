///|
test "regression sigmoid grid" {
  for x in [-20.0, -10.0, -5.0, -1.0, 0.0, 1.0, 5.0, 10.0, 20.0] {
    let value = sigmoid(x)
    assert_true(value >= 0.0 && value <= 1.0)
    assert_close(sigmoid(-x), 1.0 - value, 1.0e-12)
  }
}

///|
test "regression probability grid" {
  for p in [0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999] {
    let bounds = rosenbaum_probability_bounds(p, 3.0)
    assert_true(bounds.lower_probability <= p)
    assert_true(bounds.upper_probability >= p)
    assert_true(bounds.amplification >= 0.0)
    assert_close(
      rosenbaum_probability_bounds(p, 1.0).lower_probability,
      p,
      1.0e-12,
    )
    assert_close(
      rosenbaum_probability_bounds(p, 1.0).upper_probability,
      p,
      1.0e-12,
    )
  }
}

///|
test "regression variance formulas" {
  let values = [2.0, 4.0, 6.0, 8.0]
  assert_close(mean(values), 5.0, 1.0e-12)
  assert_close(variance(values, sample=false), 5.0, 1.0e-12)
  assert_close(variance(values, sample=true), 20.0 / 3.0, 1.0e-12)
  assert_close(std_dev(values, sample=false), 5.0.sqrt(), 1.0e-12)
  assert_close(sum(scale_vector(values, 2.0)), 40.0, 1.0e-12)
}

///|
test "regression matrix systems" {
  for value in [1.0, 2.0, 3.0, 5.0] {
    let solution = solve_linear_system([[value]], [value * 2.0])
    assert_close(solution[0], 2.0, 1.0e-12)
  }
  let systems = [
    ([[1.0, 0.0], [0.0, 1.0]], [2.0, 3.0]),
    ([[3.0, 1.0], [1.0, 3.0]], [7.0, 5.0]),
  ]
  for system in systems {
    let result = solve_linear_system(system.0, system.1)
    assert_eq(result.length(), 2)
    assert_true(is_finite(result[0]) && is_finite(result[1]))
  }
}

///|
test "regression dataset selections" {
  let data = CausalDataset::new(
    [[0.0], [1.0], [2.0], [3.0], [4.0]],
    [false, true, false, true, false],
    [10.0, 11.0, 12.0, 13.0, 14.0],
  )
  for indices in [[], [0], [1, 3], [4, 3, 2, 1, 0], [-1, 99]] {
    let selected = data.select(indices)
    assert_true(selected.n() <= data.n())
    assert_true(selected.is_valid())
  }
}

///|
test "regression matching distances" {
  let covariates = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
  let treatment = [true, true, false, false]
  let pairs = nearest_neighbor_match(covariates, treatment, replace=false)
  assert_eq(pairs.length(), 2)
  for pair in pairs {
    assert_true(pair.treated_index < 2)
    assert_true(pair.control_index >= 2)
    assert_true(pair.distance >= 0.0)
    assert_true(pair.weight > 0.0)
  }
}

///|
test "regression fold coverage" {
  let folds = make_folds(100, 10, 42UL)
  let counts = Array::make(100, 0)
  for fold in folds {
    assert_eq(fold.validation_indices.length(), 10)
    for index in fold.validation_indices {
      counts[index] += 1
    }
  }
  for count in counts {
    assert_eq(count, 1)
  }
}

///|
test "regression simulation seeds" {
  let first = simulate_observational_data(SyntheticConfig::new(40, 2, 1.0, 1UL))
  let second = simulate_observational_data(
    SyntheticConfig::new(40, 2, 1.0, 1UL),
  )
  let third = simulate_observational_data(SyntheticConfig::new(40, 2, 1.0, 2UL))
  assert_eq(first.dataset.outcome, second.dataset.outcome)
  assert_true(first.dataset.outcome != third.dataset.outcome)
}

///|
test "regression graph ordering" {
  let graph = CausalGraph::new(["a", "b", "c", "d"])
    .add_edge("a", "b")
    .add_edge("a", "c")
    .add_edge("b", "d")
    .add_edge("c", "d")
  let order = graph.topological_order()
  assert_eq(order.length(), 4)
  assert_true(order[0] == "a")
  assert_true(order[3] == "d")
  assert_eq(graph.ancestors("d").length(), 3)
}

///|
test "regression balance thresholds" {
  let x = [[1.0], [2.0], [3.0], [4.0]]
  let t = [true, false, true, false]
  let strict = balance_table(x, t, ["x"], threshold=0.01)
  let loose = balance_table(x, t, ["x"], threshold=2.0)
  assert_true(
    balanced_covariate_count(loose) >= balanced_covariate_count(strict),
  )
}

///|
test "regression risk conversion" {
  for threshold in [-1.0, 0.0, 0.5, 1.0, 2.0] {
    let outcomes = threshold_outcome([0.0, 0.5, 1.0], threshold)
    assert_eq(outcomes.length(), 3)
    for value in binary_outcome_to_double(outcomes) {
      assert_true(value == 0.0 || value == 1.0)
    }
  }
}

///|
test "regression survival repeated events" {
  let curve = kaplan_meier([1.0, 1.0, 1.0, 2.0, 2.0], [
    true, false, true, true, false,
  ])
  assert_eq(curve.length(), 2)
  assert_true(curve[0].events == 2)
  assert_true(curve[0].survival < 1.0)
  assert_true(curve[1].survival <= curve[0].survival)
}
