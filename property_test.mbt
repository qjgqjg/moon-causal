///|
test "deterministic random state" {
  let first = RandomState::new(123UL)
  let second = RandomState::new(123UL)
  for _ in 0..<20 {
    assert_close(first.uniform(), second.uniform(), 1.0e-12)
  }
  let different = RandomState::new(124UL)
  assert_true(first.uniform() != different.uniform())
  assert_true(is_finite(first.normal()))
}

///|
test "quantile monotonicity" {
  let values = [4.0, 1.0, 9.0, 2.0, 3.0]
  let mut previous = quantile(values, 0.0)
  for i in 1..<=20 {
    let current = quantile(values, i.to_double() / 20.0)
    assert_true(current >= previous)
    previous = current
  }
  assert_close(quantile(values, -1.0), 1.0, 1.0e-12)
  assert_close(quantile(values, 2.0), 9.0, 1.0e-12)
}

///|
test "weight invariants" {
  let values = [1.0, 2.0, 5.0, 10.0]
  let weights = [0.5, 1.0, 2.0, 0.5]
  let center = weighted_mean(values, weights)
  assert_true(center >= 1.0 && center <= 10.0)
  assert_true(weighted_variance(values, weights) >= 0.0)
  assert_true(effective_sample_size(weights) <= values.length().to_double())
  let normalized = normalized_weights(weights)
  assert_close(sum(normalized), values.length().to_double(), 1.0e-12)
  assert_close(weighted_mean(values, normalized), center, 1.0e-12)
}

///|
test "matrix identity invariants" {
  for size in 1..<=5 {
    let identity = identity_matrix(size)
    assert_eq(matrix_multiply(identity, identity), identity)
  }
  let matrix = [[2.0, 1.0], [1.0, 2.0]]
  let product = matrix_multiply(matrix, matrix_inverse(matrix))
  assert_close(product[0][0], 1.0, 1.0e-10)
  assert_close(product[1][1], 1.0, 1.0e-10)
  assert_close(product[0][1], 0.0, 1.0e-10)
  assert_close(product[1][0], 0.0, 1.0e-10)
}

///|
test "design matrix dimensions" {
  for width in [0, 1, 2, 3] {
    let row = Array::make(width, 1.0)
    let with_intercept = make_design_matrix([row, row], true)
    let without_intercept = make_design_matrix([row, row], false)
    assert_eq(with_intercept.values[0].length(), width + 1)
    assert_eq(without_intercept.values[0].length(), width)
  }
  let polynomial = polynomial_design([[1.0, 2.0]], 3)
  assert_eq(polynomial.values[0].length(), 6)
}

///|
test "model prediction consistency" {
  let x = [[-2.0, 1.0], [-1.0, 0.0], [1.0, 0.0], [2.0, 1.0]]
  let t = [false, false, true, true]
  let model = fit_logistic_regression(x, t, max_iterations=400)
  let direct = Array::new()
  for row in x {
    direct.push(safe_probability(sigmoid(model_predict(model, row))))
  }
  let batched = predict_propensity(model, x)
  for i in 0..<direct.length() {
    assert_close(direct[i], batched[i], 1.0e-12)
  }
  let outcome_model = fit_linear_outcome_model(x, [1.0, 2.0, 4.0, 5.0])
  assert_eq(predict_outcomes(outcome_model, x).length(), x.length())
  assert_eq(residuals(outcome_model, x, [1.0, 2.0, 4.0, 5.0]).length(), 4)
}

///|
test "cross validation partition invariants" {
  for fold_count in 2..<=5 {
    let folds = make_folds(17, fold_count, 99UL)
    assert_eq(folds.length(), fold_count)
    let seen = Array::make(17, 0)
    for fold in folds {
      assert_true(fold.validation_indices.length() > 0)
      for index in fold.validation_indices {
        seen[index] += 1
      }
      assert_eq(
        fold.training_indices.length() + fold.validation_indices.length(),
        17,
      )
    }
    for count in seen {
      assert_eq(count, 1)
    }
  }
}

///|
test "matching replacement invariant" {
  let x = [[0.0], [0.1], [0.2], [0.3]]
  let t = [true, true, true, false]
  let without = nearest_neighbor_match(x, t, replace=false)
  let with_replacement = nearest_neighbor_match(x, t, replace=true)
  assert_true(with_replacement.length() >= without.length())
  for pair in with_replacement {
    assert_eq(pair.control_index, 3)
  }
  for effect in matching_effects([4.0, 4.0, 4.0, 1.0], with_replacement) {
    assert_close(effect, 3.0, 1.0e-12)
  }
}

///|
test "stratified aggregate invariant" {
  let y = [5.0, 4.0, 7.0, 5.0, 9.0, 6.0]
  let t = [true, false, true, false, true, false]
  let effects = stratified_effects(y, t, [0, 0, 1, 1, 2, 2])
  assert_eq(effects.length(), 3)
  for effect in effects {
    assert_eq(effect.treated_count, 1)
    assert_eq(effect.control_count, 1)
  }
  assert_close(aggregate_stratified_effects(effects).estimate, 2.0, 1.0e-12)
}
