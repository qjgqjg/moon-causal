///|
fn feature_statistics(
  x : Array[Array[Double]],
) -> (Array[Double], Array[Double]) {
  if x.length() == 0 {
    return ([], [])
  }
  let width = x[0].length()
  let means = Array::make(width, 0.0)
  for row in x {
    for j in 0..<width {
      means[j] += row[j]
    }
  }
  for j in 0..<width {
    means[j] /= x.length().to_double()
  }
  let scales = Array::make(width, 0.0)
  for row in x {
    for j in 0..<width {
      let delta = row[j] - means[j]
      scales[j] += delta * delta
    }
  }
  for j in 0..<width {
    scales[j] = (scales[j] / x.length().to_double()).sqrt()
    if scales[j] < 1.0e-12 {
      scales[j] = 1.0
    }
  }
  (means, scales)
}

///|
fn standardized_row(
  row : Array[Double],
  means : Array[Double],
  scales : Array[Double],
) -> Array[Double] {
  let result = Array::new(capacity=row.length())
  for j in 0..<row.length() {
    result.push((row[j] - means[j]) / scales[j])
  }
  result
}

///|
fn model_linear_predict(row : Array[Double], model : ModelFit) -> Double {
  let standardized = standardized_row(
    row,
    model.feature_means,
    model.feature_scales,
  )
  model.intercept + dot(standardized, model.coefficients)
}

///|
/// Fits a binary logistic model with standardized features and L2 regularization.
pub fn fit_logistic_regression(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  learning_rate? : Double = 0.15,
  max_iterations? : Int = 600,
  l2? : Double = 1.0e-6,
) -> ModelFit {
  let n = covariates.length()
  let width = if n == 0 { 0 } else { covariates[0].length() }
  let (means, scales) = feature_statistics(covariates)
  let coefficients = Array::make(width, 0.0)
  let mut intercept = 0.0
  let mut converged = false
  let mut iterations = 0
  for iteration in 0..<max_iterations {
    let mut gradient_intercept = 0.0
    let gradients = Array::make(width, 0.0)
    for i in 0..<n {
      let row = standardized_row(covariates[i], means, scales)
      let linear = intercept + dot(row, coefficients)
      let prediction = sigmoid(linear)
      let target = if treatment[i] { 1.0 } else { 0.0 }
      let error = prediction - target
      gradient_intercept += error
      for j in 0..<width {
        gradients[j] += error * row[j]
      }
    }
    let denominator = if n == 0 { 1.0 } else { n.to_double() }
    let intercept_step = learning_rate * gradient_intercept / denominator
    intercept -= intercept_step
    let mut max_step = intercept_step.abs()
    for j in 0..<width {
      let gradient = gradients[j] / denominator + l2 * coefficients[j]
      let step = learning_rate * gradient
      coefficients[j] -= step
      if step.abs() > max_step {
        max_step = step.abs()
      }
    }
    iterations = iteration + 1
    if max_step < 1.0e-8 {
      converged = true
      break
    }
  }
  {
    coefficients,
    intercept,
    iterations,
    converged,
    loss: logistic_log_loss(
      covariates, treatment, coefficients, intercept, means, scales, l2,
    ),
    feature_means: means,
    feature_scales: scales,
  }
}

///|
fn logistic_log_loss(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  coefficients : Array[Double],
  intercept : Double,
  means : Array[Double],
  scales : Array[Double],
  l2 : Double,
) -> Double {
  if covariates.length() == 0 {
    return 0.0
  }
  let mut loss = 0.0
  for i in 0..<covariates.length() {
    let row = standardized_row(covariates[i], means, scales)
    let probability = safe_probability(
      sigmoid(intercept + dot(row, coefficients)),
    )
    let target = if treatment[i] { 1.0 } else { 0.0 }
    loss -= target * @math.ln(probability) +
      (1.0 - target) * @math.ln(1.0 - probability)
  }
  for coefficient in coefficients {
    loss += 0.5 *
      l2 *
      coefficient *
      coefficient *
      covariates.length().to_double()
  }
  loss / covariates.length().to_double()
}

///|
/// Predicts treatment probabilities from a fitted logistic model.
pub fn predict_propensity(
  model : ModelFit,
  covariates : Array[Array[Double]],
) -> Array[Double] {
  let result = Array::new(capacity=covariates.length())
  for row in covariates {
    result.push(safe_probability(sigmoid(model_linear_predict(row, model))))
  }
  result
}

///|
pub fn model_predict(model : ModelFit, row : Array[Double]) -> Double {
  model_linear_predict(row, model)
}

///|
/// Fits an ordinary or ridge-regularized linear outcome model.
pub fn fit_linear_outcome_model(
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
  ridge? : Double = 1.0e-8,
) -> ModelFit {
  let n = covariates.length()
  let width = if n == 0 { 0 } else { covariates[0].length() }
  let (means, scales) = feature_statistics(covariates)
  let design = matrix_zeros(n, width + 1)
  for i in 0..<n {
    design[i][0] = 1.0
    let row = standardized_row(covariates[i], means, scales)
    for j in 0..<width {
      design[i][j + 1] = row[j]
    }
  }
  let transposed = matrix_transpose(design)
  let normal = add_ridge(matrix_multiply(transposed, design), ridge)
  if normal.length() > 0 {
    normal[0][0] -= ridge
  }
  let rhs = matrix_vector_multiply(transposed, outcomes)
  let solution = solve_linear_system(normal, rhs)
  let coefficients = Array::new(capacity=width)
  let mut intercept = 0.0
  if solution.length() == width + 1 {
    intercept = solution[0]
    for j in 0..<width {
      coefficients.push(solution[j + 1])
    }
  } else {
    for _ in 0..<width {
      coefficients.push(0.0)
    }
  }
  let mut loss = 0.0
  for i in 0..<n {
    let error = if solution.length() == width + 1 {
      design[i][0] * intercept +
      dot(design[i][1:].to_owned(), coefficients) -
      outcomes[i]
    } else {
      -outcomes[i]
    }
    loss += error * error
  }
  {
    coefficients,
    intercept,
    iterations: 1,
    converged: solution.length() == width + 1,
    loss: if n == 0 {
      0.0
    } else {
      loss / n.to_double()
    },
    feature_means: means,
    feature_scales: scales,
  }
}

///|
pub fn predict_outcomes(
  model : ModelFit,
  covariates : Array[Array[Double]],
) -> Array[Double] {
  let result = Array::new(capacity=covariates.length())
  for row in covariates {
    result.push(model_linear_predict(row, model))
  }
  result
}

///|
pub fn residuals(
  model : ModelFit,
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
) -> Array[Double] {
  let predictions = predict_outcomes(model, covariates)
  let result = Array::new(capacity=outcomes.length())
  for i in 0..<outcomes.length() {
    result.push(outcomes[i] - predictions[i])
  }
  result
}

///|
pub fn classification_accuracy(
  predictions : Array[Double],
  treatment : Array[Bool],
  threshold? : Double = 0.5,
) -> Double {
  if predictions.length() == 0 || predictions.length() != treatment.length() {
    return 0.0
  }
  let mut correct = 0
  for i in 0..<predictions.length() {
    let predicted = predictions[i] >= threshold
    if predicted == treatment[i] {
      correct += 1
    }
  }
  correct.to_double() / predictions.length().to_double()
}
