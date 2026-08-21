///|
/// Configuration shared by deterministic first-order optimizers.
pub struct OptimizerConfig {
  learning_rate : Double
  tolerance : Double
  max_iterations : Int
  l2 : Double
  l1 : Double
  fit_intercept : Bool
}

///|
/// Iteration-level optimization trace.
pub struct OptimizationTrace {
  iterations : Int
  converged : Bool
  initial_loss : Double
  final_loss : Double
  gradient_norm : Double
  step_norm : Double
  losses : Array[Double]
}

///|
/// Generic linear optimization result with standardized feature metadata.
pub struct OptimizedLinearModel {
  coefficients : Array[Double]
  intercept : Double
  feature_means : Array[Double]
  feature_scales : Array[Double]
  trace : OptimizationTrace
}

///|
/// Binary logistic optimization result and calibration diagnostics.
pub struct LogisticOptimizationResult {
  model : OptimizedLinearModel
  probabilities : Array[Double]
  log_loss : Double
  accuracy : Double
  brier_score : Double
}

///|
/// A one-dimensional quadratic objective.
pub struct QuadraticObjective {
  center : Double
  curvature : Double
  linear_term : Double
}

///|
/// Creates a scalar quadratic objective.
pub fn quadratic_objective(
  center : Double,
  curvature : Double,
  linear_term? : Double = 0.0,
) -> QuadraticObjective {
  { center, curvature, linear_term }
}

///|
/// Summary of a backtracking line search.
pub struct LineSearchResult {
  step : Double
  objective : Double
  evaluations : Int
  accepted : Bool
}

///|
fn opt_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    sum(values) / values.length().to_double()
  }
}

///|
fn opt_dot(first : Array[Double], second : Array[Double]) -> Double {
  let n = if first.length() < second.length() {
    first.length()
  } else {
    second.length()
  }
  let mut result = 0.0
  for i in 0..<n {
    result += first[i] * second[i]
  }
  result
}

///|
fn opt_norm(values : Array[Double]) -> Double {
  let mut result = 0.0
  for value in values {
    result += value * value
  }
  result.sqrt()
}

///|
fn opt_standardization(
  matrix : Array[Array[Double]],
) -> (Array[Double], Array[Double]) {
  if matrix.length() == 0 {
    return ([], [])
  }
  let width = matrix[0].length()
  let means = Array::make(width, 0.0)
  let scales = Array::make(width, 1.0)
  for row in matrix {
    for j in 0..<width {
      means[j] += row[j]
    }
  }
  for j in 0..<width {
    means[j] /= matrix.length().to_double()
  }
  for row in matrix {
    for j in 0..<width {
      let difference = row[j] - means[j]
      scales[j] += difference * difference
    }
  }
  for j in 0..<width {
    scales[j] = (scales[j] / matrix.length().to_double()).sqrt()
    if scales[j] < 1.0e-12 {
      scales[j] = 1.0
    }
  }
  (means, scales)
}

///|
fn opt_standardize(
  row : Array[Double],
  means : Array[Double],
  scales : Array[Double],
) -> Array[Double] {
  let result = Array::new(capacity=row.length())
  for j in 0..<row.length() {
    let mean_value = if j < means.length() { means[j] } else { 0.0 }
    let scale = if j < scales.length() { scales[j] } else { 1.0 }
    result.push((row[j] - mean_value) / scale)
  }
  result
}

///|
/// Returns a conservative optimizer configuration.
pub fn default_optimizer_config() -> OptimizerConfig {
  {
    learning_rate: 0.08,
    tolerance: 1.0e-7,
    max_iterations: 800,
    l2: 1.0e-6,
    l1: 0.0,
    fit_intercept: true,
  }
}

///|
/// Creates an optimizer configuration with explicit regularization.
pub fn optimizer_config(
  learning_rate? : Double = 0.08,
  tolerance? : Double = 1.0e-7,
  max_iterations? : Int = 800,
  l2? : Double = 1.0e-6,
  l1? : Double = 0.0,
  fit_intercept? : Bool = true,
) -> OptimizerConfig {
  {
    learning_rate: if learning_rate > 0.0 {
      learning_rate
    } else {
      0.08
    },
    tolerance: if tolerance > 0.0 {
      tolerance
    } else {
      1.0e-7
    },
    max_iterations: if max_iterations > 0 {
      max_iterations
    } else {
      800
    },
    l2: if l2 >= 0.0 {
      l2
    } else {
      0.0
    },
    l1: if l1 >= 0.0 {
      l1
    } else {
      0.0
    },
    fit_intercept,
  }
}

///|
fn soft_threshold(value : Double, penalty : Double) -> Double {
  if value > penalty {
    value - penalty
  } else if value < -penalty {
    value + penalty
  } else {
    0.0
  }
}

///|
fn bounded_sigmoid(value : Double) -> Double {
  let limited = clamp(value, -40.0, 40.0)
  if limited >= 0.0 {
    1.0 / (1.0 + @math.exp(-limited))
  } else {
    let exponential = @math.exp(limited)
    exponential / (1.0 + exponential)
  }
}

///|
fn binary_target(value : Bool) -> Double {
  if value {
    1.0
  } else {
    0.0
  }
}

///|
fn logistic_loss_internal(
  matrix : Array[Array[Double]],
  labels : Array[Bool],
  coefficients : Array[Double],
  intercept : Double,
  config : OptimizerConfig,
) -> Double {
  let n = if matrix.length() < labels.length() {
    matrix.length()
  } else {
    labels.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut loss = 0.0
  for i in 0..<n {
    let probability = safe_probability(
      bounded_sigmoid(intercept + opt_dot(matrix[i], coefficients)),
    )
    let target = if labels[i] { 1.0 } else { 0.0 }
    loss -= target * @math.ln(probability) +
      (1.0 - target) * @math.ln(1.0 - probability)
  }
  for coefficient in coefficients {
    loss += 0.5 * config.l2 * coefficient * coefficient
    loss += config.l1 * coefficient.abs()
  }
  loss / n.to_double()
}

///|
fn linear_loss_internal(
  matrix : Array[Array[Double]],
  outcomes : Array[Double],
  coefficients : Array[Double],
  intercept : Double,
  config : OptimizerConfig,
) -> Double {
  let n = if matrix.length() < outcomes.length() {
    matrix.length()
  } else {
    outcomes.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut loss = 0.0
  for i in 0..<n {
    let error = intercept + opt_dot(matrix[i], coefficients) - outcomes[i]
    loss += error * error * 0.5
  }
  for coefficient in coefficients {
    loss += 0.5 * config.l2 * coefficient * coefficient +
      config.l1 * coefficient.abs()
  }
  loss / n.to_double()
}

///|
fn append_trace(
  losses : Array[Double],
  iterations : Int,
  converged : Bool,
  gradient : Double,
  step : Double,
) -> OptimizationTrace {
  let initial = if losses.length() == 0 { 0.0 } else { losses[0] }
  let last_loss = if losses.length() == 0 {
    0.0
  } else {
    losses[losses.length() - 1]
  }
  {
    iterations,
    converged,
    initial_loss: initial,
    final_loss: last_loss,
    gradient_norm: gradient,
    step_norm: step,
    losses,
  }
}

///|
/// Fits a regularized logistic model with stable batch gradient descent.
pub fn optimize_logistic(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  config? : OptimizerConfig = default_optimizer_config(),
) -> LogisticOptimizationResult {
  let n = if covariates.length() < treatment.length() {
    covariates.length()
  } else {
    treatment.length()
  }
  if n == 0 {
    let trace = append_trace([], 0, false, 0.0, 0.0)
    let model = {
      coefficients: [],
      intercept: 0.0,
      feature_means: [],
      feature_scales: [],
      trace,
    }
    return {
      model,
      probabilities: [],
      log_loss: 0.0,
      accuracy: 0.0,
      brier_score: 0.0,
    }
  }
  let (means, scales) = opt_standardization(covariates[:n].to_owned())
  let standardized : Array[Array[Double]] = Array::new(capacity=n)
  for i in 0..<n {
    standardized.push(opt_standardize(covariates[i], means, scales))
  }
  let width = standardized[0].length()
  let coefficients = Array::make(width, 0.0)
  let mut intercept = 0.0
  let losses : Array[Double] = Array::new(capacity=config.max_iterations)
  let mut converged = false
  let mut last_gradient = 0.0
  let mut last_step = 0.0
  let mut iterations = 0
  for iteration in 0..<config.max_iterations {
    let gradient = Array::make(width, 0.0)
    let mut intercept_gradient = 0.0
    for i in 0..<n {
      let probability = bounded_sigmoid(
        intercept + opt_dot(standardized[i], coefficients),
      )
      let error = probability - binary_target(treatment[i])
      if config.fit_intercept {
        intercept_gradient += error
      }
      for j in 0..<width {
        gradient[j] += error * standardized[i][j]
      }
    }
    let scale = n.to_double()
    intercept_gradient /= scale
    let mut step_norm = if config.fit_intercept {
      config.learning_rate * intercept_gradient
    } else {
      0.0
    }
    if config.fit_intercept {
      intercept -= config.learning_rate * intercept_gradient
    }
    for j in 0..<width {
      let raw = gradient[j] / scale + config.l2 * coefficients[j]
      let updated = soft_threshold(
        coefficients[j] - config.learning_rate * raw,
        config.learning_rate * config.l1,
      )
      let step = updated - coefficients[j]
      coefficients[j] = updated
      step_norm += step * step
    }
    step_norm = step_norm.sqrt()
    let current_loss = logistic_loss_internal(
      standardized, treatment, coefficients, intercept, config,
    )
    losses.push(current_loss)
    last_gradient = opt_norm(gradient) / scale + intercept_gradient.abs()
    last_step = step_norm
    iterations = iteration + 1
    if last_gradient < config.tolerance || step_norm < config.tolerance {
      converged = true
      break
    }
  }
  let trace = append_trace(
    losses, iterations, converged, last_gradient, last_step,
  )
  let model = {
    coefficients,
    intercept,
    feature_means: means,
    feature_scales: scales,
    trace,
  }
  let probabilities = optimized_predict(model, covariates[:n].to_owned())
  let mut correct = 0
  let mut brier = 0.0
  for i in 0..<n {
    let target = if treatment[i] { 1.0 } else { 0.0 }
    if (probabilities[i] >= 0.5) == treatment[i] {
      correct += 1
    }
    let error = probabilities[i] - target
    brier += error * error
  }
  {
    model,
    probabilities,
    log_loss: trace.final_loss,
    accuracy: correct.to_double() / n.to_double(),
    brier_score: brier / n.to_double(),
  }
}

///|
/// Predicts probabilities from an optimized linear model.
pub fn optimized_predict(
  model : OptimizedLinearModel,
  covariates : Array[Array[Double]],
) -> Array[Double] {
  let result = Array::new(capacity=covariates.length())
  for row in covariates {
    let standardized = opt_standardize(
      row,
      model.feature_means,
      model.feature_scales,
    )
    result.push(
      safe_probability(
        bounded_sigmoid(
          model.intercept + opt_dot(standardized, model.coefficients),
        ),
      ),
    )
  }
  result
}

///|
/// Predicts the linear response of an optimized model.
pub fn optimized_linear_predict(
  model : OptimizedLinearModel,
  covariates : Array[Array[Double]],
) -> Array[Double] {
  let result = Array::new(capacity=covariates.length())
  for row in covariates {
    let standardized = opt_standardize(
      row,
      model.feature_means,
      model.feature_scales,
    )
    result.push(model.intercept + opt_dot(standardized, model.coefficients))
  }
  result
}

///|
/// Fits a ridge linear outcome model by stable gradient descent.
pub fn optimize_ridge(
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
  config? : OptimizerConfig = optimizer_config(learning_rate=0.05, l2=1.0e-3),
) -> OptimizedLinearModel {
  let n = if covariates.length() < outcomes.length() {
    covariates.length()
  } else {
    outcomes.length()
  }
  if n == 0 {
    return {
      coefficients: [],
      intercept: 0.0,
      feature_means: [],
      feature_scales: [],
      trace: append_trace([], 0, false, 0.0, 0.0),
    }
  }
  let (means, scales) = opt_standardization(covariates[:n].to_owned())
  let design : Array[Array[Double]] = Array::new(capacity=n)
  for i in 0..<n {
    design.push(opt_standardize(covariates[i], means, scales))
  }
  let width = design[0].length()
  let coefficients = Array::make(width, 0.0)
  let mut intercept = opt_mean(outcomes[:n].to_owned())
  let losses : Array[Double] = Array::new(capacity=config.max_iterations)
  let mut converged = false
  let mut gradient_norm = 0.0
  let mut step_norm = 0.0
  let mut iterations = 0
  for iteration in 0..<config.max_iterations {
    let gradient = Array::make(width, 0.0)
    let mut intercept_gradient = 0.0
    for i in 0..<n {
      let error = intercept + opt_dot(design[i], coefficients) - outcomes[i]
      if config.fit_intercept {
        intercept_gradient += error
      }
      for j in 0..<width {
        gradient[j] += error * design[i][j]
      }
    }
    let denominator = n.to_double()
    intercept_gradient /= denominator
    if config.fit_intercept {
      intercept -= config.learning_rate * intercept_gradient
    }
    let mut local_step = if config.fit_intercept {
      config.learning_rate * intercept_gradient
    } else {
      0.0
    }
    for j in 0..<width {
      let raw = gradient[j] / denominator + config.l2 * coefficients[j]
      let updated = soft_threshold(
        coefficients[j] - config.learning_rate * raw,
        config.learning_rate * config.l1,
      )
      let step = updated - coefficients[j]
      coefficients[j] = updated
      local_step += step * step
    }
    step_norm = local_step.sqrt()
    gradient_norm = opt_norm(gradient) / denominator + intercept_gradient.abs()
    losses.push(
      linear_loss_internal(
        design,
        outcomes[:n].to_owned(),
        coefficients,
        intercept,
        config,
      ),
    )
    iterations = iteration + 1
    if gradient_norm < config.tolerance || step_norm < config.tolerance {
      converged = true
      break
    }
  }
  {
    coefficients,
    intercept,
    feature_means: means,
    feature_scales: scales,
    trace: append_trace(losses, iterations, converged, gradient_norm, step_norm),
  }
}

///|
/// Fits a lasso model with coordinate descent and soft thresholding.
pub fn optimize_lasso(
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
  l1? : Double = 0.05,
  max_iterations? : Int = 500,
  tolerance? : Double = 1.0e-7,
) -> OptimizedLinearModel {
  let n = if covariates.length() < outcomes.length() {
    covariates.length()
  } else {
    outcomes.length()
  }
  if n == 0 {
    return {
      coefficients: [],
      intercept: 0.0,
      feature_means: [],
      feature_scales: [],
      trace: append_trace([], 0, false, 0.0, 0.0),
    }
  }
  let (means, scales) = opt_standardization(covariates[:n].to_owned())
  let design : Array[Array[Double]] = Array::new(capacity=n)
  for i in 0..<n {
    design.push(opt_standardize(covariates[i], means, scales))
  }
  let width = design[0].length()
  let coefficients = Array::make(width, 0.0)
  let mut intercept = opt_mean(outcomes[:n].to_owned())
  let losses : Array[Double] = Array::new(capacity=max_iterations)
  let mut converged = false
  let mut max_change = 0.0
  let mut iterations = 0
  for iteration in 0..<max_iterations {
    let previous_intercept = intercept
    intercept = opt_mean(outcomes[:n].to_owned())
    let predictions = optimized_linear_predict(
      {
        coefficients: coefficients.copy(),
        intercept: previous_intercept,
        feature_means: means,
        feature_scales: scales,
        trace: append_trace([], 0, false, 0.0, 0.0),
      },
      covariates[:n].to_owned(),
    )
    for prediction in predictions {
      intercept += prediction * 0.0
    }
    for j in 0..<width {
      let mut numerator = 0.0
      let mut denominator = 0.0
      for i in 0..<n {
        let mut partial = intercept
        for k in 0..<width {
          if k != j {
            partial += design[i][k] * coefficients[k]
          }
        }
        numerator += design[i][j] * (outcomes[i] - partial)
        denominator += design[i][j] * design[i][j]
      }
      let updated = if denominator == 0.0 {
        0.0
      } else {
        soft_threshold(numerator / n.to_double(), l1) /
        (denominator / n.to_double())
      }
      let change = (updated - coefficients[j]).abs()
      if change > max_change {
        max_change = change
      }
      coefficients[j] = updated
    }
    let config = optimizer_config(l1~, tolerance~, max_iterations=1)
    losses.push(
      linear_loss_internal(
        design,
        outcomes[:n].to_owned(),
        coefficients,
        intercept,
        config,
      ),
    )
    iterations = iteration + 1
    if max_change < tolerance {
      converged = true
      break
    }
    max_change = 0.0
  }
  {
    coefficients,
    intercept,
    feature_means: means,
    feature_scales: scales,
    trace: append_trace(losses, iterations, converged, max_change, max_change),
  }
}

///|
/// Calculates residuals from a fitted optimized linear model.
pub fn optimized_residuals(
  model : OptimizedLinearModel,
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
) -> Array[Double] {
  let predictions = optimized_linear_predict(model, covariates)
  let n = if predictions.length() < outcomes.length() {
    predictions.length()
  } else {
    outcomes.length()
  }
  let result = Array::new(capacity=n)
  for i in 0..<n {
    result.push(outcomes[i] - predictions[i])
  }
  result
}

///|
/// Computes the coefficient of determination for an optimized model.
pub fn optimized_r_squared(
  model : OptimizedLinearModel,
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
) -> Double {
  let residuals = optimized_residuals(model, covariates, outcomes)
  let observed = outcomes[:residuals.length()].to_owned()
  if observed.length() < 2 {
    return 0.0
  }
  let average = opt_mean(observed)
  let mut total = 0.0
  for value in observed {
    total += (value - average) * (value - average)
  }
  let mut unexplained = 0.0
  for residual in residuals {
    unexplained += residual * residual
  }
  if total == 0.0 {
    0.0
  } else {
    1.0 - unexplained / total
  }
}

///|
/// Computes mean absolute error from an optimized model.
pub fn optimized_mae(
  model : OptimizedLinearModel,
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
) -> Double {
  let residuals = optimized_residuals(model, covariates, outcomes)
  if residuals.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  for residual in residuals {
    total += residual.abs()
  }
  total / residuals.length().to_double()
}

///|
/// Computes root mean squared error from an optimized model.
pub fn optimized_rmse(
  model : OptimizedLinearModel,
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
) -> Double {
  let residuals = optimized_residuals(model, covariates, outcomes)
  if residuals.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  for residual in residuals {
    total += residual * residual
  }
  (total / residuals.length().to_double()).sqrt()
}

///|
/// Computes gradients for a logistic loss at a supplied parameter vector.
pub fn logistic_gradient(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  coefficients : Array[Double],
  intercept : Double,
) -> (Array[Double], Double) {
  let n = if covariates.length() < treatment.length() {
    covariates.length()
  } else {
    treatment.length()
  }
  let gradient = Array::make(coefficients.length(), 0.0)
  let mut intercept_gradient = 0.0
  if n == 0 {
    return (gradient, intercept_gradient)
  }
  for i in 0..<n {
    let error = bounded_sigmoid(
        intercept + opt_dot(covariates[i], coefficients),
      ) -
      binary_target(treatment[i])
    intercept_gradient += error
    for j in 0..<coefficients.length() {
      gradient[j] += error * covariates[i][j]
    }
  }
  for j in 0..<gradient.length() {
    gradient[j] /= n.to_double()
  }
  (gradient, intercept_gradient / n.to_double())
}

///|
/// Computes a diagonal logistic Hessian approximation.
pub fn logistic_hessian_diagonal(
  covariates : Array[Array[Double]],
  coefficients : Array[Double],
  intercept : Double,
) -> Array[Double] {
  let result = Array::make(coefficients.length(), 0.0)
  let n = covariates.length()
  if n == 0 {
    return result
  }
  for row in covariates {
    let probability = bounded_sigmoid(intercept + opt_dot(row, coefficients))
    let variance = probability * (1.0 - probability)
    for j in 0..<coefficients.length() {
      result[j] += row[j] * row[j] * variance
    }
  }
  for j in 0..<result.length() {
    result[j] /= n.to_double()
  }
  result
}

///|
/// Fits a one-dimensional quadratic objective analytically.
pub fn solve_quadratic(
  objective : QuadraticObjective,
  initial? : Double = 0.0,
) -> Double {
  if objective.curvature <= 0.0 {
    initial
  } else {
    objective.center - objective.linear_term / objective.curvature
  }
}

///|
/// Evaluates a one-dimensional quadratic objective.
pub fn quadratic_value(
  objective : QuadraticObjective,
  value : Double,
) -> Double {
  0.5 *
  objective.curvature *
  (value - objective.center) *
  (value - objective.center) +
  objective.linear_term * value
}

///|
/// Performs Armijo backtracking for a scalar objective.
pub fn backtracking_line_search(
  value : Double,
  direction : Double,
  objective : (Double) -> Double,
  derivative : Double,
  initial_step? : Double = 1.0,
  shrink? : Double = 0.5,
  armijo? : Double = 1.0e-4,
  maximum_evaluations? : Int = 30,
) -> LineSearchResult {
  let baseline = objective(value)
  let mut step = initial_step
  let mut evaluations = 0
  let mut accepted = false
  let mut candidate_value = baseline
  while evaluations < maximum_evaluations {
    let candidate = value + step * direction
    candidate_value = objective(candidate)
    evaluations += 1
    if candidate_value <= baseline + armijo * step * derivative * direction {
      accepted = true
      break
    }
    step *= shrink
  }
  { step, objective: candidate_value, evaluations, accepted }
}

///|
/// Applies deterministic gradient descent to a scalar objective.
pub fn minimize_scalar(
  initial : Double,
  objective : (Double) -> Double,
  derivative : (Double) -> Double,
  config? : OptimizerConfig = default_optimizer_config(),
) -> (Double, OptimizationTrace) {
  let mut value = initial
  let losses : Array[Double] = Array::new(capacity=config.max_iterations)
  let mut gradient_norm = 0.0
  let mut step_norm = 0.0
  let mut converged = false
  let mut iterations = 0
  for iteration in 0..<config.max_iterations {
    let gradient = derivative(value)
    let step = -config.learning_rate * gradient
    value += step
    losses.push(objective(value))
    gradient_norm = gradient.abs()
    step_norm = step.abs()
    iterations = iteration + 1
    if gradient_norm < config.tolerance || step_norm < config.tolerance {
      converged = true
      break
    }
  }
  (value, append_trace(losses, iterations, converged, gradient_norm, step_norm))
}

///|
/// Returns a feature-wise variance score for screening.
pub fn feature_variances(covariates : Array[Array[Double]]) -> Array[Double] {
  if covariates.length() == 0 {
    return []
  }
  let width = covariates[0].length()
  let result = Array::make(width, 0.0)
  let means = Array::make(width, 0.0)
  for row in covariates {
    for j in 0..<width {
      means[j] += row[j]
    }
  }
  for j in 0..<width {
    means[j] /= covariates.length().to_double()
  }
  for row in covariates {
    for j in 0..<width {
      let d = row[j] - means[j]
      result[j] += d * d
    }
  }
  for j in 0..<width {
    result[j] /= covariates.length().to_double()
  }
  result
}

///|
/// Returns indexes of features passing a variance threshold.
pub fn select_variance_features(
  covariates : Array[Array[Double]],
  threshold? : Double = 1.0e-8,
) -> Array[Int] {
  let variances = feature_variances(covariates)
  let result : Array[Int] = Array::new()
  for j in 0..<variances.length() {
    if variances[j] > threshold {
      result.push(j)
    }
  }
  result
}

///|
/// Selects columns from a row-major matrix.
pub fn select_matrix_columns(
  matrix : Array[Array[Double]],
  indexes : Array[Int],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let selected = Array::new(capacity=indexes.length())
    for index in indexes {
      if index >= 0 && index < row.length() {
        selected.push(row[index])
      }
    }
    result.push(selected)
  }
  result
}

///|
/// Adds a deterministic intercept column to a matrix.
pub fn add_intercept_column(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=matrix.length())
  for row in matrix {
    let expanded = Array::new(capacity=row.length() + 1)
    expanded.push(1.0)
    for value in row {
      expanded.push(value)
    }
    result.push(expanded)
  }
  result
}

///|
/// Computes a finite-difference derivative for a scalar objective.
pub fn numerical_derivative(
  objective : (Double) -> Double,
  value : Double,
  epsilon? : Double = 1.0e-6,
) -> Double {
  (objective(value + epsilon) - objective(value - epsilon)) / (2.0 * epsilon)
}

///|
/// Computes a finite-difference second derivative.
pub fn numerical_second_derivative(
  objective : (Double) -> Double,
  value : Double,
  epsilon? : Double = 1.0e-4,
) -> Double {
  (
    objective(value + epsilon) -
    2.0 * objective(value) +
    objective(value - epsilon)
  ) /
  (epsilon * epsilon)
}

///|
/// Returns a stable binary threshold from a score and desired positive rate.
pub fn threshold_for_rate(scores : Array[Double], rate : Double) -> Double {
  if scores.length() == 0 {
    return 0.5
  }
  let bounded = clamp(rate, 0.0, 1.0)
  quantile(scores, 1.0 - bounded)
}

///|
/// Converts probabilities into treatment labels.
pub fn threshold_probabilities(
  probabilities : Array[Double],
  threshold : Double,
) -> Array[Bool] {
  let result = Array::new(capacity=probabilities.length())
  for probability in probabilities {
    result.push(probability >= threshold)
  }
  result
}

///|
/// Computes a confusion matrix in the order true-negative, false-positive, false-negative, true-positive.
pub fn confusion_counts(
  actual : Array[Bool],
  predicted : Array[Bool],
) -> Array[Int] {
  let n = if actual.length() < predicted.length() {
    actual.length()
  } else {
    predicted.length()
  }
  let result = Array::make(4, 0)
  for i in 0..<n {
    if !actual[i] && !predicted[i] {
      result[0] += 1
    } else if !actual[i] && predicted[i] {
      result[1] += 1
    } else if actual[i] && !predicted[i] {
      result[2] += 1
    } else {
      result[3] += 1
    }
  }
  result
}

///|
/// Computes precision, recall, specificity and F1 in a fixed order.
pub fn classification_summary(
  actual : Array[Bool],
  predicted : Array[Bool],
) -> Array[Double] {
  let counts = confusion_counts(actual, predicted)
  let tn = counts[0].to_double()
  let fp = counts[1].to_double()
  let false_negative = counts[2].to_double()
  let tp = counts[3].to_double()
  let precision = if tp + fp == 0.0 { 0.0 } else { tp / (tp + fp) }
  let recall = if tp + false_negative == 0.0 {
    0.0
  } else {
    tp / (tp + false_negative)
  }
  let specificity = if tn + fp == 0.0 { 0.0 } else { tn / (tn + fp) }
  let f1 = if precision + recall == 0.0 {
    0.0
  } else {
    2.0 * precision * recall / (precision + recall)
  }
  [precision, recall, specificity, f1]
}

///|
/// Returns a calibration table with equal-width score bins.
pub fn calibration_bins(
  probabilities : Array[Double],
  treatment : Array[Bool],
  bins? : Int = 10,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new(capacity=bins)
  let counts = Array::make(bins, 0)
  let predicted = Array::make(bins, 0.0)
  let observed = Array::make(bins, 0.0)
  let n = if probabilities.length() < treatment.length() {
    probabilities.length()
  } else {
    treatment.length()
  }
  for i in 0..<n {
    let index = (clamp(probabilities[i], 0.0, 0.999999) * bins.to_double()).to_int()
    counts[index] += 1
    predicted[index] += probabilities[i]
    if treatment[i] {
      observed[index] += 1.0
    }
  }
  for i in 0..<bins {
    let count = counts[i].to_double()
    result.push([
      if count == 0.0 {
        0.0
      } else {
        predicted[i] / count
      },
      if count == 0.0 {
        0.0
      } else {
        observed[i] / count
      },
      count,
    ])
  }
  result
}

///|
/// Computes expected calibration error from score bins.
pub fn expected_calibration_error(
  probabilities : Array[Double],
  treatment : Array[Bool],
  bins? : Int = 10,
) -> Double {
  let table = calibration_bins(probabilities, treatment, bins~)
  let n = probabilities.length().to_double()
  if n == 0.0 {
    return 0.0
  }
  let mut result = 0.0
  for row in table {
    result += row[2] / n * (row[0] - row[1]).abs()
  }
  result
}

///|
/// Computes a weighted ridge penalty excluding an intercept coefficient.
pub fn ridge_penalty(coefficients : Array[Double], strength : Double) -> Double {
  let mut result = 0.0
  for coefficient in coefficients {
    result += 0.5 * strength * coefficient * coefficient
  }
  result
}

///|
/// Computes a weighted lasso penalty.
pub fn lasso_penalty(coefficients : Array[Double], strength : Double) -> Double {
  let mut result = 0.0
  for coefficient in coefficients {
    result += strength * coefficient.abs()
  }
  result
}

///|
/// Computes an elastic-net penalty.
pub fn elastic_net_penalty(
  coefficients : Array[Double],
  l1 : Double,
  l2 : Double,
) -> Double {
  lasso_penalty(coefficients, l1) + ridge_penalty(coefficients, l2)
}

///|
/// Returns the number of effectively non-zero coefficients.
pub fn nonzero_coefficient_count(
  coefficients : Array[Double],
  tolerance? : Double = 1.0e-8,
) -> Int {
  let mut result = 0
  for coefficient in coefficients {
    if coefficient.abs() > tolerance {
      result += 1
    }
  }
  result
}

///|
/// Creates a compact optimization report vector.
pub fn optimization_summary(model : OptimizedLinearModel) -> Array[Double] {
  [
    model.trace.iterations.to_double(),
    if model.trace.converged {
      1.0
    } else {
      0.0
    },
    model.trace.initial_loss,
    model.trace.final_loss,
    model.trace.gradient_norm,
    model.trace.step_norm,
    model.coefficients.length().to_double(),
    nonzero_coefficient_count(model.coefficients).to_double(),
  ]
}
