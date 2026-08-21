///|
/// Robust regression fit with Huber loss.
pub struct RobustRegressionFit {
  coefficients : Array[Double]
  intercept : Double
  scale : Double
  iterations : Int
  converged : Bool
  objective : Double
  residuals : Array[Double]
}

///|
/// Residual influence diagnostics.
pub struct ResidualDiagnostics {
  residuals : Array[Double]
  standardized : Array[Double]
  leverage : Array[Double]
  cooks_distance : Array[Double]
  outlier_count : Int
  maximum_leverage : Double
}

///|
/// Robust location and scale estimate.
pub struct RobustLocationScale {
  location : Double
  scale : Double
  q25 : Double
  median : Double
  q75 : Double
  iterations : Int
}

///|
fn rr_dot(first : Array[Double], second : Array[Double]) -> Double {
  let n = first.length().min(second.length())
  let mut result = 0.0
  for i in 0..<n {
    result += first[i] * second[i]
  }
  result
}

///|
/// Computes a robust median/MAD location-scale summary.
pub fn robust_location_scale(values : Array[Double]) -> RobustLocationScale {
  let observed = values.filter(fn(value) { is_finite(value) })
  let center = median(observed)
  let deviations = Array::new(capacity=observed.length())
  for value in observed {
    deviations.push((value - center).abs())
  }
  {
    location: center,
    scale: (1.4826 * median(deviations)).max(1.0e-12),
    q25: quantile(observed, 0.25),
    median: center,
    q75: quantile(observed, 0.75),
    iterations: 1,
  }
}

///|
/// Huber loss at a residual and tuning constant.
pub fn huber_loss(residual : Double, tuning? : Double = 1.345) -> Double {
  let magnitude = residual.abs()
  if magnitude <= tuning {
    0.5 * residual * residual
  } else {
    tuning * (magnitude - 0.5 * tuning)
  }
}

///|
/// Huber score weight.
pub fn huber_score_weight(
  residual : Double,
  tuning? : Double = 1.345,
) -> Double {
  let magnitude = residual.abs()
  if magnitude == 0.0 || magnitude <= tuning {
    1.0
  } else {
    tuning / magnitude
  }
}

///|
/// Fits a robust linear model with iteratively reweighted least squares.
pub fn fit_robust_regression(
  covariates : Array[Array[Double]],
  outcomes : Array[Double],
  tuning? : Double = 1.345,
  iterations? : Int = 40,
  tolerance? : Double = 1.0e-7,
) -> RobustRegressionFit {
  let n = covariates.length().min(outcomes.length())
  if n == 0 {
    return {
      coefficients: [],
      intercept: 0.0,
      scale: 0.0,
      iterations: 0,
      converged: false,
      objective: 0.0,
      residuals: [],
    }
  }
  let width = covariates[0].length()
  let coefficients = Array::make(width, 0.0)
  let mut intercept = mean(outcomes[:n].to_owned())
  let mut scale = std_dev(outcomes[:n].to_owned()).max(1.0e-12)
  let mut converged = false
  let mut used_iterations = 0
  for iteration in 0..<iterations {
    let residuals = Array::new(capacity=n)
    let weights = Array::new(capacity=n)
    for i in 0..<n {
      let residual = outcomes[i] -
        intercept -
        rr_dot(covariates[i], coefficients)
      residuals.push(residual)
      weights.push(huber_score_weight(residual / scale, tuning~))
    }
    let mut weighted_residual = 0.0
    let mut weight_total = 0.0
    for i in 0..<n {
      weighted_residual += weights[i] * residuals[i]
      weight_total += weights[i]
    }
    let intercept_step = if weight_total == 0.0 {
      0.0
    } else {
      weighted_residual / weight_total
    }
    intercept += intercept_step
    let mut max_step = intercept_step.abs()
    for j in 0..<width {
      let mut numerator = 0.0
      let mut denominator = 0.0
      for i in 0..<n {
        numerator += weights[i] *
          covariates[i][j] *
          (outcomes[i] - intercept - rr_dot(covariates[i], coefficients))
        denominator += weights[i] * covariates[i][j] * covariates[i][j]
      }
      let step = if denominator == 0.0 { 0.0 } else { numerator / denominator }
      coefficients[j] += step
      if step.abs() > max_step {
        max_step = step.abs()
      }
    }
    scale = median_absolute_deviation(residuals).max(1.0e-12)
    used_iterations = iteration + 1
    if max_step < tolerance {
      converged = true
      break
    }
  }
  let final_residuals = Array::new(capacity=n)
  let mut objective = 0.0
  for i in 0..<n {
    let residual = outcomes[i] - intercept - rr_dot(covariates[i], coefficients)
    final_residuals.push(residual)
    objective += huber_loss(residual / scale, tuning~)
  }
  {
    coefficients,
    intercept,
    scale,
    iterations: used_iterations,
    converged,
    objective,
    residuals: final_residuals,
  }
}

///|
/// Predicts outcomes from a robust regression fit.
pub fn robust_regression_predict(
  model : RobustRegressionFit,
  covariates : Array[Array[Double]],
) -> Array[Double] {
  let result = Array::new(capacity=covariates.length())
  for row in covariates {
    result.push(model.intercept + rr_dot(row, model.coefficients))
  }
  result
}

///|
/// Computes robust regression residual diagnostics.
pub fn robust_residual_diagnostics(
  model : RobustRegressionFit,
  covariates : Array[Array[Double]],
) -> ResidualDiagnostics {
  let residuals = model.residuals.copy()
  let scale = model.scale.max(1.0e-12)
  let standardized = Array::new(capacity=residuals.length())
  for residual in residuals {
    standardized.push(residual / scale)
  }
  let leverage = Array::new(capacity=covariates.length())
  for row in covariates {
    let norm = rr_dot(row, row)
    leverage.push(norm / (1.0 + norm))
  }
  let cooks = Array::new(capacity=residuals.length())
  let mut outliers = 0
  let mut maximum = 0.0
  for i in 0..<residuals.length() {
    let leverage_value = if i < leverage.length() { leverage[i] } else { 0.0 }
    let value = standardized[i] * standardized[i] * leverage_value
    cooks.push(value)
    if standardized[i].abs() > 3.0 {
      outliers += 1
    }
    if leverage_value > maximum {
      maximum = leverage_value
    }
  }
  {
    residuals,
    standardized,
    leverage,
    cooks_distance: cooks,
    outlier_count: outliers,
    maximum_leverage: maximum,
  }
}

///|
/// Computes a Huber-weighted mean.
pub fn huber_mean(
  values : Array[Double],
  tuning? : Double = 1.345,
  iterations? : Int = 20,
) -> Double {
  let mut center = median(values)
  let scale = median_absolute_deviation(values).max(1.0e-12)
  for _ in 0..<iterations {
    let mut numerator = 0.0
    let mut denominator = 0.0
    for value in values {
      let weight = huber_score_weight((value - center) / scale, tuning~)
      numerator += weight * value
      denominator += weight
    }
    if denominator > 0.0 {
      center = numerator / denominator
    }
  }
  center
}

///|
/// Computes a robust weighted covariance matrix.
pub fn robust_covariance_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  if matrix.length() == 0 {
    return []
  }
  let width = matrix[0].length()
  let centers = Array::make(width, 0.0)
  for j in 0..<width {
    centers[j] = huber_mean(matrix.map(fn(row) { row[j] }))
  }
  let result = Array::make(width, Array::make(width, 0.0))
  for row in matrix {
    for j in 0..<width {
      for k in 0..<width {
        result[j][k] += (row[j] - centers[j]) * (row[k] - centers[k])
      }
    }
  }
  for j in 0..<width {
    for k in 0..<width {
      result[j][k] /= matrix.length().to_double().max(1.0)
    }
  }
  result
}

///|
/// Computes an outlier-resistant standard error of the mean.
pub fn robust_mean_standard_error(values : Array[Double]) -> Double {
  let center = huber_mean(values)
  let residuals = Array::new(capacity=values.length())
  for value in values {
    residuals.push(value - center)
  }
  let scale = median_absolute_deviation(values)
  if values.length() == 0 {
    0.0
  } else {
    scale / values.length().to_double().sqrt()
  }
}

///|
/// Computes a robust effect size between two groups.
pub fn robust_group_effect(
  first : Array[Double],
  second : Array[Double],
) -> AdvancedEffect {
  let estimate = huber_mean(first) - huber_mean(second)
  let standard_error = (robust_mean_standard_error(first) *
  robust_mean_standard_error(first) +
  robust_mean_standard_error(second) * robust_mean_standard_error(second)).sqrt()
  let interval = normal_confidence_interval(estimate, standard_error)
  {
    estimate,
    standard_error,
    lower: interval[0],
    upper: interval[1],
    effective_sample_size: first.length().to_double() +
    second.length().to_double(),
    estimand: "robust group effect",
    passes: first.length() > 2 && second.length() > 2,
  }
}

///|
/// Returns a compact robust regression summary vector.
pub fn robust_regression_summary(model : RobustRegressionFit) -> Array[Double] {
  [
    model.intercept,
    model.scale,
    model.iterations.to_double(),
    if model.converged {
      1.0
    } else {
      0.0
    },
    model.objective,
    model.coefficients.length().to_double(),
  ]
}
