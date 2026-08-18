///|
pub fn covariance(a : Array[Double], b : Array[Double]) -> Double {
  if a.length() == 0 || a.length() != b.length() {
    return 0.0
  }
  let mut result = 0.0
  let mean_a = mean(a)
  let mean_b = mean(b)
  for i in 0..<a.length() {
    result += (a[i] - mean_a) * (b[i] - mean_b)
  }
  let denominator = if a.length() <= 1 {
    1.0
  } else {
    (a.length() - 1).to_double()
  }
  result / denominator
}

///|
pub fn correlation(a : Array[Double], b : Array[Double]) -> Double {
  let denominator = std_dev(a) * std_dev(b)
  if denominator == 0.0 {
    0.0
  } else {
    covariance(a, b) / denominator
  }
}

///|
pub fn r_squared(
  outcomes : Array[Double],
  predictions : Array[Double],
) -> Double {
  let n = if outcomes.length() < predictions.length() {
    outcomes.length()
  } else {
    predictions.length()
  }
  if n == 0 {
    return 0.0
  }
  let center = mean(outcomes)
  let mut total = 0.0
  let mut residual = 0.0
  for i in 0..<n {
    let total_difference = outcomes[i] - center
    let residual_difference = outcomes[i] - predictions[i]
    total += total_difference * total_difference
    residual += residual_difference * residual_difference
  }
  if total == 0.0 {
    0.0
  } else {
    1.0 - residual / total
  }
}

///|
pub fn log_loss(
  probabilities : Array[Double],
  treatment : Array[Bool],
) -> Double {
  let n = if probabilities.length() < treatment.length() {
    probabilities.length()
  } else {
    treatment.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut result = 0.0
  for i in 0..<n {
    let probability = safe_probability(probabilities[i])
    let target = if treatment[i] { 1.0 } else { 0.0 }
    result -= target * @math.ln(probability) +
      (1.0 - target) * @math.ln(1.0 - probability)
  }
  result / n.to_double()
}

///|
pub fn brier_score(
  probabilities : Array[Double],
  treatment : Array[Bool],
) -> Double {
  let n = if probabilities.length() < treatment.length() {
    probabilities.length()
  } else {
    treatment.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut result = 0.0
  for i in 0..<n {
    let target = if treatment[i] { 1.0 } else { 0.0 }
    let difference = probabilities[i] - target
    result += difference * difference
  }
  result / n.to_double()
}

///|
pub fn root_mean_squared_error_raw(
  outcomes : Array[Double],
  predictions : Array[Double],
) -> Double {
  let n = if outcomes.length() < predictions.length() {
    outcomes.length()
  } else {
    predictions.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut result = 0.0
  for i in 0..<n {
    let difference = outcomes[i] - predictions[i]
    result += difference * difference
  }
  (result / n.to_double()).sqrt()
}

///|
pub fn mean_absolute_error_raw(
  outcomes : Array[Double],
  predictions : Array[Double],
) -> Double {
  let n = if outcomes.length() < predictions.length() {
    outcomes.length()
  } else {
    predictions.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut result = 0.0
  for i in 0..<n {
    result += (outcomes[i] - predictions[i]).abs()
  }
  result / n.to_double()
}

///|
pub fn treatment_t_statistic(
  estimate : Double,
  standard_error : Double,
) -> Double {
  if standard_error == 0.0 {
    0.0
  } else {
    estimate / standard_error
  }
}

///|
pub fn calibration_error(
  probabilities : Array[Double],
  treatment : Array[Bool],
  bins : Int,
) -> Double {
  let calibration = propensity_calibration(probabilities, treatment, bins)
  let mut result = 0.0
  for value in calibration {
    result += value.abs()
  }
  if calibration.length() == 0 {
    0.0
  } else {
    result / calibration.length().to_double()
  }
}

///|
pub fn standardized_effect_size(
  estimate : Double,
  outcome : Array[Double],
) -> Double {
  let scale = std_dev(outcome, sample=false)
  if scale == 0.0 {
    0.0
  } else {
    estimate / scale
  }
}
