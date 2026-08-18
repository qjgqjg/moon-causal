///|
pub struct OverlapReport {
  minimum : Double
  maximum : Double
  treated_minimum : Double
  treated_maximum : Double
  control_minimum : Double
  control_maximum : Double
  common_support_lower : Double
  common_support_upper : Double
  overlap_fraction : Double
}

///|
pub fn overlap_report(
  propensity_scores : Array[Double],
  treatment : Array[Bool],
) -> OverlapReport {
  let treated = group_values(propensity_scores, treatment, true)
  let control = group_values(propensity_scores, treatment, false)
  let minimum = quantile(propensity_scores, 0.0)
  let maximum = quantile(propensity_scores, 1.0)
  let treated_minimum = quantile(treated, 0.0)
  let treated_maximum = quantile(treated, 1.0)
  let control_minimum = quantile(control, 0.0)
  let control_maximum = quantile(control, 1.0)
  let lower = if treated_minimum > control_minimum {
    treated_minimum
  } else {
    control_minimum
  }
  let upper = if treated_maximum < control_maximum {
    treated_maximum
  } else {
    control_maximum
  }
  let mut inside = 0
  for score in propensity_scores {
    if score >= lower && score <= upper {
      inside += 1
    }
  }
  {
    minimum,
    maximum,
    treated_minimum,
    treated_maximum,
    control_minimum,
    control_maximum,
    common_support_lower: lower,
    common_support_upper: upper,
    overlap_fraction: if propensity_scores.length() == 0 {
      0.0
    } else {
      inside.to_double() / propensity_scores.length().to_double()
    },
  }
}

///|
pub fn positivity_violations(
  propensity_scores : Array[Double],
  epsilon : Double,
) -> Int {
  let threshold = clamp(epsilon, 0.0, 0.5)
  let mut violations = 0
  for score in propensity_scores {
    if score <= threshold || score >= 1.0 - threshold {
      violations += 1
    }
  }
  violations
}

///|
pub fn propensity_calibration(
  propensity_scores : Array[Double],
  treatment : Array[Bool],
  bins : Int,
) -> Array[Double] {
  let result : Array[Double] = Array::new(
    capacity=if bins > 0 { bins } else { 0 },
  )
  if bins <= 0 {
    return result
  }
  for bin in 0..<bins {
    let lower = bin.to_double() / bins.to_double()
    let upper = (bin + 1).to_double() / bins.to_double()
    let mut expected = 0.0
    let mut observed = 0.0
    let mut count = 0
    for i in 0..<propensity_scores.length() {
      if propensity_scores[i] >= lower &&
        (propensity_scores[i] < upper || bin == bins - 1) {
        expected += propensity_scores[i]
        if treatment[i] {
          observed += 1.0
        }
        count += 1
      }
    }
    if count == 0 {
      result.push(0.0)
    } else {
      result.push((observed - expected) / count.to_double())
    }
  }
  result
}

///|
pub fn standardized_residuals(
  outcomes : Array[Double],
  predictions : Array[Double],
) -> Array[Double] {
  let residual = Array::new(capacity=outcomes.length())
  let n = if outcomes.length() < predictions.length() {
    outcomes.length()
  } else {
    predictions.length()
  }
  for i in 0..<n {
    residual.push(outcomes[i] - predictions[i])
  }
  let scale = std_dev(residual, sample=false)
  if scale == 0.0 {
    return Array::make(residual.length(), 0.0)
  }
  for i in 0..<residual.length() {
    residual[i] /= scale
  }
  residual
}

///|
pub fn mean_absolute_error(
  outcomes : Array[Double],
  predictions : Array[Double],
) -> Double {
  let residuals = standardized_residuals(outcomes, predictions)
  let mut result = 0.0
  for value in residuals {
    result += value.abs()
  }
  if residuals.length() == 0 {
    0.0
  } else {
    result / residuals.length().to_double()
  }
}

///|
pub fn root_mean_squared_error(
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
