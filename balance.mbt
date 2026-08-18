///|
fn variance_ratio(treated : Array[Double], control : Array[Double]) -> Double {
  let denominator = variance(control, sample=true)
  if denominator == 0.0 {
    1.0
  } else {
    variance(treated, sample=true) / denominator
  }
}

///|
fn balance_metric_for_column(
  covariate : Array[Double],
  treatment : Array[Bool],
  name : String,
  threshold : Double,
) -> BalanceMetric {
  let treated = group_values(covariate, treatment, true)
  let control = group_values(covariate, treatment, false)
  let treated_mean = mean(treated)
  let control_mean = mean(control)
  let deviation = smd(covariate, treatment)
  let absolute_difference = (treated_mean - control_mean).abs()
  {
    name,
    treated_mean,
    control_mean,
    standardized_difference: deviation,
    variance_ratio: variance_ratio(treated, control),
    absolute_difference,
    balanced: deviation.abs() <= threshold,
  }
}

///|
/// Computes covariate balance diagnostics before or after an adjustment.
pub fn balance_table(
  covariates : Array[Array[Double]],
  treatment : Array[Bool],
  names : Array[String],
  threshold? : Double = 0.1,
) -> Array[BalanceMetric] {
  if covariates.length() == 0 {
    return []
  }
  let width = covariates[0].length()
  let result : Array[BalanceMetric] = Array::new(capacity=width)
  for j in 0..<width {
    let column = Array::new(capacity=covariates.length())
    for row in covariates {
      column.push(row[j])
    }
    let name = if j < names.length() { names[j] } else { "x{j + 1}" }
    result.push(balance_metric_for_column(column, treatment, name, threshold))
  }
  result
}

///|
/// Returns the largest absolute standardized difference in a balance table.
pub fn maximum_absolute_smd(metrics : Array[BalanceMetric]) -> Double {
  let mut result = 0.0
  for metric in metrics {
    if metric.standardized_difference.abs() > result {
      result = metric.standardized_difference.abs()
    }
  }
  result
}

///|
pub fn balanced_covariate_count(metrics : Array[BalanceMetric]) -> Int {
  let mut result = 0
  for metric in metrics {
    if metric.balanced {
      result += 1
    }
  }
  result
}

///|
pub fn weighted_smd(
  covariate : Array[Double],
  treatment : Array[Bool],
  weights : Array[Double],
) -> Double {
  let treated = Array::new()
  let controls = Array::new()
  let treated_weights = Array::new()
  let control_weights = Array::new()
  for i in 0..<covariate.length() {
    if treatment[i] {
      treated.push(covariate[i])
      treated_weights.push(weights[i])
    } else {
      controls.push(covariate[i])
      control_weights.push(weights[i])
    }
  }
  let difference = weighted_mean(treated, treated_weights) -
    weighted_mean(controls, control_weights)
  let pooled = (
      weighted_variance(treated, treated_weights) +
      weighted_variance(controls, control_weights)
    ) /
    2.0
  if pooled <= 0.0 {
    0.0
  } else {
    difference / pooled.sqrt()
  }
}
