///|
pub fn bootstrap_ipw_ate(
  outcomes : Array[Double],
  treatment : Array[Bool],
  propensity_scores : Array[Double],
  replicates : Int,
  seed : UInt64,
) -> BootstrapSummary {
  let point = estimate_ipw_ate(outcomes, treatment, propensity_scores).estimate
  let rng = RandomState::new(seed)
  let samples : Array[Double] = Array::new(capacity=replicates)
  let mut successful = 0
  for _ in 0..<replicates {
    let indices = resample_indices(rng, outcomes.length())
    let y = resampled_values(outcomes, indices)
    let t = resampled_bool(treatment, indices)
    let p = resampled_values(propensity_scores, indices)
    let mut treated = 0
    let mut control = 0
    for value in t {
      if value {
        treated += 1
      } else {
        control += 1
      }
    }
    if treated > 0 && control > 0 {
      samples.push(estimate_ipw_ate(y, t, p).estimate)
      successful += 1
    }
  }
  {
    point_estimate: point,
    standard_error: std_dev(samples),
    lower: quantile(samples, 0.025),
    upper: quantile(samples, 0.975),
    replicates,
    successful_replicates: successful,
  }
}

///|
pub fn bootstrap_percentile_interval(
  samples : Array[Double],
  level : Double,
) -> (Double, Double) {
  let alpha = (1.0 - clamp(level, 0.0, 1.0)) / 2.0
  (quantile(samples, alpha), quantile(samples, 1.0 - alpha))
}

///|
pub fn bootstrap_bias(
  samples : Array[Double],
  point_estimate : Double,
) -> Double {
  mean(samples) - point_estimate
}

///|
pub fn bootstrap_relative_bias(
  samples : Array[Double],
  point_estimate : Double,
) -> Double {
  if point_estimate == 0.0 {
    0.0
  } else {
    bootstrap_bias(samples, point_estimate) / point_estimate
  }
}
