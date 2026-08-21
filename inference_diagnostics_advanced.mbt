///|
/// A confidence interval with the critical value used to construct it.
pub struct InferenceInterval {
  estimate : Double
  standard_error : Double
  critical_value : Double
  lower : Double
  upper : Double
  width : Double
  confidence_level : Double
  passes : Bool
}

///|
/// One named contrast with a multiplicity-adjusted decision.
pub struct InferenceContrast {
  name : String
  estimate : Double
  standard_error : Double
  z_score : Double
  p_value : Double
  interval : InferenceInterval
  significant : Bool
}

///|
/// Collection-level inference diagnostics.
pub struct InferenceSummary {
  contrasts : Int
  finite_contrasts : Int
  significant_contrasts : Int
  mean_standard_error : Double
  mean_interval_width : Double
  minimum_p_value : Double
  passes : Bool
}

///|
/// A power calculation at one sample size and effect size.
pub struct InferencePowerPoint {
  sample_size : Int
  effect_size : Double
  standard_error : Double
  power : Double
  alpha : Double
  passes : Bool
}

///|
/// A finite-safe normal cumulative distribution approximation.
pub fn inference_normal_cdf(value : Double) -> Double {
  if !is_finite(value) {
    if value < 0.0 {
      0.0
    } else {
      1.0
    }
  } else {
    let sign = if value < 0.0 { -1.0 } else { 1.0 }
    let x = value.abs() / 1.4142135623730951
    let t = 1.0 / (1.0 + 0.3275911 * x)
    let polynomial = 1.061405429 * t -
      1.453152027 * t * t +
      1.421413741 * t * t * t -
      0.284496736 * t * t * t * t +
      0.254829592 * t * t * t * t * t
    let erf = 1.0 - polynomial * @math.exp(-x * x)
    0.5 * (1.0 + sign * erf)
  }
}

///|
/// Computes a normal quantile by bounded binary search.
pub fn inference_normal_quantile(probability : Double) -> Double {
  let target = clamp(probability, 1.0e-9, 1.0 - 1.0e-9)
  let mut lower = -9.0
  let mut upper = 9.0
  for _ in 0..<80 {
    let middle = (lower + upper) / 2.0
    if inference_normal_cdf(middle) < target {
      lower = middle
    } else {
      upper = middle
    }
  }
  (lower + upper) / 2.0
}

///|
/// Returns a critical value for a two-sided normal interval.
pub fn inference_critical_value(confidence_level? : Double = 0.95) -> Double {
  let confidence = clamp(confidence_level, 1.0e-6, 1.0 - 1.0e-6)
  inference_normal_quantile((1.0 + confidence) / 2.0)
}

///|
/// Constructs a normal-approximation interval.
pub fn inference_interval(
  estimate : Double,
  standard_error : Double,
  confidence_level? : Double = 0.95,
) -> InferenceInterval {
  let confidence = clamp(confidence_level, 1.0e-6, 1.0 - 1.0e-6)
  let error = if is_finite(standard_error) { standard_error.abs() } else { 0.0 }
  let critical = inference_critical_value(confidence_level=confidence)
  let margin = critical * error
  {
    estimate,
    standard_error: error,
    critical_value: critical,
    lower: estimate - margin,
    upper: estimate + margin,
    width: 2.0 * margin,
    confidence_level: confidence,
    passes: is_finite(estimate) && is_finite(error),
  }
}

///|
/// Computes a standardized contrast score.
pub fn inference_z_score(estimate : Double, standard_error : Double) -> Double {
  if !is_finite(estimate) ||
    !is_finite(standard_error) ||
    standard_error.abs() <= 1.0e-12 {
    0.0
  } else {
    estimate / standard_error.abs()
  }
}

///|
/// Computes a two-sided normal p-value.
pub fn inference_two_sided_p(z_score : Double) -> Double {
  let value = inference_normal_cdf(-z_score.abs()) * 2.0
  clamp(value, 0.0, 1.0)
}

///|
/// Creates a named contrast and applies an alpha threshold.
pub fn inference_contrast(
  name : String,
  estimate : Double,
  standard_error : Double,
  confidence_level? : Double = 0.95,
  alpha? : Double = 0.05,
) -> InferenceContrast {
  let interval = inference_interval(estimate, standard_error, confidence_level~)
  let z_score = inference_z_score(estimate, standard_error)
  let p_value = inference_two_sided_p(z_score)
  {
    name,
    estimate,
    standard_error: standard_error.abs(),
    z_score,
    p_value,
    interval,
    significant: p_value <= clamp(alpha, 1.0e-9, 1.0),
  }
}

///|
/// Creates contrasts with a Bonferroni-adjusted alpha level.
pub fn inference_contrasts(
  names : Array[String],
  estimates : Array[Double],
  standard_errors : Array[Double],
  confidence_level? : Double = 0.95,
  alpha? : Double = 0.05,
) -> Array[InferenceContrast] {
  let n = names.length().min(estimates.length()).min(standard_errors.length())
  let adjusted = clamp(alpha, 1.0e-9, 1.0) / n.to_double().max(1.0)
  let result : Array[InferenceContrast] = Array::new(capacity=n)
  for i in 0..<n {
    result.push(
      inference_contrast(
        names[i],
        estimates[i],
        standard_errors[i],
        confidence_level~,
        alpha=adjusted,
      ),
    )
  }
  result
}

///|
/// Summarizes a set of named contrasts.
pub fn inference_summary(
  contrasts : Array[InferenceContrast],
) -> InferenceSummary {
  let mut finite = 0
  let mut significant = 0
  let mut error_sum = 0.0
  let mut width_sum = 0.0
  let mut minimum_p = 1.0
  for contrast in contrasts {
    if is_finite(contrast.estimate) && is_finite(contrast.standard_error) {
      finite += 1
      error_sum += contrast.standard_error
      width_sum += contrast.interval.width
      if contrast.p_value < minimum_p {
        minimum_p = contrast.p_value
      }
    }
    if contrast.significant {
      significant += 1
    }
  }
  {
    contrasts: contrasts.length(),
    finite_contrasts: finite,
    significant_contrasts: significant,
    mean_standard_error: if finite == 0 {
      0.0
    } else {
      error_sum / finite.to_double()
    },
    mean_interval_width: if finite == 0 {
      0.0
    } else {
      width_sum / finite.to_double()
    },
    minimum_p_value: minimum_p,
    passes: contrasts.length() > 0 && finite == contrasts.length(),
  }
}

///|
/// Calculates empirical coverage for intervals against one truth value.
pub fn inference_coverage(
  truth : Double,
  intervals : Array[InferenceInterval],
) -> Double {
  if intervals.length() == 0 {
    0.0
  } else {
    let covered = intervals.fold(init=0, fn(total, interval) {
      if interval.lower <= truth && truth <= interval.upper {
        total + 1
      } else {
        total
      }
    })
    covered.to_double() / intervals.length().to_double()
  }
}

///|
/// Computes bias relative to a known target.
pub fn inference_bias(truth : Double, estimates : Array[Double]) -> Double {
  mean_or(estimates, truth) - truth
}

///|
/// Computes root mean squared error relative to a known target.
pub fn inference_rmse(truth : Double, estimates : Array[Double]) -> Double {
  if estimates.length() == 0 {
    0.0
  } else {
    let mut squared = 0.0
    for estimate in estimates {
      let error = estimate - truth
      squared += error * error
    }
    (squared / estimates.length().to_double()).sqrt()
  }
}

///|
/// Computes Monte Carlo standard error for repeated estimates.
pub fn inference_monte_carlo_se(estimates : Array[Double]) -> Double {
  if estimates.length() == 0 {
    0.0
  } else {
    std_dev(estimates) / estimates.length().to_double().sqrt()
  }
}

///|
/// Computes a finite-sample standard error for a mean.
pub fn inference_mean_standard_error(
  values : Array[Double],
  weights? : Array[Double] = [],
) -> Double {
  if weights.length() == 0 {
    if values.length() == 0 {
      0.0
    } else {
      std_dev(values) / values.length().to_double().sqrt()
    }
  } else {
    weighted_variance(values, weights).sqrt() /
    effective_sample_size(weights).sqrt().max(1.0)
  }
}

///|
/// Computes a cluster-robust design effect.
pub fn inference_design_effect(
  cluster_sizes : Array[Int],
  intraclass_correlation : Double,
) -> Double {
  if cluster_sizes.length() == 0 {
    1.0
  } else {
    let mean_size = cluster_sizes
      .fold(init=0, fn(total, value) { total + value })
      .to_double() /
      cluster_sizes.length().to_double()
    1.0 + (mean_size - 1.0).max(0.0) * clamp(intraclass_correlation, 0.0, 1.0)
  }
}

///|
/// Returns the standard error inflated by a design effect.
pub fn inference_cluster_standard_error(
  standard_error : Double,
  cluster_sizes : Array[Int],
  intraclass_correlation : Double,
) -> Double {
  standard_error.abs() *
  inference_design_effect(cluster_sizes, intraclass_correlation).sqrt()
}

///|
/// Estimates the required sample per arm for a normal mean contrast.
pub fn inference_required_sample_size(
  standard_deviation : Double,
  minimum_effect : Double,
  alpha? : Double = 0.05,
  power? : Double = 0.8,
  design_effect? : Double = 1.0,
) -> Int {
  let effect = minimum_effect.abs().max(1.0e-12)
  let deviation = standard_deviation.abs()
  let critical = inference_normal_quantile(
    1.0 - clamp(alpha, 1.0e-9, 0.5) / 2.0,
  )
  let power_critical = inference_normal_quantile(clamp(power, 0.5, 0.999999))
  let numerator = (critical + power_critical) * deviation
  let ratio = numerator / effect
  let raw = 2.0 * ratio * ratio * design_effect.max(1.0)
  raw.ceil().to_int().max(2)
}

///|
/// Computes normal-approximation power for a two-arm contrast.
pub fn inference_power(
  sample_size_per_arm : Int,
  effect : Double,
  standard_deviation : Double,
  alpha? : Double = 0.05,
  design_effect? : Double = 1.0,
) -> Double {
  let n = sample_size_per_arm.max(1).to_double()
  let standard_error = standard_deviation.abs() *
    design_effect.max(1.0).sqrt() *
    (2.0 / n).sqrt()
  let critical = inference_normal_quantile(
    1.0 - clamp(alpha, 1.0e-9, 0.5) / 2.0,
  )
  let noncentral = effect.abs() / standard_error.max(1.0e-12)
  clamp(
    inference_normal_cdf(-critical - noncentral) +
    1.0 -
    inference_normal_cdf(critical - noncentral),
    0.0,
    1.0,
  )
}

///|
/// Creates one power curve point.
pub fn inference_power_point(
  sample_size_per_arm : Int,
  effect : Double,
  standard_deviation : Double,
  alpha? : Double = 0.05,
  design_effect? : Double = 1.0,
) -> InferencePowerPoint {
  let n = sample_size_per_arm.max(1)
  let power = inference_power(
    n,
    effect,
    standard_deviation,
    alpha~,
    design_effect~,
  )
  let standard_error = standard_deviation.abs() * (2.0 / n.to_double()).sqrt()
  {
    sample_size: n,
    effect_size: effect,
    standard_error,
    power,
    alpha: clamp(alpha, 1.0e-9, 1.0),
    passes: power >= 0.8,
  }
}

///|
/// Computes an evenly spaced power curve.
pub fn inference_power_curve(
  sample_sizes : Array[Int],
  effect : Double,
  standard_deviation : Double,
  alpha? : Double = 0.05,
) -> Array[InferencePowerPoint] {
  let result : Array[InferencePowerPoint] = Array::new(
    capacity=sample_sizes.length(),
  )
  for sample_size in sample_sizes {
    result.push(
      inference_power_point(sample_size, effect, standard_deviation, alpha~),
    )
  }
  result
}

///|
/// Returns the first sample size on a power curve meeting a target.
pub fn inference_first_power_sample(
  curve : Array[InferencePowerPoint],
  target_power? : Double = 0.8,
) -> Int {
  let target = clamp(target_power, 0.0, 1.0)
  for point in curve {
    if point.power >= target {
      return point.sample_size
    }
  }
  0
}

///|
/// Converts a bootstrap summary into a reusable interval.
pub fn inference_bootstrap_interval(
  summary : BootstrapSummary,
  confidence_level? : Double = 0.95,
) -> InferenceInterval {
  let standard_error = summary.standard_error.abs()
  let interval = inference_interval(
    summary.point_estimate,
    standard_error,
    confidence_level~,
  )
  {
    estimate: summary.point_estimate,
    standard_error,
    critical_value: interval.critical_value,
    lower: summary.lower,
    upper: summary.upper,
    width: (summary.upper - summary.lower).abs(),
    confidence_level: interval.confidence_level,
    passes: summary.successful_replicates > 0 &&
    summary.successful_replicates <= summary.replicates,
  }
}

///|
/// Computes a leave-one-out sensitivity curve for a mean.
pub fn inference_leave_one_out(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new(capacity=values.length())
  let total = values.fold(init=0.0, fn(sum, value) { sum + value })
  for i in 0..<values.length() {
    let remaining = values.length() - 1
    result.push(
      if remaining == 0 {
        0.0
      } else {
        (total - values[i]) / remaining.to_double()
      },
    )
  }
  result
}

///|
/// Computes the maximum absolute leave-one-out change.
pub fn inference_leave_one_out_sensitivity(values : Array[Double]) -> Double {
  let baseline = mean_or(values, 0.0)
  let curve = inference_leave_one_out(values)
  let mut maximum = 0.0
  for value in curve {
    let difference = (value - baseline).abs()
    if difference > maximum {
      maximum = difference
    }
  }
  maximum
}

///|
/// Computes a robust sandwich-style standard error from influence values.
pub fn inference_influence_standard_error(
  influence : Array[Double],
  sample_size : Int,
) -> Double {
  let n = sample_size.max(1).to_double()
  (variance(influence) / n).sqrt()
}

///|
/// Computes a finite-safe effect precision score.
pub fn inference_precision_score(
  estimate : Double,
  standard_error : Double,
) -> Double {
  let denominator = estimate.abs().max(1.0e-12)
  if !is_finite(estimate) || !is_finite(standard_error) {
    0.0
  } else {
    clamp(estimate.abs() / denominator / (1.0 + standard_error.abs()), 0.0, 1.0)
  }
}

///|
/// Returns a stable summary vector for one interval.
pub fn inference_interval_summary(
  interval : InferenceInterval,
) -> Array[Double] {
  [
    interval.estimate,
    interval.standard_error,
    interval.critical_value,
    interval.lower,
    interval.upper,
    interval.width,
    interval.confidence_level,
    if interval.passes {
      1.0
    } else {
      0.0
    },
  ]
}

///|
/// Returns a stable fingerprint for a set of contrasts.
pub fn inference_contrast_fingerprint(
  contrasts : Array[InferenceContrast],
) -> UInt64 {
  let rows : Array[Array[Double]] = Array::new(capacity=contrasts.length())
  for contrast in contrasts {
    rows.push([
      contrast.estimate,
      contrast.standard_error,
      contrast.z_score,
      contrast.p_value,
      contrast.interval.lower,
      contrast.interval.upper,
    ])
  }
  matrix_checksum(rows)
}

///|
/// Serializes inference results as a compact text artifact.
pub fn inference_summary_text(
  summary : InferenceSummary,
  contrasts : Array[InferenceContrast],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("contrasts=")
  builder.write_string(summary.contrasts.to_string())
  builder.write_string("\nfinite=")
  builder.write_string(summary.finite_contrasts.to_string())
  builder.write_string("\nsignificant=")
  builder.write_string(summary.significant_contrasts.to_string())
  builder.write_string("\n")
  for contrast in contrasts {
    builder.write_string(contrast.name)
    builder.write_string(";estimate=")
    builder.write_string(contrast.estimate.to_string())
    builder.write_string(";p=")
    builder.write_string(contrast.p_value.to_string())
    builder.write_string(";significant=")
    builder.write_string(contrast.significant.to_string())
    builder.write_string("\n")
  }
  builder.to_string()
}

///|
/// Returns a compact power-curve summary vector.
pub fn inference_power_summary(
  curve : Array[InferencePowerPoint],
) -> Array[Double] {
  let first = inference_first_power_sample(curve)
  let mut maximum = 0.0
  for point in curve {
    if point.power > maximum {
      maximum = point.power
    }
  }
  [
    curve.length().to_double(),
    first.to_double(),
    maximum,
    if first > 0 {
      1.0
    } else {
      0.0
    },
  ]
}

///|
/// Computes an interval for a difference in two independent means.
pub fn inference_two_mean_contrast(
  treated_mean : Double,
  control_mean : Double,
  treated_standard_error : Double,
  control_standard_error : Double,
  confidence_level? : Double = 0.95,
) -> InferenceInterval {
  let estimate = treated_mean - control_mean
  let standard_error = (treated_standard_error * treated_standard_error +
  control_standard_error * control_standard_error).sqrt()
  inference_interval(estimate, standard_error, confidence_level~)
}

///|
/// Computes Welch's degrees-of-freedom approximation.
pub fn inference_welch_degrees_of_freedom(
  treated_variance : Double,
  treated_size : Int,
  control_variance : Double,
  control_size : Int,
) -> Double {
  let treated_term = treated_variance.abs() / treated_size.max(1).to_double()
  let control_term = control_variance.abs() / control_size.max(1).to_double()
  let numerator = (treated_term + control_term) * (treated_term + control_term)
  let denominator = treated_term *
    treated_term /
    (treated_size.max(1) - 1).max(1).to_double() +
    control_term * control_term / (control_size.max(1) - 1).max(1).to_double()
  if denominator <= 1.0e-12 {
    1.0
  } else {
    numerator / denominator
  }
}

///|
/// Computes a standardized mean contrast using a pooled deviation.
pub fn inference_standardized_contrast(
  estimate : Double,
  treated_variance : Double,
  control_variance : Double,
) -> Double {
  estimate /
  ((treated_variance + control_variance) / 2.0).abs().sqrt().max(1.0e-12)
}

///|
/// Computes a finite-safe ratio of two effects and a delta-method standard error.
pub fn inference_ratio_effect(
  numerator : Double,
  denominator : Double,
  numerator_standard_error : Double,
  denominator_standard_error : Double,
  confidence_level? : Double = 0.95,
) -> InferenceInterval {
  let safe_denominator = if denominator.abs() <= 1.0e-12 {
    if denominator < 0.0 {
      -1.0e-12
    } else {
      1.0e-12
    }
  } else {
    denominator
  }
  let ratio = numerator / safe_denominator
  let first = numerator_standard_error / safe_denominator
  let second = numerator *
    denominator_standard_error /
    (safe_denominator * safe_denominator)
  let standard_error = (first * first + second * second).sqrt()
  inference_interval(ratio, standard_error, confidence_level~)
}

///|
/// Computes heterogeneity across subgroup estimates.
pub fn inference_heterogeneity(
  estimates : Array[Double],
  standard_errors : Array[Double],
) -> Array[Double] {
  let n = estimates.length().min(standard_errors.length())
  if n == 0 {
    return [0.0, 0.0, 0.0]
  }
  let weights = Array::new(capacity=n)
  for i in 0..<n {
    weights.push(1.0 / @math.pow(standard_errors[i].abs().max(1.0e-12), 2.0))
  }
  let center = weighted_mean(estimates[:n].to_owned(), weights)
  let mut statistic = 0.0
  for i in 0..<n {
    statistic += weights[i] * (estimates[i] - center) * (estimates[i] - center)
  }
  let degrees = (n - 1).max(1).to_double()
  [center, statistic, statistic / degrees]
}

///|
/// Converts a list of intervals into their widths.
pub fn inference_interval_widths(
  intervals : Array[InferenceInterval],
) -> Array[Double] {
  intervals.map(fn(interval) { interval.width })
}

///|
/// Counts intervals that contain a specified target.
pub fn inference_containing_count(
  intervals : Array[InferenceInterval],
  target : Double,
) -> Int {
  intervals.fold(init=0, fn(total, interval) {
    if interval.lower <= target && target <= interval.upper {
      total + 1
    } else {
      total
    }
  })
}
