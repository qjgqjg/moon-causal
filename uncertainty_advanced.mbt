///|
/// Multiple-testing adjustment result.
pub struct MultipleTestingResult {
  raw_p_values : Array[Double]
  adjusted_p_values : Array[Double]
  rejected : Array[Bool]
  alpha : Double
  procedure : String
}

///|
/// Bootstrap interval with bias and Monte Carlo diagnostics.
pub struct BootstrapInterval {
  estimate : Double
  lower : Double
  upper : Double
  bias : Double
  standard_error : Double
  replicates : Int
  confidence_level : Double
  passes : Bool
}

///|
/// Cluster bootstrap summary.
pub struct ClusterBootstrapResult {
  estimate : Double
  replicates : Array[Double]
  interval : BootstrapInterval
  cluster_count : Int
  average_cluster_size : Double
}

///|
/// Influence-function summary for a scalar estimand.
pub struct InfluenceSummary {
  estimate : Double
  influence : Array[Double]
  standard_error : Double
  maximum_absolute_influence : Double
  effective_sample_size : Double
  high_influence_count : Int
}

///|
/// Sandwich covariance result for a coefficient vector.
pub struct SandwichResult {
  covariance : Array[Array[Double]]
  standard_errors : Array[Double]
  condition_proxy : Double
  robust : Bool
}

///|
/// Delta-method result for a transformed scalar estimand.
pub struct DeltaMethodResult {
  estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
  derivative : Double
  confidence_level : Double
}

///|
/// Simultaneous confidence band for a curve.
pub struct SimultaneousBand {
  estimate : Array[Double]
  lower : Array[Double]
  upper : Array[Double]
  critical_value : Double
  confidence_level : Double
}

///|
fn ua_normal_cdf(value : Double) -> Double {
  let sign = if value < 0.0 { -1.0 } else { 1.0 }
  let x = value.abs()
  let t = 1.0 / (1.0 + 0.2316419 * x)
  let polynomial = t *
    (
      0.319381530 +
      t *
      (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
    )
  let density = @math.exp(-0.5 * x * x) / 2.5066282746310002
  0.5 + sign * (0.5 - density * polynomial)
}

///|
fn ua_normal_quantile(probability : Double) -> Double {
  let p = clamp(probability, 1.0e-7, 1.0 - 1.0e-7)
  let mut low = -9.0
  let mut high = 9.0
  for _ in 0..<70 {
    let middle = (low + high) / 2.0
    if ua_normal_cdf(middle) < p {
      low = middle
    } else {
      high = middle
    }
  }
  (low + high) / 2.0
}

///|
fn ua_sorted(values : Array[Double]) -> Array[Double] {
  let result = values.copy()
  for i in 1..<result.length() {
    let value = result[i]
    let mut cursor = i
    while cursor > 0 && result[cursor - 1] > value {
      result[cursor] = result[cursor - 1]
      cursor -= 1
    }
    result[cursor] = value
  }
  result
}

///|
fn ua_quantile(values : Array[Double], probability : Double) -> Double {
  let sorted = ua_sorted(values)
  if sorted.length() == 0 {
    return 0.0
  }
  let position = clamp(probability, 0.0, 1.0) *
    (sorted.length() - 1).to_double()
  let lower = position.to_int()
  let upper = if lower + 1 < sorted.length() { lower + 1 } else { lower }
  sorted[lower] +
  (sorted[upper] - sorted[lower]) * (position - lower.to_double())
}

///|
/// Converts a z statistic to a two-sided normal p-value.
pub fn normal_two_sided_p_value(z_statistic : Double) -> Double {
  2.0 * (1.0 - ua_normal_cdf(z_statistic.abs()))
}

///|
/// Computes a normal confidence interval for a scalar estimate.
pub fn normal_confidence_interval(
  estimate : Double,
  standard_error : Double,
  confidence_level? : Double = 0.95,
) -> Array[Double] {
  let level = clamp(confidence_level, 0.5, 0.999999)
  let critical = ua_normal_quantile(0.5 + level / 2.0)
  [
    estimate - critical * standard_error.abs(),
    estimate + critical * standard_error.abs(),
  ]
}

///|
/// Applies Bonferroni adjustment.
pub fn bonferroni_adjust(
  p_values : Array[Double],
  alpha? : Double = 0.05,
) -> MultipleTestingResult {
  let count = p_values.length().to_double().max(1.0)
  let adjusted = Array::new(capacity=p_values.length())
  let rejected = Array::new(capacity=p_values.length())
  let level = clamp(alpha, 1.0e-8, 1.0)
  for p_value in p_values {
    let adjusted_value = clamp(p_value, 0.0, 1.0) * count
    adjusted.push(adjusted_value.min(1.0))
    rejected.push(adjusted_value <= level)
  }
  {
    raw_p_values: p_values.copy(),
    adjusted_p_values: adjusted,
    rejected,
    alpha: level,
    procedure: "bonferroni",
  }
}

///|
/// Applies Holm step-down adjustment while preserving original ordering.
pub fn holm_adjust(
  p_values : Array[Double],
  alpha? : Double = 0.05,
) -> MultipleTestingResult {
  let n = p_values.length()
  let order : Array[Int] = Array::new(capacity=n)
  for i in 0..<n {
    order.push(i)
  }
  for i in 1..<n {
    let value = order[i]
    let mut cursor = i
    while cursor > 0 && p_values[order[cursor - 1]] > p_values[value] {
      order[cursor] = order[cursor - 1]
      cursor -= 1
    }
    order[cursor] = value
  }
  let adjusted = Array::make(n, 0.0)
  let level = clamp(alpha, 1.0e-8, 1.0)
  let mut running = 0.0
  for position in 0..<n {
    let index = order[position]
    let candidate = clamp(p_values[index], 0.0, 1.0) *
      (n - position).to_double()
    running = running.max(candidate)
    adjusted[index] = running.min(1.0)
  }
  let rejected = Array::new(capacity=n)
  for value in adjusted {
    rejected.push(value <= level)
  }
  {
    raw_p_values: p_values.copy(),
    adjusted_p_values: adjusted,
    rejected,
    alpha: level,
    procedure: "holm",
  }
}

///|
/// Applies Benjamini-Hochberg false-discovery-rate adjustment.
pub fn benjamini_hochberg(
  p_values : Array[Double],
  alpha? : Double = 0.05,
) -> MultipleTestingResult {
  let n = p_values.length()
  let order : Array[Int] = Array::new(capacity=n)
  for i in 0..<n {
    order.push(i)
  }
  for i in 1..<n {
    let value = order[i]
    let mut cursor = i
    while cursor > 0 && p_values[order[cursor - 1]] > p_values[value] {
      order[cursor] = order[cursor - 1]
      cursor -= 1
    }
    order[cursor] = value
  }
  let adjusted = Array::make(n, 1.0)
  let mut running = 1.0
  for reverse in 0..<n {
    let position = n - 1 - reverse
    let index = order[position]
    let candidate = clamp(p_values[index], 0.0, 1.0) *
      n.to_double() /
      (position + 1).to_double()
    running = running.min(candidate)
    adjusted[index] = running.min(1.0)
  }
  let level = clamp(alpha, 1.0e-8, 1.0)
  let rejected = Array::new(capacity=n)
  for value in adjusted {
    rejected.push(value <= level)
  }
  {
    raw_p_values: p_values.copy(),
    adjusted_p_values: adjusted,
    rejected,
    alpha: level,
    procedure: "benjamini-hochberg",
  }
}

///|
/// Counts discoveries in a multiple-testing result.
pub fn rejected_count(result : MultipleTestingResult) -> Int {
  let mut count = 0
  for value in result.rejected {
    if value {
      count += 1
    }
  }
  count
}

///|
/// Computes an empirical bootstrap interval.
pub fn bootstrap_interval(
  estimate : Double,
  replicates : Array[Double],
  confidence_level? : Double = 0.95,
) -> BootstrapInterval {
  let level = clamp(confidence_level, 0.5, 0.999999)
  let lower_probability = (1.0 - level) / 2.0
  let upper_probability = 1.0 - lower_probability
  let lower = ua_quantile(replicates, lower_probability)
  let upper = ua_quantile(replicates, upper_probability)
  let center = mean_or(replicates, estimate)
  let bias = center - estimate
  let standard_error = std_dev(replicates)
  {
    estimate,
    lower,
    upper,
    bias,
    standard_error,
    replicates: replicates.length(),
    confidence_level: level,
    passes: replicates.length() >= 50 && is_finite(lower) && is_finite(upper),
  }
}

///|
/// Computes a bias-corrected percentile interval.
pub fn bias_corrected_interval(
  estimate : Double,
  replicates : Array[Double],
  confidence_level? : Double = 0.95,
) -> BootstrapInterval {
  let below = replicates.filter(fn(value) { value < estimate }).length()
  let probability = if replicates.length() == 0 {
    0.5
  } else {
    below.to_double() / replicates.length().to_double()
  }
  let correction = ua_normal_quantile(probability)
  let level = clamp(confidence_level, 0.5, 0.999999)
  let alpha = ua_normal_quantile((1.0 - level) / 2.0)
  let upper_alpha = ua_normal_quantile(1.0 - (1.0 - level) / 2.0)
  let lower_probability = ua_normal_cdf(2.0 * correction + alpha)
  let upper_probability = ua_normal_cdf(2.0 * correction + upper_alpha)
  let lower = ua_quantile(replicates, lower_probability)
  let upper = ua_quantile(replicates, upper_probability)
  {
    estimate,
    lower,
    upper,
    bias: mean_or(replicates, estimate) - estimate,
    standard_error: std_dev(replicates),
    replicates: replicates.length(),
    confidence_level: level,
    passes: replicates.length() >= 50,
  }
}

///|
/// Resamples whole clusters and computes a mean-effect bootstrap distribution.
pub fn cluster_bootstrap_mean(
  treatment : Array[Bool],
  outcomes : Array[Double],
  cluster_ids : Array[Int],
  replicates : Int,
  seed : UInt64,
) -> ClusterBootstrapResult {
  let n = treatment.length().min(outcomes.length()).min(cluster_ids.length())
  let clusters : Array[Int] = Array::new()
  for cluster in cluster_ids[:n] {
    if !clusters.contains(cluster) {
      clusters.push(cluster)
    }
  }
  let rng = RandomState::new(seed)
  let draws : Array[Double] = Array::new(
    capacity=if replicates > 0 { replicates } else { 0 },
  )
  let original = estimate_difference_in_means(
      outcomes[:n].to_owned(),
      treatment[:n].to_owned(),
    ).estimate
  for _ in 0..<replicates {
    let sampled_outcomes : Array[Double] = Array::new()
    let sampled_treatment : Array[Bool] = Array::new()
    for _ in 0..<clusters.length() {
      if clusters.length() > 0 {
        let cluster = clusters[(rng.uniform() * clusters.length().to_double()).to_int()]
        for i in 0..<n {
          if cluster_ids[i] == cluster {
            sampled_outcomes.push(outcomes[i])
            sampled_treatment.push(treatment[i])
          }
        }
      }
    }
    draws.push(
      estimate_difference_in_means(sampled_outcomes, sampled_treatment).estimate,
    )
  }
  let interval = bootstrap_interval(original, draws)
  {
    estimate: original,
    replicates: draws,
    interval,
    cluster_count: clusters.length(),
    average_cluster_size: if clusters.length() == 0 {
      0.0
    } else {
      n.to_double() / clusters.length().to_double()
    },
  }
}

///|
/// Builds an influence summary for a vector of individual contributions.
pub fn influence_summary(
  estimate : Double,
  influence : Array[Double],
  threshold? : Double = 4.0,
) -> InfluenceSummary {
  let standard_error = if influence.length() == 0 {
    0.0
  } else {
    std_dev(influence) / influence.length().to_double().sqrt()
  }
  let mut maximum = 0.0
  let mut high = 0
  let cutoff = threshold * standard_error.max(1.0e-12)
  for value in influence {
    if value.abs() > maximum {
      maximum = value.abs()
    }
    if value.abs() > cutoff {
      high += 1
    }
  }
  let weights = Array::make(influence.length(), 1.0)
  {
    estimate,
    influence,
    standard_error,
    maximum_absolute_influence: maximum,
    effective_sample_size: effective_sample_size(weights),
    high_influence_count: high,
  }
}

///|
/// Computes a diagonal sandwich covariance from score rows.
pub fn diagonal_sandwich(
  scores : Array[Array[Double]],
  bread_diagonal : Array[Double],
  robust? : Bool = true,
) -> SandwichResult {
  let width = bread_diagonal.length()
  let covariance = Array::make(width, Array::make(width, 0.0))
  if scores.length() == 0 {
    return {
      covariance,
      standard_errors: Array::make(width, 0.0),
      condition_proxy: 0.0,
      robust,
    }
  }
  for score in scores {
    for j in 0..<width.min(score.length()) {
      covariance[j][j] += score[j] * score[j]
    }
  }
  let standard_errors = Array::new(capacity=width)
  let mut condition_proxy = 0.0
  for j in 0..<width {
    let scale = bread_diagonal[j].abs().max(1.0e-12)
    covariance[j][j] /= scores.length().to_double() * scale * scale
    standard_errors.push(covariance[j][j].max(0.0).sqrt())
    if scale > condition_proxy {
      condition_proxy = scale
    }
  }
  { covariance, standard_errors, condition_proxy, robust }
}

///|
/// Computes a coefficient covariance from residuals and a design matrix.
pub fn regression_sandwich(
  design : Array[Array[Double]],
  residuals : Array[Double],
  ridge? : Double = 1.0e-8,
) -> SandwichResult {
  let n = design.length().min(residuals.length())
  if n == 0 {
    return {
      covariance: [],
      standard_errors: [],
      condition_proxy: 0.0,
      robust: true,
    }
  }
  let width = design[0].length()
  let bread = add_ridge(
    matrix_multiply(matrix_transpose(design), design),
    ridge,
  )
  let inverse = matrix_inverse(bread)
  if inverse.length() == 0 {
    return {
      covariance: [],
      standard_errors: [],
      condition_proxy: 0.0,
      robust: true,
    }
  }
  let meat = Array::make(width, Array::make(width, 0.0))
  for i in 0..<n {
    for j in 0..<width {
      for k in 0..<width {
        meat[j][k] += design[i][j] * residuals[i] * design[i][k] * residuals[i]
      }
    }
  }
  let covariance = matrix_multiply(matrix_multiply(inverse, meat), inverse)
  let errors = Array::new(capacity=width)
  for j in 0..<width {
    errors.push(covariance[j][j].max(0.0).sqrt())
  }
  { covariance, standard_errors: errors, condition_proxy: 1.0, robust: true }
}

///|
/// Computes a scalar delta-method confidence interval.
pub fn delta_method(
  estimate : Double,
  variance_value : Double,
  derivative : Double,
  confidence_level? : Double = 0.95,
) -> DeltaMethodResult {
  let standard_error = (variance_value.max(0.0) * derivative * derivative).sqrt()
  let interval = normal_confidence_interval(
    estimate,
    standard_error,
    confidence_level~,
  )
  {
    estimate,
    standard_error,
    lower: interval[0],
    upper: interval[1],
    derivative,
    confidence_level: clamp(confidence_level, 0.5, 0.999999),
  }
}

///|
/// Computes a log-ratio delta-method result from two independent estimates.
pub fn log_ratio_delta(
  numerator : Double,
  denominator : Double,
  numerator_variance : Double,
  denominator_variance : Double,
  confidence_level? : Double = 0.95,
) -> DeltaMethodResult {
  let safe_numerator = numerator.abs().max(1.0e-12)
  let safe_denominator = denominator.abs().max(1.0e-12)
  let value = @math.ln(safe_numerator / safe_denominator)
  let variance_value = numerator_variance / (safe_numerator * safe_numerator) +
    denominator_variance / (safe_denominator * safe_denominator)
  delta_method(value, variance_value, 1.0, confidence_level~)
}

///|
/// Builds a simultaneous band by a conservative max-z critical value.
pub fn simultaneous_band(
  estimates : Array[Double],
  standard_errors : Array[Double],
  confidence_level? : Double = 0.95,
) -> SimultaneousBand {
  let level = clamp(confidence_level, 0.5, 0.999999)
  let width = estimates.length().min(standard_errors.length())
  let critical = ua_normal_quantile(
    0.5 + level / (2.0 * width.max(1).to_double()),
  )
  let lower = Array::new(capacity=width)
  let upper = Array::new(capacity=width)
  let estimate = estimates[:width].to_owned()
  for i in 0..<width {
    let radius = critical * standard_errors[i].abs()
    lower.push(estimates[i] - radius)
    upper.push(estimates[i] + radius)
  }
  { estimate, lower, upper, critical_value: critical, confidence_level: level }
}

///|
/// Computes a wild bootstrap distribution using Rademacher multipliers.
pub fn wild_bootstrap_mean(
  values : Array[Double],
  replicates : Int,
  seed : UInt64,
) -> Array[Double] {
  let rng = RandomState::new(seed)
  let center = mean_or(values, 0.0)
  let residuals = Array::new(capacity=values.length())
  for value in values {
    residuals.push(value - center)
  }
  let result = Array::new(capacity=if replicates > 0 { replicates } else { 0 })
  for _ in 0..<replicates {
    let mut total = 0.0
    for residual in residuals {
      total += if rng.uniform() < 0.5 { -residual } else { residual }
    }
    let contribution = if residuals.length() == 0 {
      0.0
    } else {
      total / residuals.length().to_double()
    }
    result.push(center + contribution)
  }
  result
}

///|
/// Computes a percentile interval for a mean under wild bootstrap.
pub fn wild_bootstrap_interval(
  values : Array[Double],
  replicates : Int,
  seed : UInt64,
  confidence_level? : Double = 0.95,
) -> BootstrapInterval {
  let estimate = mean_or(values, 0.0)
  bootstrap_interval(
    estimate,
    wild_bootstrap_mean(values, replicates, seed),
    confidence_level~,
  )
}

///|
/// Computes a leave-one-out jackknife estimate of a mean.
pub fn jackknife_means(values : Array[Double]) -> Array[Double] {
  let n = values.length()
  let result = Array::new(capacity=n)
  if n == 0 {
    return result
  }
  let total = sum(values)
  for i in 0..<n {
    result.push(
      if n == 1 {
        0.0
      } else {
        (total - values[i]) / (n - 1).to_double()
      },
    )
  }
  result
}

///|
/// Computes the jackknife acceleration constant.
pub fn jackknife_acceleration(values : Array[Double]) -> Double {
  let jackknife = jackknife_means(values)
  if jackknife.length() < 2 {
    return 0.0
  }
  let center = mean(jackknife)
  let mut numerator = 0.0
  let mut denominator = 0.0
  for value in jackknife {
    let difference = center - value
    numerator += difference * difference * difference
    denominator += difference * difference
  }
  if denominator == 0.0 {
    0.0
  } else {
    numerator / (6.0 * @math.pow(denominator, 1.5))
  }
}

///|
/// Computes a robust median absolute deviation scale.
fn ua_median_absolute_deviation(values : Array[Double]) -> Double {
  let center = ua_quantile(values, 0.5)
  let deviations = Array::new(capacity=values.length())
  for value in values {
    deviations.push((value - center).abs())
  }
  1.4826 * ua_quantile(deviations, 0.5)
}

///|
/// Computes Huber weights for robust standard-error monitoring.
pub fn huber_weights(
  residuals : Array[Double],
  tuning? : Double = 1.345,
) -> Array[Double] {
  let scale = ua_median_absolute_deviation(residuals).max(1.0e-12)
  let cutoff = tuning * scale
  let result = Array::new(capacity=residuals.length())
  for residual in residuals {
    let magnitude = residual.abs()
    result.push(if magnitude <= cutoff { 1.0 } else { cutoff / magnitude })
  }
  result
}

///|
/// Produces a compact uncertainty summary vector.
pub fn uncertainty_summary(result : BootstrapInterval) -> Array[Double] {
  [
    result.estimate,
    result.lower,
    result.upper,
    result.bias,
    result.standard_error,
    result.replicates.to_double(),
    result.confidence_level,
    if result.passes {
      1.0
    } else {
      0.0
    },
  ]
}
