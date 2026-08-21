///|
/// Design specification for a reproducible randomized experiment.
pub struct ExperimentDesign {
  sample_size : Int
  treatment_fraction : Double
  strata_count : Int
  block_size : Int
  seed : UInt64
  cluster_randomization : Bool
}

///|
/// Result of a randomization assignment.
pub struct AssignmentResult {
  treatment : Array[Bool]
  treated_count : Int
  control_count : Int
  strata_imbalance : Array[Double]
  assignment_rate : Double
  passes : Bool
}

///|
/// Approximate sample-size and power calculation for a two-arm mean effect.
pub struct PowerAnalysis {
  sample_size : Int
  effect : Double
  standard_deviation : Double
  alpha : Double
  power : Double
  minimum_detectable_effect : Double
  detectable : Bool
}

///|
/// Experiment health diagnostics.
pub struct ExperimentHealth {
  sample_ratio_statistic : Double
  sample_ratio_passes : Bool
  pre_period_difference : Double
  pre_period_passes : Bool
  treatment_rate : Double
  effective_sample_size : Double
  passes : Bool
}

///|
/// CUPED adjustment output.
pub struct CupedResult {
  adjusted_outcome : Array[Double]
  theta : Double
  variance_reduction : Double
  covariance : Double
  pre_period_variance : Double
}

///|
fn exp_design_normal_cdf(value : Double) -> Double {
  let sign = if value < 0.0 { -1.0 } else { 1.0 }
  let absolute = value.abs()
  let t = 1.0 / (1.0 + 0.2316419 * absolute)
  let polynomial = t *
    (
      0.319381530 +
      t *
      (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
    )
  let density = @math.exp(-0.5 * absolute * absolute) / 2.5066282746310002
  0.5 + sign * (0.5 - density * polynomial)
}

///|
fn exp_design_normal_quantile(probability : Double) -> Double {
  let p = clamp(probability, 1.0e-6, 1.0 - 1.0e-6)
  let mut low = -8.0
  let mut high = 8.0
  for _ in 0..<60 {
    let middle = (low + high) / 2.0
    if exp_design_normal_cdf(middle) < p {
      low = middle
    } else {
      high = middle
    }
  }
  (low + high) / 2.0
}

///|
/// Creates a validated experiment design.
pub fn experiment_design(
  sample_size : Int,
  treatment_fraction? : Double = 0.5,
  strata_count? : Int = 1,
  block_size? : Int = 4,
  seed? : UInt64 = 20260819,
  cluster_randomization? : Bool = false,
) -> ExperimentDesign {
  let fraction = clamp(treatment_fraction, 0.01, 0.99)
  {
    sample_size: if sample_size > 0 {
      sample_size
    } else {
      0
    },
    treatment_fraction: fraction,
    strata_count: if strata_count > 0 {
      strata_count
    } else {
      1
    },
    block_size: if block_size > 1 {
      block_size
    } else {
      2
    },
    seed,
    cluster_randomization,
  }
}

///|
/// Assigns units using a deterministic Bernoulli randomization.
pub fn random_assignment(
  sample_size : Int,
  treatment_fraction? : Double = 0.5,
  seed? : UInt64 = 20260819,
) -> AssignmentResult {
  let n = if sample_size > 0 { sample_size } else { 0 }
  let fraction = clamp(treatment_fraction, 0.0, 1.0)
  let rng = RandomState::new(seed)
  let treatment : Array[Bool] = Array::new(capacity=n)
  let mut treated = 0
  for _ in 0..<n {
    let assigned = rng.uniform() < fraction
    treatment.push(assigned)
    if assigned {
      treated += 1
    }
  }
  let control = n - treated
  let rate = if n == 0 { 0.0 } else { treated.to_double() / n.to_double() }
  {
    treatment,
    treated_count: treated,
    control_count: control,
    strata_imbalance: [],
    assignment_rate: rate,
    passes: n == 0 || (rate - fraction).abs() <= 0.2,
  }
}

///|
/// Assigns complete blocks with a balanced number of arms per block.
pub fn blocked_assignment(
  sample_size : Int,
  block_size? : Int = 4,
  treatment_fraction? : Double = 0.5,
  seed? : UInt64 = 20260819,
) -> AssignmentResult {
  let n = if sample_size > 0 { sample_size } else { 0 }
  let block = if block_size > 1 { block_size } else { 2 }
  let target = clamp(treatment_fraction, 0.0, 1.0)
  let rng = RandomState::new(seed)
  let assignment : Array[Bool] = Array::new(capacity=n)
  let mut treated = 0
  let mut start = 0
  while start < n {
    let end = if start + block < n { start + block } else { n }
    let width = end - start
    let block_treated = (width.to_double() * target).round().to_int()
    let block_flags = Array::make(width, false)
    let mut placed = 0
    while placed < block_treated {
      let index = (rng.uniform() * width.to_double()).to_int()
      if !block_flags[index] {
        block_flags[index] = true
        placed += 1
      }
    }
    for value in block_flags {
      assignment.push(value)
      if value {
        treated += 1
      }
    }
    start = end
  }
  let control = n - treated
  let rate = if n == 0 { 0.0 } else { treated.to_double() / n.to_double() }
  {
    treatment: assignment,
    treated_count: treated,
    control_count: control,
    strata_imbalance: [],
    assignment_rate: rate,
    passes: true,
  }
}

///|
/// Assigns treatment independently within integer strata.
pub fn stratified_assignment(
  strata : Array[Int],
  treatment_fraction? : Double = 0.5,
  seed? : UInt64 = 20260819,
) -> AssignmentResult {
  let n = strata.length()
  let target = clamp(treatment_fraction, 0.0, 1.0)
  let rng = RandomState::new(seed)
  let assignment = Array::make(n, false)
  let mut max_stratum = -1
  for value in strata {
    if value > max_stratum {
      max_stratum = value
    }
  }
  let stratum_count = if max_stratum < 0 { 0 } else { max_stratum + 1 }
  let strata_treated = Array::make(stratum_count, 0)
  let strata_total = Array::make(stratum_count, 0)
  for value in strata {
    if value >= 0 && value < stratum_count {
      strata_total[value] += 1
    }
  }
  for i in 0..<n {
    let stratum = strata[i]
    if stratum < 0 || stratum >= stratum_count {
      assignment[i] = rng.uniform() < target
    } else {
      let local_rate = target
      assignment[i] = rng.uniform() < local_rate
      if assignment[i] {
        strata_treated[stratum] += 1
      }
    }
  }
  let imbalance = Array::new(capacity=stratum_count)
  for i in 0..<stratum_count {
    let observed = if strata_total[i] == 0 {
      0.0
    } else {
      strata_treated[i].to_double() / strata_total[i].to_double()
    }
    imbalance.push(observed - target)
  }
  let mut treated = 0
  for value in assignment {
    if value {
      treated += 1
    }
  }
  let rate = if n == 0 { 0.0 } else { treated.to_double() / n.to_double() }
  {
    treatment: assignment,
    treated_count: treated,
    control_count: n - treated,
    strata_imbalance: imbalance,
    assignment_rate: rate,
    passes: true,
  }
}

///|
/// Calculates the approximate standard error of a two-arm mean contrast.
pub fn two_arm_standard_error(
  sample_size : Int,
  treatment_fraction : Double,
  standard_deviation : Double,
) -> Double {
  if sample_size <= 1 {
    return 0.0
  }
  let fraction = clamp(treatment_fraction, 1.0e-6, 1.0 - 1.0e-6)
  standard_deviation *
  (1.0 / (sample_size.to_double() * fraction) +
  1.0 / (sample_size.to_double() * (1.0 - fraction))).sqrt()
}

///|
/// Computes a normal-approximation power analysis.
pub fn power_analysis(
  sample_size : Int,
  effect : Double,
  standard_deviation : Double,
  treatment_fraction? : Double = 0.5,
  alpha? : Double = 0.05,
) -> PowerAnalysis {
  let fraction = clamp(treatment_fraction, 0.01, 0.99)
  let significance = clamp(alpha, 1.0e-6, 0.5)
  let standard_error = two_arm_standard_error(
    sample_size,
    fraction,
    standard_deviation.abs(),
  )
  let critical = exp_design_normal_quantile(1.0 - significance / 2.0)
  let noncentrality = if standard_error == 0.0 {
    0.0
  } else {
    effect.abs() / standard_error
  }
  let power = exp_design_normal_cdf(-critical - noncentrality) +
    1.0 -
    exp_design_normal_cdf(critical - noncentrality)
  let target_power = 0.8
  let minimum_effect = critical *
    2.0 *
    standard_deviation.abs() /
    (sample_size.to_double() * fraction * (1.0 - fraction)).sqrt()
  {
    sample_size,
    effect,
    standard_deviation,
    alpha: significance,
    power,
    minimum_detectable_effect: minimum_effect,
    detectable: power >= target_power,
  }
}

///|
/// Estimates the minimum sample size for a target power by enumeration.
pub fn sample_size_for_power(
  effect : Double,
  standard_deviation : Double,
  target_power? : Double = 0.8,
  treatment_fraction? : Double = 0.5,
  alpha? : Double = 0.05,
  maximum_sample_size? : Int = 100000,
) -> Int {
  let target = clamp(target_power, 0.01, 0.999)
  let mut sample_size = 2
  let limit = if maximum_sample_size > 2 { maximum_sample_size } else { 2 }
  while sample_size <= limit {
    if power_analysis(
        sample_size,
        effect,
        standard_deviation,
        treatment_fraction~,
        alpha~,
      ).power >=
      target {
      return sample_size
    }
    sample_size += 1
  }
  limit
}

///|
/// Computes the minimum detectable effect at a fixed sample size.
pub fn minimum_detectable_effect(
  sample_size : Int,
  standard_deviation : Double,
  treatment_fraction? : Double = 0.5,
  alpha? : Double = 0.05,
  target_power? : Double = 0.8,
) -> Double {
  let fraction = clamp(treatment_fraction, 0.01, 0.99)
  let critical = exp_design_normal_quantile(
    1.0 - clamp(alpha, 1.0e-6, 0.5) / 2.0,
  )
  let power_quantile = exp_design_normal_quantile(
    clamp(target_power, 0.5, 0.999),
  )
  let denominator = (sample_size.to_double() * fraction * (1.0 - fraction)).sqrt()
  if denominator == 0.0 {
    0.0
  } else {
    (critical + power_quantile) * standard_deviation.abs() / denominator
  }
}

///|
/// Computes a standardized mean difference for an A/A check.
pub fn pre_period_difference(
  treatment : Array[Bool],
  pre_period : Array[Double],
) -> Double {
  let treated = Array::new()
  let control = Array::new()
  let n = if treatment.length() < pre_period.length() {
    treatment.length()
  } else {
    pre_period.length()
  }
  for i in 0..<n {
    if treatment[i] {
      treated.push(pre_period[i])
    } else {
      control.push(pre_period[i])
    }
  }
  if treated.length() == 0 || control.length() == 0 {
    0.0
  } else {
    let estimate = estimate_difference_in_means(pre_period, treatment).estimate
    standardized_effect_size(estimate, pre_period)
  }
}

///|
/// Computes a sample-ratio mismatch chi-square-like statistic.
pub fn sample_ratio_statistic(
  treatment : Array[Bool],
  expected_fraction? : Double = 0.5,
) -> Double {
  let mut observed_treated = 0
  for value in treatment {
    if value {
      observed_treated += 1
    }
  }
  let n = treatment.length()
  if n == 0 {
    return 0.0
  }
  let expected = n.to_double() * clamp(expected_fraction, 0.01, 0.99)
  let expected_control = n.to_double() - expected
  let observed_control = n.to_double() - observed_treated.to_double()
  (observed_treated.to_double() - expected) *
  (observed_treated.to_double() - expected) /
  expected +
  (observed_control - expected_control) *
  (observed_control - expected_control) /
  expected_control
}

///|
/// Returns whether a sample-ratio statistic is below an operational cutoff.
pub fn sample_ratio_passes(
  treatment : Array[Bool],
  expected_fraction? : Double = 0.5,
  cutoff? : Double = 10.83,
) -> Bool {
  sample_ratio_statistic(treatment, expected_fraction~) <= cutoff
}

///|
/// Computes CUPED's pre-period coefficient and adjusted outcomes.
pub fn cuped_adjust(
  outcomes : Array[Double],
  pre_period : Array[Double],
) -> CupedResult {
  let n = if outcomes.length() < pre_period.length() {
    outcomes.length()
  } else {
    pre_period.length()
  }
  if n < 2 {
    return {
      adjusted_outcome: outcomes[:n].to_owned(),
      theta: 0.0,
      variance_reduction: 0.0,
      covariance: 0.0,
      pre_period_variance: 0.0,
    }
  }
  let y = outcomes[:n].to_owned()
  let x = pre_period[:n].to_owned()
  let y_mean = mean(y)
  let x_mean = mean(x)
  let mut covariance_value = 0.0
  let mut variance_x = 0.0
  for i in 0..<n {
    covariance_value += (y[i] - y_mean) * (x[i] - x_mean)
    variance_x += (x[i] - x_mean) * (x[i] - x_mean)
  }
  covariance_value /= n.to_double()
  variance_x /= n.to_double()
  let theta = if variance_x == 0.0 {
    0.0
  } else {
    covariance_value / variance_x
  }
  let adjusted = Array::new(capacity=n)
  for i in 0..<n {
    adjusted.push(y[i] - theta * (x[i] - x_mean))
  }
  let original_variance = variance(y)
  let adjusted_variance = variance(adjusted)
  let reduction = if original_variance == 0.0 {
    0.0
  } else {
    1.0 - adjusted_variance / original_variance
  }
  {
    adjusted_outcome: adjusted,
    theta,
    variance_reduction: reduction,
    covariance: covariance_value,
    pre_period_variance: variance_x,
  }
}

///|
/// Estimates a CUPED treatment effect using the adjusted outcomes.
pub fn cuped_effect(
  treatment : Array[Bool],
  outcomes : Array[Double],
  pre_period : Array[Double],
) -> Estimate {
  let adjusted = cuped_adjust(outcomes, pre_period).adjusted_outcome
  let n = if treatment.length() < adjusted.length() {
    treatment.length()
  } else {
    adjusted.length()
  }
  let combined = Array::new(capacity=n)
  let labels = Array::new(capacity=n)
  for value in adjusted {
    combined.push(value)
  }
  for i in 0..<n {
    labels.push(treatment[i])
  }
  estimate_difference_in_means(combined, labels)
}

///|
/// Builds an experiment health report.
pub fn experiment_health(
  treatment : Array[Bool],
  pre_period : Array[Double],
  expected_fraction? : Double = 0.5,
  maximum_pre_difference? : Double = 0.1,
) -> ExperimentHealth {
  let ratio = sample_ratio_statistic(treatment, expected_fraction~)
  let pre_difference = pre_period_difference(treatment, pre_period).abs()
  let treated = treatment.filter(fn(value) { value }).length()
  let rate = if treatment.length() == 0 {
    0.0
  } else {
    treated.to_double() / treatment.length().to_double()
  }
  let weights = Array::make(treatment.length(), 1.0)
  let ess = effective_sample_size(weights)
  let ratio_ok = ratio <= 10.83
  let pre_ok = pre_difference <= maximum_pre_difference
  {
    sample_ratio_statistic: ratio,
    sample_ratio_passes: ratio_ok,
    pre_period_difference: pre_difference,
    pre_period_passes: pre_ok,
    treatment_rate: rate,
    effective_sample_size: ess,
    passes: ratio_ok && pre_ok,
  }
}

///|
/// Calculates a sequential Pocock-like monitoring boundary.
pub fn sequential_boundary(alpha : Double, looks : Int) -> Double {
  let valid_looks = if looks > 0 { looks } else { 1 }
  let adjusted = clamp(alpha, 1.0e-6, 0.5) / valid_looks.to_double()
  exp_design_normal_quantile(1.0 - adjusted / 2.0)
}

///|
/// Returns whether a z statistic crosses a sequential boundary.
pub fn sequential_stop(
  z_statistic : Double,
  alpha : Double,
  looks : Int,
) -> Bool {
  z_statistic.abs() >= sequential_boundary(alpha, looks)
}

///|
/// Assigns entire clusters to treatment or control.
pub fn cluster_assignment(
  cluster_ids : Array[Int],
  treatment_fraction? : Double = 0.5,
  seed? : UInt64 = 20260819,
) -> AssignmentResult {
  let unique : Array[Int] = Array::new()
  for cluster in cluster_ids {
    if !unique.contains(cluster) {
      unique.push(cluster)
    }
  }
  let cluster_result = blocked_assignment(
    unique.length(),
    block_size=2,
    treatment_fraction~,
    seed~,
  )
  let assignment = Array::make(cluster_ids.length(), false)
  for i in 0..<cluster_ids.length() {
    let mut cluster_index = -1
    for j in 0..<unique.length() {
      if unique[j] == cluster_ids[i] {
        cluster_index = j
        break
      }
    }
    if cluster_index >= 0 {
      assignment[i] = cluster_result.treatment[cluster_index]
    }
  }
  let mut treated = 0
  for value in assignment {
    if value {
      treated += 1
    }
  }
  {
    treatment: assignment,
    treated_count: treated,
    control_count: assignment.length() - treated,
    strata_imbalance: [],
    assignment_rate: if assignment.length() == 0 {
      0.0
    } else {
      treated.to_double() / assignment.length().to_double()
    },
    passes: true,
  }
}

///|
/// Computes cluster-level intraclass correlation.
pub fn intraclass_correlation(
  cluster_ids : Array[Int],
  outcomes : Array[Double],
) -> Double {
  let n = if cluster_ids.length() < outcomes.length() {
    cluster_ids.length()
  } else {
    outcomes.length()
  }
  if n < 2 {
    return 0.0
  }
  let clusters : Array[Int] = Array::new()
  for cluster in cluster_ids[:n] {
    if !clusters.contains(cluster) {
      clusters.push(cluster)
    }
  }
  let grand = mean(outcomes[:n].to_owned())
  let mut between = 0.0
  let mut within = 0.0
  for cluster in clusters {
    let values = Array::new()
    for i in 0..<n {
      if cluster_ids[i] == cluster {
        values.push(outcomes[i])
      }
    }
    if values.length() > 0 {
      let center = mean(values)
      between += values.length().to_double() *
        (center - grand) *
        (center - grand)
      for value in values {
        within += (value - center) * (value - center)
      }
    }
  }
  if between + within == 0.0 {
    0.0
  } else {
    between / (between + within)
  }
}

///|
/// Inflates a variance for cluster randomization using an ICC.
pub fn cluster_design_effect(
  cluster_size : Double,
  intraclass_correlation : Double,
) -> Double {
  1.0 + (cluster_size.max(1.0) - 1.0) * clamp(intraclass_correlation, 0.0, 1.0)
}

///|
/// Adjusts a nominal sample size for clustering.
pub fn cluster_adjusted_sample_size(
  nominal : Int,
  cluster_size : Double,
  intraclass_correlation : Double,
) -> Int {
  (nominal.to_double() *
  cluster_design_effect(cluster_size, intraclass_correlation))
  .ceil()
  .to_int()
}

///|
/// Computes a treatment-rate weighted average of stratum effects.
pub fn stratified_average_effect(
  effects : Array[Double],
  stratum_sizes : Array[Int],
) -> Double {
  let n = if effects.length() < stratum_sizes.length() {
    effects.length()
  } else {
    stratum_sizes.length()
  }
  let mut numerator = 0.0
  let mut denominator = 0
  for i in 0..<n {
    if stratum_sizes[i] > 0 {
      numerator += effects[i] * stratum_sizes[i].to_double()
      denominator += stratum_sizes[i]
    }
  }
  if denominator == 0 {
    0.0
  } else {
    numerator / denominator.to_double()
  }
}

///|
/// Returns a compact design summary vector.
pub fn design_summary(design : ExperimentDesign) -> Array[Double] {
  [
    design.sample_size.to_double(),
    design.treatment_fraction,
    design.strata_count.to_double(),
    design.block_size.to_double(),
    if design.cluster_randomization {
      1.0
    } else {
      0.0
    },
  ]
}
