///|
/// Survey design metadata.
pub struct SurveyDesign {
  weights : Array[Double]
  strata : Array[Int]
  clusters : Array[Int]
  population_size : Int
  sample_size : Int
  with_replacement : Bool
}

///|
/// Survey-weighted causal estimate.
pub struct SurveyEffect {
  estimate : Double
  standard_error : Double
  effective_sample_size : Double
  design_effect : Double
  finite_population_correction : Double
  passes : Bool
}

///|
/// Post-stratification cell summary.
pub struct PostStratum {
  key : Int
  sample_count : Int
  population_count : Int
  sample_mean : Double
  weighted_contribution : Double
}

///|
/// Calibration-weight result.
pub struct CalibrationResult {
  weights : Array[Double]
  iterations : Int
  maximum_margin_error : Double
  effective_sample_size : Double
  converged : Bool
}

///|
/// Creates survey design metadata with safe population defaults.
pub fn survey_design(
  weights : Array[Double],
  strata? : Array[Int] = [],
  clusters? : Array[Int] = [],
  population_size? : Int = 0,
  with_replacement? : Bool = false,
) -> SurveyDesign {
  {
    weights: weights.copy(),
    strata: strata.copy(),
    clusters: clusters.copy(),
    population_size: population_size.max(weights.length()),
    sample_size: weights.length(),
    with_replacement,
  }
}

///|
/// Computes Horvitz-Thompson total.
pub fn horvitz_thompson_total(
  values : Array[Double],
  inclusion_probabilities : Array[Double],
) -> Double {
  let n = values.length().min(inclusion_probabilities.length())
  let mut total = 0.0
  for i in 0..<n {
    total += values[i] / safe_probability(inclusion_probabilities[i])
  }
  total
}

///|
/// Computes Horvitz-Thompson mean under a supplied population size.
pub fn horvitz_thompson_mean(
  values : Array[Double],
  inclusion_probabilities : Array[Double],
  population_size : Int,
) -> Double {
  if population_size <= 0 {
    return 0.0
  }
  horvitz_thompson_total(values, inclusion_probabilities) /
  population_size.to_double()
}

///|
/// Computes a Hájek normalized weighted mean.
pub fn hajek_mean(values : Array[Double], weights : Array[Double]) -> Double {
  weighted_mean(values, weights)
}

///|
/// Computes a Hájek weighted treatment effect.
pub fn survey_weighted_effect(
  outcomes : Array[Double],
  treatment : Array[Bool],
  weights : Array[Double],
  design : SurveyDesign,
) -> SurveyEffect {
  let n = outcomes.length().min(treatment.length()).min(weights.length())
  let treated_values = Array::new()
  let control_values = Array::new()
  let treated_weights = Array::new()
  let control_weights = Array::new()
  for i in 0..<n {
    if treatment[i] {
      treated_values.push(outcomes[i])
      treated_weights.push(weights[i])
    } else {
      control_values.push(outcomes[i])
      control_weights.push(weights[i])
    }
  }
  let estimate = weighted_mean(treated_values, treated_weights) -
    weighted_mean(control_values, control_weights)
  let treated_ess = effective_sample_size(treated_weights)
  let control_ess = effective_sample_size(control_weights)
  let treated_var = weighted_variance(treated_values, treated_weights)
  let control_var = weighted_variance(control_values, control_weights)
  let treated_component = if treated_ess == 0.0 {
    0.0
  } else {
    treated_var / treated_ess
  }
  let control_component = if control_ess == 0.0 {
    0.0
  } else {
    control_var / control_ess
  }
  let standard_error = (treated_component + control_component).sqrt()
  let ess = treated_ess + control_ess
  let deff = if ess == 0.0 { 0.0 } else { n.to_double() / ess }
  let fpc = if design.population_size <= 0 || design.population_size <= n {
    1.0
  } else {
    ((design.population_size - n).to_double() /
    (design.population_size - 1).to_double()).sqrt()
  }
  {
    estimate,
    standard_error: standard_error * fpc,
    effective_sample_size: ess,
    design_effect: deff,
    finite_population_correction: fpc,
    passes: treated_ess > 1.0 && control_ess > 1.0,
  }
}

///|
/// Calculates a post-stratified mean using population cell counts.
pub fn poststratified_mean(
  keys : Array[Int],
  values : Array[Double],
  weights : Array[Int],
) -> (Double, Array[PostStratum]) {
  let n = keys.length().min(values.length())
  let groups : Array[Int] = Array::new()
  for key in keys[:n] {
    if !groups.contains(key) {
      groups.push(key)
    }
  }
  let summaries : Array[PostStratum] = Array::new(capacity=groups.length())
  let mut result = 0.0
  let mut population_total = 0
  for key in groups {
    let selected = Array::new()
    for i in 0..<n {
      if keys[i] == key {
        selected.push(values[i])
      }
    }
    let population = if key >= 0 && key < weights.length() {
      weights[key]
    } else {
      selected.length()
    }
    summaries.push({
      key,
      sample_count: selected.length(),
      population_count: population,
      sample_mean: mean_or(selected, 0.0),
      weighted_contribution: 0.0,
    })
    population_total += population
  }
  for summary in summaries {
    let contribution = if population_total == 0 {
      0.0
    } else {
      summary.sample_mean *
      summary.population_count.to_double() /
      population_total.to_double()
    }
    result += contribution
  }
  (result, summaries)
}

///|
/// Adjusts weights toward a target total using a bounded ratio.
pub fn calibrate_weight_totals(
  weights : Array[Double],
  target_total : Double,
  maximum_ratio? : Double = 10.0,
) -> CalibrationResult {
  let current = sum(weights)
  if current == 0.0 {
    return {
      weights: Array::make(weights.length(), 0.0),
      iterations: 0,
      maximum_margin_error: target_total.abs(),
      effective_sample_size: 0.0,
      converged: false,
    }
  }
  let ratio = clamp(target_total / current, 0.0, maximum_ratio)
  let adjusted = Array::new(capacity=weights.length())
  for weight in weights {
    adjusted.push(weight * ratio)
  }
  {
    weights: adjusted,
    iterations: 1,
    maximum_margin_error: (sum(adjusted) - target_total).abs(),
    effective_sample_size: effective_sample_size(adjusted),
    converged: true,
  }
}

///|
/// Performs iterative proportional fitting on integer cells.
pub fn rake_integer_cells(
  keys : Array[Int],
  weights : Array[Double],
  target_totals : Array[Double],
  iterations : Int,
) -> CalibrationResult {
  let result = weights.copy()
  let groups : Array[Int] = Array::new()
  for key in keys {
    if !groups.contains(key) {
      groups.push(key)
    }
  }
  let rounds = if iterations > 0 { iterations } else { 1 }
  let mut margin_error = 0.0
  for iteration in 0..<rounds {
    margin_error = 0.0
    for key in groups {
      let mut current = 0.0
      for i in 0..<keys.length().min(result.length()) {
        if keys[i] == key {
          current += result[i]
        }
      }
      let target = if key >= 0 && key < target_totals.length() {
        target_totals[key]
      } else {
        current
      }
      let ratio = if current == 0.0 { 1.0 } else { target / current }
      for i in 0..<keys.length().min(result.length()) {
        if keys[i] == key {
          result[i] *= ratio
        }
      }
      let error = (current - target).abs()
      if error > margin_error {
        margin_error = error
      }
    }
    if margin_error < 1.0e-8 {
      return {
        weights: result,
        iterations: iteration + 1,
        maximum_margin_error: margin_error,
        effective_sample_size: effective_sample_size(result),
        converged: true,
      }
    }
  }
  {
    weights: result,
    iterations: rounds,
    maximum_margin_error: margin_error,
    effective_sample_size: effective_sample_size(result),
    converged: margin_error < 1.0e-6,
  }
}

///|
/// Computes a weighted quantile without relying on sampling assumptions.
pub fn survey_weighted_quantile(
  values : Array[Double],
  weights : Array[Double],
  probability : Double,
) -> Double {
  let n = values.length().min(weights.length())
  if n == 0 {
    return 0.0
  }
  let order : Array[Int] = Array::new(capacity=n)
  for i in 0..<n {
    order.push(i)
  }
  for i in 1..<n {
    let value = order[i]
    let mut cursor = i
    while cursor > 0 && values[order[cursor - 1]] > values[value] {
      order[cursor] = order[cursor - 1]
      cursor -= 1
    }
    order[cursor] = value
  }
  let total = sum(weights[:n].to_owned())
  let target = clamp(probability, 0.0, 1.0) * total
  let mut cumulative = 0.0
  for index in order {
    cumulative += weights[index]
    if cumulative >= target {
      return values[index]
    }
  }
  values[order[n - 1]]
}

///|
/// Computes a weighted variance with a finite-population adjustment.
pub fn survey_weighted_variance(
  values : Array[Double],
  weights : Array[Double],
  population_size : Int,
) -> Double {
  let base = weighted_variance(values, weights)
  if population_size <= values.length() || population_size <= 1 {
    base
  } else {
    base * population_size.to_double() / (population_size - 1).to_double()
  }
}

///|
/// Computes a cluster-level survey design effect.
pub fn survey_cluster_design_effect(
  clusters : Array[Int],
  weights : Array[Double],
) -> Double {
  let groups : Array[Int] = Array::new()
  for cluster in clusters {
    if !groups.contains(cluster) {
      groups.push(cluster)
    }
  }
  if groups.length() == 0 {
    return 0.0
  }
  let sizes = Array::new()
  for cluster in groups {
    let mut total = 0.0
    for i in 0..<clusters.length().min(weights.length()) {
      if clusters[i] == cluster {
        total += weights[i]
      }
    }
    sizes.push(total)
  }
  if sizes.length() == 0 {
    0.0
  } else {
    sum(sizes) / effective_sample_size(sizes).max(1.0)
  }
}

///|
/// Generates delete-one-cluster replicate weights.
pub fn delete_one_cluster_weights(
  clusters : Array[Int],
  weights : Array[Double],
) -> Array[Array[Double]] {
  let groups : Array[Int] = Array::new()
  for cluster in clusters {
    if !groups.contains(cluster) {
      groups.push(cluster)
    }
  }
  let result : Array[Array[Double]] = Array::new(capacity=groups.length())
  for removed in groups {
    let replicate = Array::new(capacity=weights.length())
    for i in 0..<weights.length() {
      replicate.push(
        if i < clusters.length() && clusters[i] == removed {
          0.0
        } else {
          weights[i]
        },
      )
    }
    result.push(replicate)
  }
  result
}

///|
/// Returns a compact survey summary vector.
pub fn survey_summary(effect : SurveyEffect) -> Array[Double] {
  [
    effect.estimate,
    effect.standard_error,
    effect.effective_sample_size,
    effect.design_effect,
    effect.finite_population_correction,
    if effect.passes {
      1.0
    } else {
      0.0
    },
  ]
}
