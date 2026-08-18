///|
pub struct StratifiedEffect {
  stratum : Int
  count : Int
  treated_count : Int
  control_count : Int
  estimate : Double
  standard_error : Double
}

///|
pub fn stratified_effects(
  outcomes : Array[Double],
  treatment : Array[Bool],
  strata : Array[Int],
) -> Array[StratifiedEffect] {
  let levels = Array::new()
  for stratum in strata {
    if !levels.contains(stratum) {
      levels.push(stratum)
    }
  }
  let result : Array[StratifiedEffect] = Array::new(capacity=levels.length())
  for level in levels {
    let y = Array::new()
    let t = Array::new()
    for i in 0..<strata.length() {
      if strata[i] == level {
        y.push(outcomes[i])
        t.push(treatment[i])
      }
    }
    let treated = group_values(y, t, true)
    let control = group_values(y, t, false)
    let estimate = mean(treated) - mean(control)
    let standard_error = if treated.length() == 0 || control.length() == 0 {
      0.0
    } else {
      (variance(treated) / treated.length().to_double() +
      variance(control) / control.length().to_double()).sqrt()
    }
    result.push({
      stratum: level,
      count: y.length(),
      treated_count: treated.length(),
      control_count: control.length(),
      estimate,
      standard_error,
    })
  }
  result
}

///|
pub fn aggregate_stratified_effects(
  effects : Array[StratifiedEffect],
) -> Estimate {
  let mut weighted_sum = 0.0
  let mut total = 0
  let mut variance_sum = 0.0
  for effect in effects {
    let weight = effect.count.to_double()
    weighted_sum += weight * effect.estimate
    total += effect.count
    variance_sum += weight *
      weight *
      effect.standard_error *
      effect.standard_error
  }
  let estimate = if total == 0 { 0.0 } else { weighted_sum / total.to_double() }
  let standard_error = if total == 0 {
    0.0
  } else {
    variance_sum.sqrt() / total.to_double()
  }
  Estimate::from_standard_error(
    estimate,
    standard_error,
    total,
    total.to_double(),
    "ATE (stratified)",
  )
}

///|
pub fn overlap_population_weights(
  propensity_scores : Array[Double],
) -> Array[Double] {
  let result = Array::new(capacity=propensity_scores.length())
  for score in propensity_scores {
    let p = safe_probability(score)
    result.push(p * (1.0 - p))
  }
  result
}

///|
pub fn entropy_like_weights(
  propensity_scores : Array[Double],
  treatment : Array[Bool],
) -> Array[Double] {
  let result = Array::new(capacity=propensity_scores.length())
  for i in 0..<propensity_scores.length() {
    let p = safe_probability(propensity_scores[i])
    let target = if treatment[i] { p } else { 1.0 - p }
    result.push(if target == 0.0 { 0.0 } else { -@math.ln(target) })
  }
  result
}

///|
pub fn truncate_weights(
  weights : Array[Double],
  maximum : Double,
) -> Array[Double] {
  let bound = if maximum < 0.0 { 0.0 } else { maximum }
  let result = Array::new(capacity=weights.length())
  for weight in weights {
    result.push(if weight > bound { bound } else { weight })
  }
  result
}

///|
pub fn normalized_weights(weights : Array[Double]) -> Array[Double] {
  let total = sum(weights)
  let result = Array::new(capacity=weights.length())
  for weight in weights {
    result.push(
      if total == 0.0 {
        0.0
      } else {
        weight * weights.length().to_double() / total
      },
    )
  }
  result
}

///|
pub fn weighted_treatment_effect(
  outcomes : Array[Double],
  treatment : Array[Bool],
  weights : Array[Double],
) -> Double {
  weighted_potential_mean(outcomes, treatment, weights, true) -
  weighted_potential_mean(outcomes, treatment, weights, false)
}
