///|
/// Calculates a conservative Gamma multiplier for the legacy p-value API.
/// Use rosenbaum_probability_bounds and sensitivity_curve for the full envelope.
pub fn rosenbaum_bounds_gamma(p_value : Double, gamma : Double) -> Double {
  // A conservative multiplicative adjustment keeps the compatibility API bounded.
  let adjusted = p_value * gamma
  if adjusted > 1.0 {
    1.0
  } else {
    adjusted
  }
}

///|
pub struct RosenbaumBounds {
  gamma : Double
  lower_probability : Double
  upper_probability : Double
  amplification : Double
}

///|
/// Computes the odds-ratio sensitivity envelope for an observed probability.
pub fn rosenbaum_probability_bounds(
  p_value : Double,
  gamma : Double,
) -> RosenbaumBounds {
  let probability = clamp(p_value, 0.0, 1.0)
  let g = if gamma < 1.0 { 1.0 } else { gamma }
  let odds = if probability == 1.0 {
    1.0e300
  } else {
    probability / (1.0 - probability)
  }
  let lower_odds = odds / g
  let upper_odds = odds * g
  let lower = if lower_odds == 1.0e300 {
    1.0
  } else {
    lower_odds / (1.0 + lower_odds)
  }
  let upper = if upper_odds == 1.0e300 {
    1.0
  } else {
    upper_odds / (1.0 + upper_odds)
  }
  {
    gamma: g,
    lower_probability: lower,
    upper_probability: upper,
    amplification: upper - lower,
  }
}

///|
pub fn sensitivity_curve(
  p_value : Double,
  maximum_gamma : Double,
  steps : Int,
) -> Array[RosenbaumBounds] {
  let count = if steps < 1 { 1 } else { steps }
  let result : Array[RosenbaumBounds] = Array::new(capacity=count)
  for i in 0..<count {
    let fraction = if count == 1 {
      1.0
    } else {
      i.to_double() / (count - 1).to_double()
    }
    let max_gamma = if maximum_gamma < 1.0 { 1.0 } else { maximum_gamma }
    let gamma = 1.0 + fraction * (max_gamma - 1.0)
    result.push(rosenbaum_probability_bounds(p_value, gamma))
  }
  result
}

///|
/// Returns the minimum confounding strength required to move a p-value above alpha.
pub fn critical_gamma(
  p_value : Double,
  alpha : Double,
  maximum_gamma : Double,
) -> Double {
  let target = clamp(alpha, 0.0, 1.0)
  if rosenbaum_bounds_gamma(p_value, 1.0) > target {
    return 1.0
  }
  let mut low = 1.0
  let mut high = if maximum_gamma < 1.0 { 1.0 } else { maximum_gamma }
  for _ in 0..<60 {
    let middle = (low + high) / 2.0
    if rosenbaum_bounds_gamma(p_value, middle) > target {
      high = middle
    } else {
      low = middle
    }
  }
  high
}

///|
/// E-value for a positive risk ratio, a compact unmeasured-confounding diagnostic.
pub fn e_value(risk_ratio : Double) -> Double {
  let ratio = if risk_ratio < 1.0 { 1.0 / risk_ratio } else { risk_ratio }
  if ratio <= 1.0 {
    1.0
  } else {
    ratio + (ratio * (ratio - 1.0)).sqrt()
  }
}
