///|
/// Calculates Inverse Probability Weights for Average Treatment Effect (ATE).
pub fn ipw_weights_ate(
  treatments : Array[Bool],
  propensity_scores : Array[Double],
) -> Array[Double] {
  let n = if treatments.length() < propensity_scores.length() {
    treatments.length()
  } else {
    propensity_scores.length()
  }
  let weights = Array::new(capacity=n)
  for i in 0..<n {
    let t = treatments[i]
    let p = safe_probability(propensity_scores[i])
    if t {
      weights.push(1.0 / p)
    } else {
      weights.push(1.0 / (1.0 - p))
    }
  }
  weights
}

///|
/// Calculates Inverse Probability Weights for Average Treatment Effect on the Treated (ATT).
pub fn ipw_weights_att(
  treatments : Array[Bool],
  propensity_scores : Array[Double],
) -> Array[Double] {
  let n = if treatments.length() < propensity_scores.length() {
    treatments.length()
  } else {
    propensity_scores.length()
  }
  let weights = Array::new(capacity=n)
  for i in 0..<n {
    let t = treatments[i]
    let p = safe_probability(propensity_scores[i])
    if t {
      weights.push(1.0)
    } else {
      weights.push(p / (1.0 - p))
    }
  }
  weights
}
