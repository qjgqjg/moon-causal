///|
/// Calculates Inverse Probability Weights for Average Treatment Effect (ATE).
pub fn ipw_weights_ate(treatments : Array[Bool], propensity_scores : Array[Double]) -> Array[Double] {
  let n = treatments.length()
  let weights = Array::new(capacity=n)
  for i in 0..<n {
    let t = treatments[i]
    let p = propensity_scores[i]
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
pub fn ipw_weights_att(treatments : Array[Bool], propensity_scores : Array[Double]) -> Array[Double] {
  let n = treatments.length()
  let weights = Array::new(capacity=n)
  for i in 0..<n {
    let t = treatments[i]
    let p = propensity_scores[i]
    if t {
      weights.push(1.0)
    } else {
      weights.push(p / (1.0 - p))
    }
  }
  weights
}
