///|
/// Calculates the Doubly Robust estimator for the Average Treatment Effect (ATE).
pub fn doubly_robust_ate(
  outcomes : Array[Double],
  treatments : Array[Bool],
  propensity_scores : Array[Double],
  predicted_outcomes_treated : Array[Double],
  predicted_outcomes_control : Array[Double]
) -> Double {
  let n = outcomes.length()
  if n == 0 {
    return 0.0
  }
  let mut mu_1_sum = 0.0
  let mut mu_0_sum = 0.0
  for i in 0..<n {
    let y = outcomes[i]
    let t = if treatments[i] { 1.0 } else { 0.0 }
    let p = propensity_scores[i]
    let y1 = predicted_outcomes_treated[i]
    let y0 = predicted_outcomes_control[i]

    let term1 = (t * y) / p - ((t - p) / p) * y1
    let term0 = ((1.0 - t) * y) / (1.0 - p) + ((t - p) / (1.0 - p)) * y0
    
    mu_1_sum += term1
    mu_0_sum += term0
  }
  (mu_1_sum - mu_0_sum) / n.to_double()
}
