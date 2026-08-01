///|
/// Applies the logistic sigmoid function to a value.
pub fn sigmoid(x : Double) -> Double {
  1.0 / (1.0 + @math.exp(-x))
}

///|
/// Calculates propensity scores using a linear combination of covariates and a sigmoid link.
/// In a real application, weights would be trained via Logistic Regression.
pub fn estimate_propensity_scores(covariates : Array[Array[Double]], weights : Array[Double]) -> Array[Double] {
  let n = covariates.length()
  let scores = Array::new(capacity=n)
  for i in 0..<n {
    let row = covariates[i]
    let mut z = 0.0
    for j in 0..<row.length() {
      if j < weights.length() {
        z += row[j] * weights[j]
      }
    }
    scores.push(sigmoid(z))
  }
  scores
}
