///|
/// Assigns each observation to a stratum based on uniformly spaced bins of propensity scores.
pub fn stratify_by_score(propensity_scores : Array[Double], num_strata : Int) -> Array[Int] {
  let n = propensity_scores.length()
  if n == 0 || num_strata <= 1 {
    return Array::make(n, 0)
  }
  let mut min_p = propensity_scores[0]
  let mut max_p = propensity_scores[0]
  for p in propensity_scores {
    if p < min_p {
      min_p = p
    }
    if p > max_p {
      max_p = p
    }
  }
  let range = max_p - min_p
  let strata = Array::new(capacity=n)
  for p in propensity_scores {
    if range == 0.0 {
      strata.push(0)
    } else {
      let mut bin = ((p - min_p) / range * num_strata.to_double()).to_int()
      if bin == num_strata {
        bin -= 1
      }
      strata.push(bin)
    }
  }
  strata
}
